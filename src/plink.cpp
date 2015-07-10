#include <iostream>
#include <math.h> 
#include "plink.hpp"

#define PACK_DENSITY 4
/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

void decode_plink(unsigned char *out,
      const unsigned char *in, const unsigned int n)
{
   unsigned int i, k;
   unsigned char tmp, geno;
   unsigned int a1, a2;

   for(i = 0 ; i < n ; ++i)
   {
      tmp = in[i];
      k = PACK_DENSITY * i;
      
      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */
      geno = (tmp & MASK0);
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : (a1 + a2);
      k++;

      geno = (tmp & MASK1) >> 2; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : (a1 + a2);
      k++;

      geno = (tmp & MASK2) >> 4; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : (a1 + a2);
      k++;

      geno = (tmp & MASK3) >> 6; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : (a1 + a2);
   }
}

void read_bed(char* data, Eigen::MatrixXd& X){
   unsigned int np;
   int N = X.rows();
   int nsnps = X.cols();
   Eigen::VectorXd tmp3(N);

   np = (unsigned int)ceil((double)N / PACK_DENSITY);
   unsigned char* tmp = new unsigned char[np];
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];

   for(unsigned int j = 0 ; j < nsnps ; j++)
      {
         // read raw genotypes
         std::copy(data + (sizeof(char) * np) * j,
                   data + (sizeof(char) * np) * (j + 1), tmp);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         for(unsigned int i = 0 ; i < N ; i++)
            tmp3(i) = (double)tmp2[i];
         X.col(j) = tmp3;
      }
   
   delete[] tmp;
   delete[] tmp2;
}

double swap_na(double x, double NA, double tol) {
         double out = (fabs(x - NA)) < tol ? 0.0 : x;
         return out;
} 

void calculate_grm(Eigen::MatrixXd& X, Eigen::MatrixXd& A, Eigen::MatrixXd& NM) {
 
   int N = X.rows();
   int nsnps = X.cols();
 
   // standardized genotype matrix
   Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(N, nsnps);
      
   count_non_missing(X, NM);
   scale_and_center_genotype(X, Z);
   A = Z * Z.transpose();
}

void update_grm(Eigen::MatrixXd& A, Eigen::MatrixXd& NM,
               Eigen::MatrixXd& Ai, Eigen::MatrixXd& NMi) {
   A  = A + Ai;
   NM = NM + NMi;
}

void scale_and_center_genotype(Eigen::MatrixXd& X, Eigen::MatrixXd& Z) {
   
   int N = X.rows();
   int nsnps = X.cols();

   // columns are genotype classes 0, 1, 2, 3(NA), frequency of the reference allele
   Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(nsnps, 5);
   
   // count number of individuals per genotype class
   for(int i = 0; i < nsnps; i++) {
      Eigen::VectorXd i_snp = X.col(i);
      for(int j = 0; j < N; j++) {
         Y(i, i_snp[j])++;
      }
   }

   // Calculate reference allele frequency
   for(int i = 0; i < nsnps; i++) {
      
      double RAF, MAF;
      RAF = (0.5 * Y(i, 1) + Y(i, 2)) / (Y(i, 0) + Y(i, 1) + Y(i, 2));
      
      if(RAF > 0.5) {
         MAF = 1.0 - RAF;
        
        // swap geno
         std::for_each(X.col(i).data(), X.col(i).data() + N, [](double &geno) {
            if(geno != 3)
               geno = 2 - geno; // otherwise do nothing and keep NA as 3
         });

      } else {
         MAF = RAF;
      }
      
      Y(i, 4) = MAF;
   }

   // Calculate standardized genotype matrix
   for(int i = 0; i < nsnps; i++) {
      Eigen::VectorXd i_snp = X.col(i);
      
      i_snp = (i_snp.array() - 2 * Y(i, 4)) / sqrt(2 * Y(i, 4) * (1 - Y(i, 4)));
      
      // NA value after centering and scaling, (before it was 3)
      double NA = (3 - 2 * Y(i, 4)) / sqrt(2 * Y(i, 4) * (1 - Y(i, 4)));
     
      // set NAs to zero
      std::for_each(i_snp.data(), i_snp.data() + N, [&](double &d) {
         d = swap_na(d, NA, 0.0000001);
      });
     
      // std::cout << i_snp.transpose() << std::endl;
      Z.col(i) = i_snp;
   }

}

void count_non_missing(Eigen::MatrixXd& X, Eigen::MatrixXd& NM) {
   int N = X.rows();
   int nsnps = X.cols();
   
   Eigen::MatrixXd NMG = Eigen::MatrixXd::Zero(N, nsnps);
   
   // missing value is zero, one otherwise.
   for(int i = 0; i < N; i++) {
      for(int j = 0; j < nsnps; j++) {
         if(X(i,j) != 3.0) {
            NMG(i,j) = 1.0;
         } else {
            NMG(i,j) = 0.0;
         }
      }
   }

   NM = NMG * NMG.transpose();
}











