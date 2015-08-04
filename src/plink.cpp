#include <iostream>
#include <fstream>
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
   int unsigned N = X.rows();
   int unsigned nsnps = X.cols();
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

void read_bed2(std::ifstream& bed, Eigen::MatrixXd& X){
   unsigned int np;
   unsigned int N = X.rows();
   unsigned int nsnps = X.cols();
   Eigen::VectorXd tmp3(N);

   np = (unsigned int)ceil((double)N / PACK_DENSITY);
   unsigned char* tmp = new unsigned char[np];
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];

   for(unsigned int j = 0 ; j < nsnps ; j++)
      {
         // read raw genotypes
         bed.read((char*)tmp, sizeof(char) * np);

         // std::copy(data + (sizeof(char) * np) * j,
         //           data + (sizeof(char) * np) * (j + 1), tmp);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         for(unsigned int i = 0 ; i < N ; i++)
            tmp3(i) = (double)tmp2[i];
         X.col(j) = tmp3;
      }
   
   delete[] tmp;
   delete[] tmp2;
}

void update_grm(Eigen::MatrixXd& A, Eigen::MatrixXd& NM,
               Eigen::MatrixXd& Ai, Eigen::MatrixXd& NMi) {
   A  = A + Ai;
   NM = NM + NMi;
}

void count_non_missing(Eigen::MatrixXd& X, Eigen::MatrixXd& NM, Eigen::MatrixXd& NMG) {
   int N = X.rows();
   int nsnps = X.cols();
   
   // NMG = Eigen::MatrixXd::Zero(N, nsnps);
   
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

void swap_na_matrix(Eigen::MatrixXd& X) {
   int nindiv = X.rows();
   for(int i = 0; i < X.cols(); i++) {
      std::for_each(X.col(i).data(), X.col(i).data() + nindiv, [](double &geno) {
                  if(geno == 3)
                     geno = 0; // swap NA representation from 3 to 0
               });
   }
}

void calculate_grm2(Eigen::MatrixXd& Xi, Eigen::MatrixXd& Ai, Eigen::MatrixXd& NMG) {
   
   // Eigen::VectorXd NM_per_column = NMG.colwise().sum

   // here is an easy way of getting bugs MAKE SURE AT LEAST ONE SNP IS NOT MISSING  
   Eigen::VectorXd Mean = Xi.colwise().sum().array() / NMG.colwise().sum().array();
   
   // center
   Xi.rowwise() -= Mean.transpose();
   // NA is not zero anymore so bring it back to zero!!!
   swap_na_mm(Xi, NMG);
   
   Eigen::VectorXd Sd = Xi.colwise().norm().array() / (NMG.colwise().sum().array() - 1).sqrt();

   // scale
   for(int i = 0; i < Xi.cols(); i++) {
      if(Sd[i] != 0.0) { // This is quite a wild assumption testing for equality, rather use greater than tolerance?
         Xi.col(i) /= Sd[i];
      }
   }

   // swap is not nessesary cause 0 divided by sd is still 0
   //swap_na_mm(Xi, NMG);

   // Eigen::MatrixXd Yi = Xi.transpose();
   Eigen::VectorXd Xi_Mean = Xi.rowwise().sum().array() / NMG.rowwise().sum().array();
   Xi.colwise() -= Xi_Mean;
   swap_na_mm(Xi, NMG);
   Ai = Xi * Xi.transpose();
}

void swap_na_mm(Eigen::MatrixXd& X, Eigen::MatrixXd& NMG) {
   for(int i = 0; i < X.cols(); i++) {
      for(int j = 0; j < X.rows(); j++) {
         if(!NMG(j, i))
            X(j, i) = 0.0;
      }
   }
}

// supposse no SNPs is missing
void calculate_grm3(Eigen::MatrixXd& Xi, Eigen::MatrixXd& Ai) {
   
   // center
   Eigen::VectorXd Mean = Xi.colwise().sum() / Xi.rows();
   Xi.rowwise() -= Mean.transpose();
   
   // Scale
   Eigen::VectorXd Sd = Xi.colwise().norm() / sqrt(Xi.rows() - 1);

   for(int i = 0; i < Xi.cols(); i++) {
      if(Sd[i] != 0.0) { // This is quite a wild assumption testing for equality, rather use greater than tolerance?
         Xi.col(i) /= Sd[i];
      }
   }


   Eigen::MatrixXd Yi = Xi.transpose();
   // Now center columns after the transpose above
   Eigen::VectorXd Yi_Mean = Yi.colwise().sum() / Yi.rows();
   Yi.rowwise() -= Yi_Mean.transpose();

   Ai = Yi.transpose() * Yi;
}



































