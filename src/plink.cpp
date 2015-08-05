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
         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         for(unsigned int i = 0 ; i < N ; i++)
            tmp3(i) = (double)tmp2[i];
         X.col(j) = tmp3;
      }
   
   delete[] tmp;
   delete[] tmp2;
}
