#define PLINK_OFFSET 3 // 3-bytes offset in plink binary format 

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
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK1) >> 2; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK2) >> 4; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK3) >> 6; 
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
   }
}


void Data::read_bed(int impute)
{
   
   in.seekg(0, std::ifstream::end);

   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
   len = (unsigned int)in.tellg() - 3;

   // size of packed data, in bytes, per SNP
   np = (unsigned int)ceil((double)N / PACK_DENSITY);
   nsnps = len / np;
   in.seekg(3, std::ifstream::beg);

   unsigned char* tmp = new unsigned char[np];

   // Allocate more than the sample size since data must take up whole bytes
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
   X = MatrixXd(N, nsnps);

   VectorXd tmp3(N);

   for(unsigned int j = 0 ; j < nsnps ; j++)
      {
         // read raw genotypes
         in.read((char*)tmp, sizeof(char) * np);

         // decode the genotypes
         decode_plink(tmp2, tmp, np);

         for(unsigned int i = 0 ; i < N ; i++)
            tmp3(i) = (double)tmp2[i];
         X.col(j) = tmp3;
      }
   
   X.col(j) = tmp3;
   p = X.cols();

   delete[] tmp;
   delete[] tmp2;
   

   in.close();
}