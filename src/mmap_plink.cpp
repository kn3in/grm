#include <iostream>
#include <fstream>
#include <numeric>
#include <typeinfo>
#include <sys/mman.h>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <Eigen/Dense>
#include "plink.hpp"
#define PLINK_OFFSET 3
#define PACK_DENSITY 4

int main () {
  int nindiv = 3925; // number of individuals as per test.fam
  int nsnps = 1000; // number of snps as per test.bim
  // Eigen::MatrixXd X(nindiv, nsnps); // Genotype
  // Eigen::MatrixXd A(nindiv, nindiv); // GRM
  // Eigen::MatrixXd NM(nindiv, nindiv); // Number of non-missing SNPs per each individual pair

  struct stat sb;
  int fd = -1; // file descriptor
  char* data = NULL;

  fd = open("../data/test.bed", O_RDONLY);
  fstat(fd, &sb);
  data = (char*)mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

  int chunk_size = 89; // SNPs to read in RAM at a time
  std::vector<int> start_index_of_chunk;
  
  int chunks = ceil( (double)nsnps / chunk_size);
  int last_chunk = nsnps % chunk_size; 
  
  std::vector<int> snps_per_chunk(chunks, chunk_size);
  if(last_chunk) snps_per_chunk.back() = last_chunk;

  for(int i = 0; i < snps_per_chunk.size(); ++i) {
    int n_snps = snps_per_chunk.front(); // last chunk start is still multiple of big chunks!!! 
    unsigned int np;
    np = (unsigned int)ceil((double)nindiv / PACK_DENSITY);
    unsigned int start_index = PLINK_OFFSET + (sizeof(char) * np) * n_snps * i;
    start_index_of_chunk.push_back(start_index);
  }
  
  off_t size = sizeof(double) * nindiv * (nindiv + 1) / 2;
  // std::cout << size << std::endl;
  int grm_file = open("a.grm.bin", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  int nm_file = open("a.grm.N.bin", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  
  int result = lseek(grm_file, size - 1, SEEK_SET);
  result = write(grm_file, "", 1);

  int result2 = lseek(nm_file, size - 1, SEEK_SET);
  result2 = write(nm_file, "", 1);


  double* grm = (double*)mmap((caddr_t)0, size, PROT_READ | PROT_WRITE, MAP_SHARED, grm_file, 0);
  double* non_missing = (double*)mmap((caddr_t)0, size, PROT_READ | PROT_WRITE, MAP_SHARED, nm_file, 0);

  // fill mmaped files with zeros, need better way of doing this!!!
  int sum_up_toi = 0;
    for(int i = 0; i < nindiv; i++) {
      sum_up_toi += i;
        for(int j = 0; j <= i; j++) {
          grm[sum_up_toi + j] = 0.0;
          non_missing[sum_up_toi + j] = 0.0;
      }
  }
  
  // int size_low_tri = nindiv * (nindiv + 1) / 2;
  
  for(int i = 0; i < snps_per_chunk.size(); ++i) {
    Eigen::MatrixXd Xi(nindiv, snps_per_chunk[i]); //current Genotype
    read_bed(data + start_index_of_chunk[i], Xi);
    
    int N = Xi.rows();
    int nsnps = Xi.cols();

    Eigen::MatrixXd Zi = Eigen::MatrixXd::Zero(N, nsnps); // current scaled-centered genotype
    Eigen::MatrixXd NMGi = Eigen::MatrixXd::Zero(N, nsnps); // current non-missing
    
    count_non_missing3(Xi, NMGi);
    scale_and_center_genotype(Xi, Zi);

    Eigen::MatrixXd Zi_t = Zi.transpose();
    Eigen::MatrixXd NMGi_t = NMGi.transpose();
    
    sum_up_toi = 0;
    // #pragma omp parallel for
    for(int i = 0; i < N; i++) {
      sum_up_toi += i;
      #pragma omp parallel for
      for(int j = 0; j <= i; j++) {
         
          grm[sum_up_toi + j] = grm[sum_up_toi + j] + (Zi_t.col(i)).dot(Zi_t.col(j));
          non_missing[sum_up_toi + j] = non_missing[sum_up_toi + j] + (NMGi_t.col(i)).dot(NMGi_t.col(j));
      
      }
   }

  }
  
  munmap(data, sb.st_size);
  close(fd);

  // scale by the number of non-missing but utilize mmap'ed arrays
  sum_up_toi = 0;
  for(int i = 0; i < nindiv; i++) {
    sum_up_toi += i;
      for(int j = 0; j <= i; j++) {
        grm[sum_up_toi + j] =  grm[sum_up_toi + j] / non_missing[sum_up_toi + j];
    }
  }
  
  munmap(grm, size);
  close(grm_file);

  return 0;
}





















