#include <iostream>
#include <fstream>
#include <sys/mman.h>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <Eigen/Dense>
#include "plink.hpp"
#define PLINK_OFFSET 3

int main () {
  int nindiv = 3925; // number of individuals as per test.fam
  int nsnps = 1000; // number of snps as per test.bim
  Eigen::MatrixXd X(nindiv, nsnps); // Genotype
  Eigen::MatrixXd A(nindiv, nindiv); // GRM
  Eigen::MatrixXd NM(nindiv, nindiv); // Number of non-missing SNPs per each individual pair

  struct stat sb;
  int fd = -1; // file descriptor
  char* data = NULL;

  fd = open("../data/test.bed", O_RDONLY);
  fstat(fd, &sb);
  data = (char*)mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

  read_bed(data + PLINK_OFFSET, X);
  
  munmap(data, sb.st_size);
  close(fd);
  
  // For testing only!
  // std::ofstream file_geno("geno.txt", std::ios::out | std::ios::trunc);
  // file_geno << X;
  // file_geno.close();
  
  calculate_grm(X, A, NM);
  
  // For testing only!
  // std::ofstream file_grm("grm.txt", std::ios::out | std::ios::trunc);
  // file_grm << A;
  // file_grm.close();
  
  // For testing only!
  // std::ofstream file_nm("non_missing.txt", std::ios::out | std::ios::trunc);
  // file_nm << NM;
  // file_nm.close();

  off_t size = sizeof(float) * nindiv * (nindiv + 1) / 2;
  // std::cout << size << std::endl;
  int grm_file = open("a.grm.bin", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  int nm_file = open("a.grm.N.bin", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  
  int result = lseek(grm_file, size - 1, SEEK_SET);
  result = write(grm_file, "", 1);

  int result2 = lseek(nm_file, size - 1, SEEK_SET);
  result2 = write(nm_file, "", 1);


  float* grm = (float*)mmap((caddr_t)0, size, PROT_READ | PROT_WRITE, MAP_SHARED, grm_file, 0);
  float* non_missing = (float*)mmap((caddr_t)0, size, PROT_READ | PROT_WRITE, MAP_SHARED, nm_file, 0);

  int sum_up_toi = 0;
  for(int i = 0; i < nindiv; i++) {
    sum_up_toi += i;
      for(int j = 0; j <= i; j++) {
        grm[sum_up_toi + j] = (float)A(i,j);
        non_missing[sum_up_toi + j] = (float)NM(i,j);
    }
  }
  
  munmap(grm, size);
  close(grm_file);

  return 0;
}





















