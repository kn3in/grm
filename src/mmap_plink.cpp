#include <iostream>
#include <sys/mman.h>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <Eigen/Dense>
#include "plink.hpp"


int main () {
  int nindiv = 3925; // number of individuals as per test.fam
  int nsnps = 1000; // number of snps as per test.bim
  Eigen::MatrixXd X(nindiv, nsnps);
  Eigen::MatrixXd A(nindiv, nindiv);
  
  struct stat sb;
  int fd = -1; // file descriptor
  char* data = NULL;

  fd = open("../data/test.bed", O_RDONLY);
  fstat(fd, &sb);
  data = (char*)mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

  read_bed(data, X);
  
  munmap(data, sb.st_size);
  close(fd);
  
  calculate_grm(X, A);
  
  off_t size = sizeof(int) * nindiv * (nindiv + 1) / 2;
  std::cout << size << std::endl;
  int grm_file = open("a.grm.bin", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  int result = lseek(grm_file, size - 1, SEEK_SET);
  result = write(grm_file, "", 1);

  float* grm = (float*)mmap((caddr_t)0, size, PROT_READ | PROT_WRITE, MAP_SHARED, grm_file, 0);
  
  for(int i = 0; i < nindiv; i++) {
    for(int j = 0; j <= i; j++) {
      grm[i + j] = A(i,j);
    }
  }

  munmap(grm, size);
  close(grm_file);
  
  return 0;
}
