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
  
  // std::cout << X << std::endl;

  calculate_grm(X, A);
  // write_grm();
  // update_grm();
  std::cout << A << std::endl;
  
  return 0;
}
