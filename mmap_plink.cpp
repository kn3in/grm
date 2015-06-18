#include <iostream>
#include <sys/mman.h>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>


using namespace std;

int main () {
  int nindiv = 3925; // number of individuals as per test.fam
  int nsnps = 1000; // number of snps as per test.bim

  //  xxd -b test.bed | head
  // 0000000: 01101100 00011011 00000001 10 10 11 11 11101111 10 01 10 11  l.....
  // a value of 00000001 indicates SNP-major
  // (i.e. list all individuals for first SNP, all individuals for second SNP, etc
  // The bed file in SNP major mode (see above)

  // size_t size = 1UL << 30;
  // cout << size << endl;
  
  struct stat sb;
  int fd = -1; // file descriptor
  int pagesize = -1;
  char* data = NULL;

  fd = open("test.bed", O_RDONLY);
  fstat(fd, &sb);
  printf("Size: %lu\n", (uint64_t)sb.st_size);
  pagesize = getpagesize();
  // data = (char*)mmap((caddr_t)0, pagesize, PROT_READ, MAP_SHARED, fd, 0);
  data = (char*)mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  // cout << "The page size is :" << pagesize << "bytes" << endl;
  
  // for(int x = 0; x < pagesize; ++x) {
  //     cerr << data[x];
  // }

  // double genotypes[nsnps][nindiv];
  
  // for(int i = 0; i < nsnps; i++) {
  //     int len ;
  //     double decoded[nindiv];
  //     snp_data = copy(data + 3 + i*nindiv, )
  //     decoded = decode_snp(snp_data, nindiv);
  //     genotypes[i][] = decoded;
  // }

  
  char ffg = data[3]; // first four genotypes;
  vector<int> ffd_decoded;

  char geno;
  int a1, a2;
  int final;
  
    geno = (ffg & MASK0);
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    final = (geno == 1) ? 3 : a1 + a2;
    ffd_decoded.push_back(final);

    geno = (ffg & MASK1) >> 2;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    final = (geno == 1) ? 3 : a1 + a2;
    ffd_decoded.push_back(final);

    geno = (ffg & MASK2) >> 4;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    final = (geno == 1) ? 3 : a1 + a2;
    ffd_decoded.push_back(final);

    geno = (ffg & MASK3) >> 6;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    final = (geno == 1) ? 3 : a1 + a2;
    ffd_decoded.push_back(final);

    //reverse(ffd_decoded.begin(), ffd_decoded.end());

  for(int i = 0; i < 4; i++) {
    cout << "The genotype of individual " << i + 1 << " is: " << ffd_decoded[i] << endl;
  }

  
  munmap(data, sb.st_size);
  close(fd);
  return 0;
}
