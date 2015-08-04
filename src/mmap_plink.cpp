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
#include "betas.hpp"
#include "plink_data.hpp"

#define PLINK_OFFSET 3
#define PACK_DENSITY 4

int main (int argc, char* argv[]) {
  
  std::string usage("USAGE: bgrm plink_data betas_data out_name");
    
  if(argc != 4) {
    std::cerr << usage << std::endl;
      return 1;
  }
  
  std::string plink_base = argv[1];
  std::string bim_file = plink_base + ".bim";
  std::string fam_file = plink_base + ".fam";
  std::string bed_file = plink_base + ".bed";

  std::string my_betas_path = argv[2];
  std::string results_file = argv[3];

  std::string grm_bin = results_file + ".grm.bin";
  std::string grm_n = results_file + ".grm.N.bin";
  std::string grm_id = results_file + ".grm.id";
  std::string g_pred = results_file + "_gpredictor.txt";

  //-----------------------------------------------------------------------------
  
  betas my_betas(my_betas_path);
  bim_data bim(bim_file);
  fam_data fam(fam_file);

  // which snps present only in the bed and not in the betas file?
  bim.setup_snps_without_betas(my_betas);
  // set up <bool> vector to decide which SNPs are in common between bed adn betas files
  bim.setup_snps_to_iterate();
  // itereate over SNPs common between bed/betas files and swap sign of the effect size if reference alleles differ 
  my_betas.swap_betas(bim);




  //-----------------------------------------------------------------------------
  
  int nindiv = fam.nindiv; // number of individuals as per test.fam
  int nsnps = bim.nsnps; // number of snps as per test.bim
  // Eigen::MatrixXd X(nindiv, nsnps); // Genotype
  Eigen::MatrixXd A(nindiv, nindiv); // GRM
  Eigen::MatrixXd NM(nindiv, nindiv); // Number of non-missing SNPs per each individual pair
  Eigen::VectorXd g_hat(nindiv); // genetic effects aka row sums

  A.setZero();
  NM.setZero();
  g_hat.setZero();


  std::ifstream bed;
  bed.open (bed_file, std::ios::in | std::ios::binary); 
  bed.seekg(3, std::ifstream::beg);
  
  // struct stat sb;
  // int fd = -1; // file descriptor
  // char* data = NULL;

  // fd = open(bed_file.c_str(), O_RDONLY);
  // fstat(fd, &sb);
  // data = (char*)mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

  int chunk_size = nsnps; // SNPs to read in RAM at a time
  std::vector<int> start_index_of_chunk;
  
  int chunks = ceil( (double)nsnps / chunk_size);
  int last_chunk = nsnps % chunk_size; 
  
  std::vector<int> snps_per_chunk(chunks, chunk_size);
  if(last_chunk) snps_per_chunk.back() = last_chunk; // Fix it with if()

  for(size_t i = 0; i < snps_per_chunk.size(); ++i) {
    int n_snps = snps_per_chunk.front(); // last chunk start is still multiple of big chunks!!! 
    unsigned int np;
    np = (unsigned int)ceil((double)nindiv / PACK_DENSITY);
    unsigned int start_index = PLINK_OFFSET + (sizeof(char) * np) * n_snps * i;
    start_index_of_chunk.push_back(start_index);
  }

  for(size_t i = 0; i < snps_per_chunk.size(); ++i) {
    Eigen::MatrixXd Xi(nindiv, snps_per_chunk[i]); //current Genotype
    Eigen::MatrixXd NMGi(nindiv, snps_per_chunk[i]); //current non-missing->1 NA->0
    Eigen::MatrixXd Ai(nindiv, nindiv); // current GRM
    Eigen::MatrixXd NMi(nindiv, nindiv); // current Number of non-missing SNPs

    Eigen::VectorXd betas(snps_per_chunk[i]); // current betas    
    betas.setZero();

    read_bed2(bed, Xi);
    
    //zeroing betas of non-overlapping SNPs
    
    int down = i * snps_per_chunk.front();
    int up = down + snps_per_chunk[i];
    
    for(int j = down; j < up; j++) {
      int local_i = j - down;
      if(bim.snps_to_use[j]) {
        int betas_ind = my_betas.rs_id2index[ bim.bim_rs_id[j] ];
        betas[local_i] = my_betas.effects[betas_ind];
      } else {
        // genotypes of all SNPs which do NOT have betas set to 3 (NA value)
        Xi.col(local_i) = Eigen::VectorXd::Constant(nindiv, 3);
      }    
    }

    count_non_missing(Xi, NMi, NMGi); // NMGi keeps track of missing values!
    swap_na_matrix(Xi); // all NA turns into zero
 
    // change Xi to genotype * effect_size
    for(int l = 0; l < Xi.cols(); l++) {
      Xi.col(l) = Xi.col(l) * betas[l]; 
    }

    // effects per individual
    g_hat += Xi.rowwise().sum(); 

    calculate_grm2(Xi, Ai, NMGi);
    update_grm(A, NM, Ai, NMi);
  }
  
  // munmap(data, sb.st_size);
  // close(fd);
  bed.close();

  // scale by number of non-missing genotypes 
  for(int i = 0; i < nindiv; i++) {
     for(int j = 0; j < nindiv; j++) {
        A(i,j) /= NM(i,j) - 1;
     }
  }

  off_t size = sizeof(float) * nindiv * (nindiv + 1) / 2;

  int grm_file = open(grm_bin.c_str(), O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  int nm_file = open(grm_n.c_str(), O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  
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

  std::ofstream id_grm(grm_id);
  std::ofstream g_eff(g_pred);


  for(size_t i = 0; i < fam.individual_id.size(); i++) {
    id_grm << i + 1 << "\t" << fam.individual_id[i] << "\n";
    g_eff << fam.family_id[i] << "\t" << fam.individual_id[i] << "\t" << g_hat[i] << "\n";
  }

  return 0;
}





















