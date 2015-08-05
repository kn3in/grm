#include <iostream>
#include <fstream>
#include <numeric>
#include <typeinfo>
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
  my_betas.order_betas(bim);




  //-----------------------------------------------------------------------------
  
  int nindiv = fam.nindiv; // number of individuals as per test.fam
  int nsnps = bim.nsnps; // number of snps as per test.bim
  
  Eigen::MatrixXd A(nindiv, nindiv); // GRM
  Eigen::MatrixXd NM(nindiv, nindiv); // Number of non-missing SNPs per each individual pair
  Eigen::VectorXd g_hat(nindiv); // genetic effects aka row sums

  A.setZero();
  NM.setZero();
  g_hat.setZero();
  
  std::ifstream bed;
  bed.open (bed_file, std::ios::in | std::ios::binary); 
  bed.seekg(3, std::ifstream::beg);
  bim.setup_snps_per_chunk(nsnps);
  
  for(size_t i = 0; i < bim.snps_per_chunk.size(); ++i) {
    Eigen::MatrixXd Xi(nindiv, bim.snps_per_chunk[i]); //current Genotype
    Eigen::MatrixXd NMGi(nindiv, bim.snps_per_chunk[i]); //current non-missing->1 NA->0
    Eigen::MatrixXd Ai(nindiv, nindiv); // current GRM
    Eigen::MatrixXd NMi(nindiv, nindiv); // current Number of non-missing SNPs

    read_bed2(bed, Xi);

    int down = i * bim.snps_per_chunk.front();
    
    for(int l = 0; l < Xi.cols(); l++) {
      if(!bim.snps_to_use[down + l])
        Xi.col(l) = Eigen::VectorXd::Constant(nindiv, 3);
    }

    count_non_missing(Xi, NMi, NMGi); // NMGi keeps track of missing values!
    swap_na_matrix(Xi); // all NA turns into zero
 
    // change Xi to genotype * effect_size
    for(int l = 0; l < Xi.cols(); l++) {
      Xi.col(l) = Xi.col(l) * my_betas.effects_in_order[down + l]; 
    }

    // effects per individual
    g_hat += Xi.rowwise().sum(); 

    calculate_grm2(Xi, Ai, NMGi);
    update_grm(A, NM, Ai, NMi);
  }

  bed.close();

  // scale by number of non-missing genotypes 
  for(int i = 0; i < nindiv; i++) {
     for(int j = 0; j < nindiv; j++) {
        A(i,j) /= NM(i,j) - 1;
     }
  }

  std::ofstream grm;
  std::ofstream nm;

  grm.open(grm_bin, std::ios::out| std::ios::trunc | std::ios::binary); 
  nm.open(grm_n, std::ios::out| std::ios::trunc | std::ios::binary); 
  
  for(int i = 0; i < nindiv; i++) {
    for(int j = 0; j <= i; j++) {
      float grm_ij = A(i,j);
      float nm_ij = NM(i,j);
      grm.write(reinterpret_cast<const char*>(&grm_ij), sizeof(float));
      nm.write(reinterpret_cast<const char*>(&nm_ij), sizeof(float));
    }
  }
  
  grm.close();
  nm.close();

  std::ofstream id_grm(grm_id);
  std::ofstream g_eff(g_pred);

  for(size_t i = 0; i < fam.individual_id.size(); i++) {
    id_grm << i + 1 << "\t" << fam.individual_id[i] << "\n";
    g_eff << fam.family_id[i] << "\t" << fam.individual_id[i] << "\t" << g_hat[i] << "\n";
  }

  return 0;
}





















