#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "plink.hpp"
#include "betas.hpp"
#include "plink_data.hpp"
#include "data.hpp"

int main (int argc, char* argv[]) {
  
  std::string usage("USAGE: bgrm plink_data betas_data out_name");
    
  if(argc != 4) {
    std::cerr << usage << std::endl;
      return 1;
  }
  
  std::string plink_base = argv[1];
  std::string my_betas_path = argv[2];
  std::string results_file = argv[3];

  std::string bim_file = plink_base + ".bim";
  std::string fam_file = plink_base + ".fam";
  std::string bed_file = plink_base + ".bed";

  betas my_betas(my_betas_path);
  bim_data bim(bim_file);
  fam_data fam(fam_file);
  data results(fam);

  bim.setup_snps_without_betas(my_betas);
  bim.setup_snps_to_iterate();

  // read in chunks but not overdo it
  bim.setup_snps_per_chunk(100000);

  my_betas.order_betas(bim);
  
  std::ifstream bed;
  bed.open (bed_file, std::ios::in | std::ios::binary); 
  bed.seekg(3, std::ifstream::beg);
  
  for(size_t i = 0; i < bim.snps_per_chunk.size(); ++i) {
    Eigen::MatrixXd Xi(fam.nindiv, bim.snps_per_chunk[i]); //current Genotype
    read_bed2(bed, Xi);
    results.update(Xi, i, bim, my_betas);
  }

  bed.close();
  results.scale_by_nm();
  results.save(results_file);
}





















