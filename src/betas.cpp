#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "betas.hpp"
#include "plink_data.hpp"

betas::betas(std::string full_path) {
   
   std::string snp_id;
   std::string ref_all;
   double eff;
   int index(0);
   
   std::ifstream eff_file(full_path.c_str());
   
   if(!eff_file) throw std::runtime_error("Runtime error: can not open file with betas [" + full_path + "]");
   
   while(eff_file) {
      
      eff_file >> snp_id;
      if(eff_file.eof()) break;
      snps.push_back(snp_id);
      rs_id2index[snp_id] = index;
      index++;

      eff_file >> ref_all;
      alleles.push_back(ref_all);

      eff_file >> eff;
      effects.push_back(eff);


   }
   
   eff_file.close();
   
   if(snps.size() != rs_id2index.size()) {
    // map size will be smaller if duplicated rs ids used as keys
      throw std::runtime_error("Runtime error: duplicated SNP ids in the betas file");

   }
}

betas::~betas() {
   
}

void betas::order_betas(const bim_data& bim) {
   effects_in_order = std::vector<double>(bim.nsnps, 0);
   std::vector<std::string> snps_with_non_matching_allels;

   for(size_t i = 0; i < bim.snps_to_use.size(); i++) {
      if(bim.snps_to_use[i]) {
         
         int beta_index = rs_id2index[ bim.bim_rs_id[i] ];
         effects_in_order[i] = effects[beta_index];

         if(bim.bim_ref_all[i] != alleles[beta_index]) {
                     
            // alternative should be the same as in betas file
            if(bim.bim_another_all[i] != alleles[beta_index]) {
               snps_with_non_matching_allels.push_back(bim.bim_rs_id[i]);
            }
            
            // swap sign for effect size
            effects_in_order[i] = -effects[beta_index];
         }

         
      }
   }
   
   if(snps_with_non_matching_allels.size() > 0) {

      std::ofstream non_matching_snp_rs_ids("rs_ids_of_snps_with_non_matching_allels.txt");

      for(size_t j = 0; j < snps_with_non_matching_allels.size(); j++) {
         non_matching_snp_rs_ids << snps_with_non_matching_allels[j] << "\n";
      }
   
      non_matching_snp_rs_ids.close();

      throw std::runtime_error("Runtime error: bad thing happend");
    }

}












