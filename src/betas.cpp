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

void betas::swap_betas(const bim_data& bim) {
   for(int i = 0; i < bim.snps_to_use.size(); i++) {
      if(bim.snps_to_use[i]) {
         int beta_index = rs_id2index[ bim.bim_rs_id[i] ];

            if(bim.bim_ref_all[i] != alleles[beta_index]) {
               // alternative should be the same as in betas file
               if(bim.bim_another_all[i] != alleles[beta_index]) {
                  throw std::runtime_error("Runtime error: neither reference nor alternative alleles for SNP: " + bim.bim_ref_all[i] + "agree with the allele provided in the betas file");
               }
               
               // swap sign for effect size otherwise do nothing
               effects[beta_index] = -effects[beta_index];
            }

      }
    }

}












