#include <fstream>
#include <stdexcept>
#include <vector>
#include <set>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <math.h>
#include "plink_data.hpp"
#include "betas.hpp"
#include "ld_data.hpp"

bim_data::bim_data(std::string path_to_bim_file) {
   
   std::ifstream bim_file(path_to_bim_file.c_str());
   if(!bim_file) throw std::runtime_error("Runtime error: can not open bim file [" + path_to_bim_file + "]");

   int chrom;
   std::string rs_id;
   double g_dist;
   int bp;
   std::string ref_all;
   std::string another_all;
   int index(0);

   while(bim_file) {
      
      bim_file >> chrom; // attempt to read
      if(bim_file.eof()) break; // are we past EOF?
      bim_chr.push_back(chrom);

      bim_file >> rs_id;
      bim_rs_id.push_back(rs_id);
      rs_id2index[rs_id] = index;
      index++;

      bim_file >> g_dist;
      bim_genetic_distance.push_back(g_dist);

      bim_file >> bp;
      bim_position.push_back(bp);

      bim_file >> ref_all;
      bim_ref_all.push_back(ref_all);

      bim_file >> another_all;
      bim_another_all.push_back(another_all);

   }
   
   bim_file.close();

   if(bim_rs_id.size() != rs_id2index.size()) {
    // map size will be smaller if duplicated rs ids used as keys
      throw std::runtime_error("Runtime error: duplicated SNP ids in the bim file");
   }
  
   nsnps = rs_id2index.size(); 
}

bim_data::~bim_data() {

}

void bim_data::find_overlap(const betas& b) {
    std::set<std::string> s1(bim_rs_id.begin(), bim_rs_id.end());
    std::set<std::string> s2(b.snps.begin(), b.snps.end());
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(snps_without_betas));
}

void bim_data::setup_snps_without_betas(const betas& b) {
 // how to make sure function has been called only once?
   bim_data::find_overlap(b);
   
   if(snps_without_betas.size() == (size_t)nsnps) 
      throw std::runtime_error("Runtime error: no overlapping SNPs in bim and betas files");
   
   std::cout << "Number of non-overlapping SNPs [present in bim, absent in betas file]: " << std::endl;
   std::cout << snps_without_betas.size() << std::endl;
}
 
void bim_data::setup_snps_to_iterate() {
   snps_to_use = std::vector<bool>(nsnps, true);
   
   for(size_t i = 0; i < snps_without_betas.size(); i++) {
      snps_to_use[ rs_id2index[ snps_without_betas[i] ] ] = false;
    }
}     

void bim_data::setup_snps_per_chunk(int chunk_size) {
    if(chunk_size > nsnps)
        chunk_size = nsnps;
    int chunks = ceil( (double)nsnps / chunk_size);
    int last_chunk = nsnps % chunk_size; 
    snps_per_chunk = std::vector<int>(chunks, chunk_size);
    if(last_chunk)
        snps_per_chunk.back() = last_chunk;
}

void bim_data::setup_snps_per_chunk2(const ld_data& ld) {
  // for now assume ld data is sorted
  // and bim data is sorted as well, too much to ask?
  
  

  for(size_t chunk = 0; chunk < ld.stop.size(); chunk++) {
    // if we pretend that we do not know anything about data being sorted
    std::vector<int> snp_index_in_a_chunk;
    
    for(int snp = 0; snp < nsnps; snp++) {
      if(bim_chr[snp] == ld.chr[chunk] &&
         bim_position[snp] >= ld.start[chunk] &&
         bim_position[snp] < ld.stop[chunk]) {
      
         snp_index_in_a_chunk.push_back(snp);
      }

    }
    
    chunks.push_back(snp_index_in_a_chunk);

  }
 
}

void bim_data::setup_number_of_snps_per_chunk() {
  for(size_t i = 0; i < chunks.size(); i++) {
    std::vector<int> chnk = chunks[i];
    if( !chnk.empty() ) {
        int n_snps = chnk.back() - chnk.front();
        snps_per_chunk.push_back(n_snps);
    }
  }
}

void bim_data::setup_running_sum_per_chunk() {
// to be used as offset for aligning snps betas (all in RAM) with the snps currently read in RAM
  running_sums.resize(snps_per_chunk.size());
  running_sums.front() = 0;
  std::partial_sum(snps_per_chunk.begin(), snps_per_chunk.end() - 1,  running_sums.begin() + 1);
}

fam_data::fam_data(std::string path_to_fam_file) {
   std::ifstream fam_file(path_to_fam_file.c_str());
   if(!fam_file) throw std::runtime_error("Runtime error: can not open fam file [" + path_to_fam_file + "]");

      std::string my_string;
      int my_int;

      int index(0);

      while(fam_file) {
      
      fam_file >> my_string; // attempt to read
      if(fam_file.eof()) break; // are we past EOF?
      family_id.push_back(my_string);

      fam_file >> my_string;
      individual_id.push_back(my_string);
      individual_id2index[my_string] = index;
      index++;

      fam_file >> my_string;
      father_id.push_back(my_string);

      fam_file >> my_string;
      mother_id.push_back(my_string);

      fam_file >> my_int;
      sex.push_back(my_int);

      fam_file >> my_int;
      phenotype.push_back(my_int);

   }
   
   fam_file.close();

   if(individual_id.size() != individual_id2index.size()) {
    // map size will be smaller if duplicated individula ids used as keys
      throw std::runtime_error("Runtime error: duplicated individual ids in the fam file");
      
   }

   nindiv = individual_id2index.size();

}

fam_data::~fam_data() {

}
































