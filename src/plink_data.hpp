#pragma once
#include <vector>
#include <string>
#include <map>

class betas;

class bim_data {
	public:
		bim_data(std::string);
	    ~bim_data();
	    void find_overlap(const betas&);
	    void setup_snps_without_betas(const betas&);
        void setup_snps_to_iterate();
		void setup_snps_per_chunk(int chunk_size);

		std::vector<int> bim_chr;
		std::vector<std::string> bim_rs_id;
		std::vector<double> bim_genetic_distance;
		std::vector<int> bim_position;
		std::vector<std::string> bim_ref_all;
		std::vector<std::string> bim_another_all;
		std::map<std::string, int> rs_id2index;
		int nsnps;
		std::vector<bool> snps_to_use;
		std::vector<std::string> snps_without_betas;
		std::vector<int> snps_per_chunk;
};

class fam_data {
	public:
		fam_data(std::string);
		~fam_data();

		std::vector<std::string> family_id;
		std::vector<std::string> individual_id;
		std::vector<std::string> father_id;
		std::vector<std::string> mother_id;
		std::vector<int> sex; // 1 male 2 female
		std::vector<int> phenotype;
		std::map<std::string, int> individual_id2index;
		int nindiv;
};

