#include <iostream>
#include <string>
#include <set>
#include <algorithm>
#include <map>
#include <stdexcept>
#include "betas.hpp"
#include "plink_data.hpp"


int main(int argc, char** argv) {

    std::string usage("USAGE: mr_grm plink_data betas_data out_name");
    
    if(argc != 4) {
    	std::cerr << usage << std::endl;
        return 1;
    }
    
	std::string plink_base = argv[1];
	std::string bim_file = plink_base + ".bim";
	std::string fam_file = plink_base + ".fam";
    std::string my_betas_path = argv[2];
	std::string results_file = argv[3];

	betas test1(my_betas_path);
	bim_data bim(bim_file);
    fam_data fam(fam_file);

    bim.setup_snps_without_betas(test1);
    bim.setup_snps_to_iterate();
    test1.swap_betas(bim);

    
}

















