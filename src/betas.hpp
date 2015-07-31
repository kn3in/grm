#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <map>
// #include "plink_data.hpp"

class bim_data;

class betas {
	public:
	// betas(); do I want default constructor w/o path?
	betas(std::string);
	~betas();
	void swap_betas(const bim_data&);

    std::vector<std::string> snps;
    std::vector<std::string> alleles;
    std::vector<double> effects;
    std::map<std::string, int> rs_id2index;

};