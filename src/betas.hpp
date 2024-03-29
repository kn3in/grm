#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <map>

class bim_data;

class betas {
	public:
	betas(std::string);
	~betas();
	void order_betas(const bim_data&);

    std::vector<std::string> snps;
    std::vector<std::string> alleles;
    std::vector<double> effects;
    std::map<std::string, int> rs_id2index;
    std::vector<double> effects_in_order;
};