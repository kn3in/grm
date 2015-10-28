#pragma once

#include <vector>

class ld_data {
	public:
		ld_data(std::string);
		~ld_data();

		std::vector<int> chr;
		std::vector<int> start;
		std::vector<int> stop;
};