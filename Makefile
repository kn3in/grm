build:
	c++ -std=c++11 -fno-inline -g -Wall src/mmap_plink.cpp src/plink.cpp -o bin/mmap_plink -I eigen325
