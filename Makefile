build:
	c++ -std=c++11 -fno-inline -g -Wall -Wextra \
	src/mmap_plink.cpp src/plink.cpp src/betas.cpp src/plink_data.cpp src/data.cpp \
	-o bin/bgrm -I eigen325
