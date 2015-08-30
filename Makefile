build:
	c++ -std=c++11 -fno-inline -g -Wall -Wextra \
	src/mmap_plink.cpp src/plink.cpp src/betas.cpp src/plink_data.cpp src/data.cpp \
	-o bin/bgrm -I eigen325

release_cluster:
	/clusterdata/apps/gcc-4.7.4/bin/c++ -std=c++11 -O3 \
	    src/mmap_plink.cpp src/plink.cpp src/betas.cpp \
	    src/plink_data.cpp src/data.cpp -o bin/bgrm2 -I eigen325