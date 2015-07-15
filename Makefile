debug:
	/usr/local/bin/c++-4.9 -std=c++11 -fno-inline -g -Wall src/mmap_plink.cpp src/plink.cpp -o bin/mmap_plink -I eigen325

release:
	/usr/local/bin/c++-4.9 -std=c++11 -O3 src/mmap_plink.cpp src/plink.cpp -o bin/mmap_plink -I eigen325

openmpd:
	/usr/local/bin/c++-4.9 -std=c++11 -fopenmp -fno-inline -g -Wall src/mmap_plink.cpp src/plink.cpp -o bin/mmap_plink -I eigen325

openmpr:
	/usr/local/bin/c++-4.9 -std=c++11 -fopenmp -O3 src/mmap_plink.cpp src/plink.cpp -o bin/mmap_plink -I eigen325