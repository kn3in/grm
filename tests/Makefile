.PHONY: test_betas
	
build_test_betas:
	g++ -Wall -Wextra -std=c++11 -I../lest -I../eigen325  \
	-o test_betas test_betas.cpp ../src/betas.cpp ../src/plink_data.cpp

test_betas:
	./test_betas -p	