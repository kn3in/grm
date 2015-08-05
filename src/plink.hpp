#pragma once

#include <fstream>
#include <Eigen/Dense>

void read_bed2(std::ifstream&, Eigen::MatrixXd&);