#include <Eigen/Dense>

void read_bed(char* , Eigen::MatrixXd&);
void calculate_grm(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);
double swap_na(double x, double NA, double tol);
void update_grm(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);
void scale_and_center_genotype(Eigen::MatrixXd&, Eigen::MatrixXd&);
void count_non_missing(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);
void swap_na_matrix(Eigen::MatrixXd&);
void calculate_grm2(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);
void swap_na_mm(Eigen::MatrixXd&, Eigen::MatrixXd&);
void calculate_grm3(Eigen::MatrixXd&, Eigen::MatrixXd&);
void read_bed2(std::ifstream&, Eigen::MatrixXd&);