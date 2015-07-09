#include <Eigen/Dense>

void read_bed(char* , Eigen::MatrixXd&);
void calculate_grm(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);
double swap_na(double x, double NA, double tol);
void update_grm(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);