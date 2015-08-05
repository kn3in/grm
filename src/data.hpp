#pragma once

#include <vector>
#include <string>
#include <Eigen/Dense>
#include "plink_data.hpp"

class data {
	public:
		data(const fam_data& fam);
		~data();

		void update(Eigen::MatrixXd& Xi, int i, const bim_data& bim, const betas& my_betas);
		void scale_by_nm();
		void save(std::string);
		void count_non_missing(Eigen::MatrixXd& Xi, Eigen::MatrixXd& NMi, Eigen::MatrixXd& NMGi);
        void swap_na_matrix(Eigen::MatrixXd& Xi);
		void calculate_grm2(Eigen::MatrixXd& Xi, Eigen::MatrixXd& Ai, Eigen::MatrixXd& NMGi);
        void swap_na_mm(Eigen::MatrixXd& Xi, Eigen::MatrixXd& NMGi);
        void update_grm(Eigen::MatrixXd& Ai, Eigen::MatrixXd& NMi);

		Eigen::MatrixXd A;       // GRM
		Eigen::MatrixXd NM;      // Pairwise non-missing
		Eigen::VectorXd g_hat;   // genetic effects
		int nindiv;
		std::vector<std::string> family_id;
		std::vector<std::string> individual_id;
};