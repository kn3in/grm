#include <string>
#include <fstream>
#include <Eigen/Dense>
#include "data.hpp"
#include "betas.hpp"

data::data(const fam_data& fam) {
	nindiv = fam.nindiv;
	family_id = fam.family_id;
	individual_id = fam.individual_id;
	A = Eigen::MatrixXd::Zero(nindiv, nindiv);
	NM = Eigen::MatrixXd::Zero(nindiv, nindiv);
	g_hat = Eigen::VectorXd::Zero(nindiv);
	
}

data::~data() {

}

void data::update(Eigen::MatrixXd& Xi, int i, const bim_data& bim, const betas& my_betas) {
	Eigen::MatrixXd NMGi(nindiv, bim.snps_per_chunk[i]); //current non-missing->1 NA->0
    Eigen::MatrixXd Ai(nindiv, nindiv); // current GRM
    Eigen::MatrixXd NMi(nindiv, nindiv); // current Number of non-missing SNPs

    int down = i * bim.snps_per_chunk.front();
    
    for(int l = 0; l < Xi.cols(); l++) {
      if(!bim.snps_to_use[down + l])
        Xi.col(l) = Eigen::VectorXd::Constant(nindiv, 3);
    }

    count_non_missing(Xi, NMi, NMGi); // NMGi keeps track of missing values!
    swap_na_matrix(Xi); // all NA turns into zero
 
    // change Xi to genotype * effect_size
    for(int l = 0; l < Xi.cols(); l++) {
      Xi.col(l) = Xi.col(l) * my_betas.effects_in_order[down + l]; 
    }

    // effects per individual
    g_hat += Xi.rowwise().sum(); 

    calculate_grm2(Xi, Ai, NMGi);
    update_grm(Ai, NMi);
}

void data::scale_by_nm() {
	for(int i = 0; i < nindiv; i++) {
		for(int j = 0; j < nindiv; j++) {
			A(i,j) /= NM(i,j) - 1;
		}
	}
}

void data::count_non_missing(Eigen::MatrixXd& Xi, Eigen::MatrixXd& NMi, Eigen::MatrixXd& NMGi) {
   int N = Xi.rows();
   int n = Xi.cols();
   
   // NMG = Eigen::MatrixXd::Zero(N, nsnps);
   
   // missing value is zero, one otherwise.
   for(int i = 0; i < N; i++) {
      for(int j = 0; j < n; j++) {
         if(Xi(i,j) != 3.0) {
            NMGi(i,j) = 1.0;
         } else {
            NMGi(i,j) = 0.0;
         }
      }
   }

   NMi = NMGi * NMGi.transpose();
}

void data::swap_na_matrix(Eigen::MatrixXd& Xi) {
   int n = Xi.rows();
   for(int i = 0; i < Xi.cols(); i++) {
      std::for_each(Xi.col(i).data(), Xi.col(i).data() + n, [](double &geno) {
                  if(geno == 3)
                     geno = 0; // swap NA representation from 3 to 0
               });
   }
}

void data::calculate_grm2(Eigen::MatrixXd& Xi, Eigen::MatrixXd& Ai, Eigen::MatrixXd& NMGi) {
   // Eigen::VectorXd NM_per_column = NMGi.colwise().sum

   // here is an easy way of getting bugs MAKE SURE AT LEAST ONE SNP IS NOT MISSING  
   Eigen::VectorXd Mean = Xi.colwise().sum().array() / NMGi.colwise().sum().array();
   
   // center
   Xi.rowwise() -= Mean.transpose();
   // NA is not zero anymore so bring it back to zero!!!
   swap_na_mm(Xi, NMGi);
   
   Eigen::VectorXd Sd = Xi.colwise().norm().array() / (NMGi.colwise().sum().array() - 1).sqrt();

   // scale
   for(int i = 0; i < Xi.cols(); i++) {
      if(Sd[i] > 1e-06) {
         Xi.col(i) /= Sd[i];
      }
   }

   // swap is not nessesary cause 0 divided by sd is still 0
   //swap_na_mm(Xi, NMGi);

   // Eigen::MatrixXd Yi = Xi.transpose();
   Eigen::VectorXd Xi_Mean = Xi.rowwise().sum().array() / NMGi.rowwise().sum().array();
   Xi.colwise() -= Xi_Mean;
   swap_na_mm(Xi, NMGi);
   Ai = Xi * Xi.transpose();
}

void data::swap_na_mm(Eigen::MatrixXd& Xi, Eigen::MatrixXd& NMGi) {
   for(int i = 0; i < Xi.cols(); i++) {
      for(int j = 0; j < Xi.rows(); j++) {
         if(!NMGi(j, i))
            Xi(j, i) = 0.0;
      }
   }
}

void data::update_grm(Eigen::MatrixXd& Ai, Eigen::MatrixXd& NMi) {
   A  += Ai;
   NM += NMi;
}

void data::save(std::string path_to) {
	std::string grm_bin = path_to + ".grm.bin";
	std::string grm_n = path_to + ".grm.N.bin";
	std::string grm_id = path_to + ".grm.id";
	std::string g_pred = path_to + "_gpredictor.txt";

	std::ofstream grm;
    std::ofstream nm;

    grm.open(grm_bin, std::ios::out| std::ios::trunc | std::ios::binary); 
    nm.open(grm_n, std::ios::out| std::ios::trunc | std::ios::binary); 
  
	for(int i = 0; i < nindiv; i++) {
		for(int j = 0; j <= i; j++) {
			float grm_ij = A(i,j);
			float nm_ij = NM(i,j);
			grm.write(reinterpret_cast<const char*>(&grm_ij), sizeof(float));
			nm.write(reinterpret_cast<const char*>(&nm_ij), sizeof(float));
		}
	}
  
    grm.close();
    nm.close();

	std::ofstream id_grm(grm_id);
	std::ofstream g_eff(g_pred);

	for(size_t i = 0; i < individual_id.size(); i++) {
		id_grm << family_id[i] << "\t" << individual_id[i] << "\n";
		g_eff << family_id[i] << "\t" << individual_id[i] << "\t" << g_hat[i] << "\n";
	}
	
	id_grm.close();
	g_eff.close();
}



























