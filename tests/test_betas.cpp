#include <fstream>
#include "lest.hpp"
#include "../src/betas.hpp"
#include "../src/plink_data.hpp"

const lest::test specification[] =
{
    CASE( "No file no results" )
    {
        EXPECT_THROWS( betas("File that does not exist") );
    },

    CASE( "Duplicated snps in the betas file get caught" )
    {
        EXPECT_THROWS( betas("duplicated_betas.txt") );
    },

    CASE( "Properly formatted file gets in" )
    {
    	SETUP() {

            std::ofstream test_betas_file("test_betas_file.txt");
            test_betas_file << "22-16105756\tA\t-0.000697339779725826"
                            << "22-16132491\tA\t0.000551860566130064"
                            << "22-16163632\tC\t-0.000715370670632225" << std::endl;
            test_betas_file.close();

            betas b("test_betas_file.txt");
            std::remove("test_betas_file.txt");

	        EXPECT(b.snps.size() == 3);
	        EXPECT(b.alleles.size() == 3);
	        EXPECT(b.effects.size() == 3);
	        EXPECT(b.rs_id2index.size() == 3);
        }
    },

    CASE( "Ref allele matching" )
    {
        SETUP() {
	        EXPECT_NO_THROW(bim_data("testing2.bim"));
            bim_data bim("testing2.bim");
	        betas b("not_duplicated_betas.txt");

            // bim.setup_snps_without_betas(b);
            // bim.setup_snps_to_iterate();
            // b.swap_betas(bim);

	        // EXPECT();


        }
    }

};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}