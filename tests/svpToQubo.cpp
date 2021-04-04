/*
 * svpToQubo.cpp
 *
 *  Created on: Apr 3, 2021
 *      Author: Milos Prokop
 *
 */
#include <vector>
#include <fstream>
#include <iostream>
#include <gtest/gtest.h>
#include "../src/vqaConfig.h"

std::vector<std::string> split(std::string const &input) {
	    std::istringstream buffer(input);
	    std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)),
	                                 std::istream_iterator<std::string>());
	    return ret;
}

bool compareHamiltonians(std::string generatedHml, std::string matlabHml){

	for(auto &token : split(matlabHml)){
		std::cout << token << "\n";

	}



	return true;
}

TEST(svpToQuboTest, binary_substitution_penalized){

	std::string config_file = "../tests/test_files/config_svpToQubo.txt";

	VqaConfig * vqaConfig;

	ASSERT_NO_THROW(vqaConfig = new VqaConfig(config_file));

	for(auto &lattice : vqaConfig->getLattices()){

		std::string generatedHamiltonian = lattice.toHamiltonianString(Lattice::x_zero_one);

		std::string matlabHamiltonian;
		std::ifstream file("../tests/test_files/lattices/matlab_output/"+lattice.name+".txt");
		std::cout << "../tests/test_files/lattices/matlab_output/"+lattice.name+".txt"<<"\n";
		ASSERT_TRUE(file.is_open());

		std::getline(file, matlabHamiltonian);
		file.close();

		EXPECT_TRUE(compareHamiltonians(generatedHamiltonian, matlabHamiltonian));

		break;

	}
}

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	return ret;

}


