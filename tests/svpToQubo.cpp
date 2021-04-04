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

bool compareHamiltonians(std::string generatedHml, std::string matlabHml, std::function<std::string(std::string)> hmlMap){

	int sign = 1;

	for(auto &token : split(matlabHml)){
		//std::cout << token << "\n";

		int constant = 0;

		if (token.find('*') != std::string::npos){

			std::cout << token << ":" << std::endl;
			std::vector<std::string> vars;
			int coeff;

			bool c = true;
			size_t pos;
			std::string parsed;
			std::string remaining = token;
			while ((pos = remaining.find('*')) != std::string::npos) {
			    parsed =  remaining.substr(0, pos);
				remaining = remaining.substr(pos+1, remaining.size()-1);
			    if(c){
			    	c = false;
					coeff = std::stoi(parsed);
			    }
			    else
			    	vars.push_back(parsed);
			}

			if(remaining.find('^') != std::string::npos){
				pos = remaining.find('^');
				vars.push_back(remaining.substr(0, pos));
				vars.push_back(remaining.substr(0, pos));
			}
			else
				vars.push_back(remaining);

			coeff = sign * coeff;
			for(auto &v:vars){
				std::cout << v << "\n";
			}

		}else{
			if(token == "-"){
				sign = -1;
			}else if(token == "+")
				sign = 1;
			else
				constant = std::stoi(token);
		}

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

		std::cout <<generatedHamiltonian;

		EXPECT_TRUE(compareHamiltonians(generatedHamiltonian, matlabHamiltonian,
				[](std::string x) {
				return "x_"+std::to_string(std::stoi(x.substr(1, x.size()-1))-1)+"_b0";
		}));

		break;

	}
}

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	return ret;

}


