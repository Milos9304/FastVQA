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
#include <iomanip>
#include <gmpxx.h>
#include "../src/vqaConfig.h"

std::vector<std::string> split(std::string const &input) {
	    std::istringstream buffer(input);
	    std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)),
	                                 std::istream_iterator<std::string>());
	    return ret;
}

std::pair<mpq_class, std::vector<std::pair<mpq_class, std::vector<std::string>>>> parseMatlabHml(std::string matlabHml){

	int sign = 1;
	mpq_class constant = 0;
	std::vector<std::pair<mpq_class, std::vector<std::string>>> result;

	for(auto &token : split(matlabHml)){
		//std::cout << token << "\n";

		if (token.find('*') != std::string::npos){

			std::vector<std::string> vars;
			mpq_class coeff = 0;

			bool c = true;
			size_t pos;
			std::string parsed;
			std::string remaining = token;
			while ((pos = remaining.find('*')) != std::string::npos) {
				parsed =  remaining.substr(0, pos);
				remaining = remaining.substr(pos+1, remaining.size()-1);
				if(c){
					c = false;
					coeff = std::stod(parsed);
				}
				else
					vars.push_back(parsed);
			}

			vars.push_back(remaining);

			coeff = sign * coeff;
			for(unsigned int i = 0; i < vars.size(); ++i){
				if(vars[i][0]!='Z')
					throw std::runtime_error("Unknown term");

				vars[i]="Z"+std::to_string(std::stoi(vars[i].substr(1, vars[i].size()-1))-1);
			}


			result.push_back(std::pair<mpq_class, std::vector<std::string>>(coeff, vars));
			//std::cout << coeff << "\n";
			//for(auto &var:vars)
			//	std::cout << var << "\n";


		}else{
			if(token == "-"){
				sign = -1;
			}else if(token == "+")
				sign = 1;
			else
				constant = std::stod(token);
		}
	}

	return std::pair<mpq_class, std::vector<std::pair<mpq_class, std::vector<std::string>>>>(constant, result);

}

std::pair<mpq_class, std::vector<std::pair<mpq_class, std::vector<std::string>>>> parseGenHml(std::string genHml){

	int sign = 1;
	bool c = true;
	mpq_class constant = 0;
	mpq_class coeff = 0;
	std::vector<std::pair<mpq_class, std::vector<std::string>>> result;

	std::vector<std::string> vars;

	for(auto &token : split(genHml)){

		//std::cout << token << "\n";

		if(token == "-"){
			if(vars.size() > 0){
				result.push_back(std::pair<mpq_class, std::vector<std::string>>(sign * coeff, vars));
				vars.clear();
			}
			sign = -1;
			continue;
		}else if(token == "+"){
			if(vars.size() > 0){
				result.push_back(std::pair<mpq_class, std::vector<std::string>>(sign * coeff, vars));
				vars.clear();
			}
			sign = 1;
			continue;
		}

		size_t pos;
		if((pos = token.find('Z')) == std::string::npos){ //coeff
			coeff = std::stod(token);
			if(c){
				constant = coeff;
				c = false;
			}
		}else{
			vars.push_back(token);
		}

	}

	if(vars.size() > 0){
		result.push_back(std::pair<mpq_class, std::vector<std::string>>(sign * coeff, vars));
		vars.clear();
	}

	return std::pair<mpq_class, std::vector<std::pair<mpq_class, std::vector<std::string>>>>(constant, result);

}

bool compareHamiltonians(std::string generatedHml, std::string matlabHml){

	std::pair<mpq_class, std::vector<std::pair<mpq_class, std::vector<std::string>>>> parsedMatlab = parseMatlabHml(matlabHml);
	std::pair<mpq_class, std::vector<std::pair<mpq_class, std::vector<std::string>>>> parsedGenHml = parseGenHml(generatedHml);

	mpq_class matlabCoeff = parsedMatlab.first;
	std::vector<std::pair<mpq_class, std::vector<std::string>>> matlabTerms = parsedMatlab.second;

	/*std::cout << matlabCoeff << "\n";
	for(auto &term:matlabTerms){
		std::cout << term.first << "* ";
		for(auto &var : term.second)
			std::cout << var << " ";
		std::cout << "\n";
	}*/

	mpq_class genCoeff = parsedGenHml.first;
	std::vector<std::pair<mpq_class, std::vector<std::string>>> genTerms = parsedGenHml.second;

	/*std::cout << genCoeff << "\n";
	for(auto &term:genTerms){
		std::cout << term.first << "* ";
		for(auto &var : term.second)
			std::cout << var << " ";
		std::cout << "\n";
	}*/

	//logd(std::to_string(matlabCoeff));
	//logd(std::to_string(genCoeff));

	//std::cout << matlabCoeff << " " << genCoeff << "\n";

	if(matlabCoeff != genCoeff){
		std::cerr<<"unequal constants\n";
		return false;
	}

	if(matlabTerms.size() != genTerms.size()){
		std::cerr<<"unequal sizes\n";\
		return false;
	}


	for(auto &genTerm : genTerms){

		mpq_class genCoeff = genTerm.first;
		bool found = false;

		mpq_class matCoeff;
		std::vector<std::string> genVars;
		std::vector<std::string> matVars;

		genVars = genTerm.second;

		for(auto &matTerm : matlabTerms){

			matCoeff = matTerm.first;
			matVars = matTerm.second;

			if(matCoeff == genCoeff){ //potential match

				if(genVars.size() != matVars.size())
					continue;

				if(genVars.size() == 1){

					if(genVars[0] == matVars[0]){
						found = true;
						break;
					}
					else
						continue;
				}

				//size=2
				std::string g0 = genVars[0];
				std::string m0 = matVars[0];
				std::string g1 = genVars[1];
				std::string m1 = matVars[1];

				if(g0 != m0){
					m0 = matVars[1];
					m1 = matVars[0];
				}

				if(g0 == m0 && g1 == m1){
					found = true;
					break;
				}else
					continue;
			}

		}
		if(!found){

			std::stringstream ss;
			ss << std::fixed << std::setprecision(2) << mpf_class(genCoeff);

			if(genVars.size() == 1)
				loge(ss.str() + " " + genVars[0] + " not found in matlab");
			else
				loge(ss.str() + " " + genVars[0] + " " + genVars[1] + " not found in matlab");
			return false;
		}
	}

	return true;
}

/*TEST(svpToQuboTest, binary_substitution_penalized){

	int penalty = 1000; //same as in matlab file

	std::string config_file = "../tests/test_files/config_svpToQubo.txt";

	VqaConfig * vqaConfig;

	ASSERT_NO_THROW(vqaConfig = new VqaConfig(config_file));

	for(auto &lattice : vqaConfig->getLattices()){

		std::string generatedHamiltonian = lattice.toHamiltonianString(new MapOptions(MapOptions::x_symmetric,
				MapOptions::naive_overapprox, MapOptions::penalty_all, penalty, 1));

		std::string matlabHamiltonian;
		std::ifstream file("../tests/test_files/lattices/matlab_output/"+lattice.name+".txt");
		std::cout << "../tests/test_files/lattices/matlab_output/"+lattice.name+".txt"<<"\n";
		ASSERT_TRUE(file.is_open());

		std::getline(file, matlabHamiltonian);
		file.close();

		//std::cout <<generatedHamiltonian << "\n";

		EXPECT_TRUE(compareHamiltonians(generatedHamiltonian, matlabHamiltonian));

	}
}*/

TEST(svpToQuboTest, two_qubit_substitution_penalized){

	int penalty = 1000; //same as in matlab file

	std::string config_file = "../tests/test_files/config_svpToQubo.txt";

	VqaConfig * vqaConfig;

	ASSERT_NO_THROW(vqaConfig = new VqaConfig(config_file));

	for(auto &lattice : vqaConfig->getLattices()){

		std::string generatedHamiltonian = lattice.toHamiltonianString(new MapOptions(MapOptions::x_symmetric,
				MapOptions::naive_overapprox, MapOptions::penalty_all, penalty, 1), false);

		std::string matlabHamiltonian;
		std::ifstream file("../tests/test_files/lattices/matlab_output/"+lattice.name+".txt");
		std::cout << "../tests/test_files/lattices/matlab_output/"+lattice.name+".txt"<<"\n";
		ASSERT_TRUE(file.is_open());

		std::getline(file, matlabHamiltonian);
		file.close();

		//std::cout <<generatedHamiltonian << "\n";

		EXPECT_TRUE(compareHamiltonians(generatedHamiltonian, matlabHamiltonian));

		break;

	}
}

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	return ret;

}


