/*
 * svpToQubo.cpp
 *
 *  Created on: Apr 3, 2021
 *      Author: Milos Prokop
 *
 */

#include <iostream>
#include <gtest/gtest.h>
#include "../src/vqaConfig.h"

TEST(svpToQuboTest, binary_substitution_penalized){

	std::string config_file = "../tests/test_files/config_svpToQubo.txt";

	VqaConfig * vqaConfig;

	ASSERT_NO_THROW(vqaConfig = new VqaConfig(config_file));

	for(auto &lattice : vqaConfig->getLattices()){

		std::string generatedHamiltonian = lattice.toHamiltonianString(Lattice::x_zero_one);

	}
}

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	return ret;

}


