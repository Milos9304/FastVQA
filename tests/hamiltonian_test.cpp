/*
 * hamiltonian_test.cpp
 *
 *  Created on: June 18, 2022
 *      Author: Milos Prokop
 */

#include <gtest/gtest.h>
#include <cmath>
#include "fastVQA.h"

#define GTEST_COUT std::cerr << "[          ] [ INFO ]"

TEST(hamiltonian_test, initializeMinusSigmaXHamiltonian){


}

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	return ret;


}
