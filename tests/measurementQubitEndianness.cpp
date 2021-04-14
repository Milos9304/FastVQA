/*
 * measurementQubitEndianness.cpp
 *
 *  Created on: Apr 14, 2021
 *      Author: Milos Prokop
 */

#include <gtest/gtest.h>
//#include "xacc_observable.hpp"
//#include "xacc_service.hpp"
#include "xacc.hpp"


TEST(measurementQubitEndianness, getMeasurementCounts_repated_measurement_mode){

	//q0 q1 q2 q3 q4
	// 1  0  0  1  1
	//
	// requires
	// output: "10011"
	//
	// i.e. the output is of form q_n-1, q_n-2, ..., 0

	xacc::setOption("quest-verbose", "true");

	auto qpu = xacc::getAccelerator("quest", {{"nbQbits", 5}, {"shots", 5},
										{"repeated_measurement_strategy", true}});

    auto provider = xacc::getIRProvider("quantum");
	auto qubitReg = xacc::qalloc(5);

    auto compiler = xacc::getCompiler("xasm");

	auto ir = compiler->compile(R"(__qpu__ void test(qbit q) {
		X(q[0]);
		X(q[3]);
		X(q[4]);
	})", qpu);

	auto program = ir->getComposite("test");
	std::vector<size_t> indices;
	for(size_t i=0; i < 5; ++i){
	  indices.push_back(i);
	}

	auto meas = provider->createInstruction("Measure", indices);
	program->addInstructions({meas});

	qpu->execute(qubitReg, program);

	ASSERT_EQ(qubitReg->getMeasurementCounts().size(), 5);

	for(std::pair<std::string, int> meas : qubitReg->getMeasurementCounts()){
		ASSERT_STREQ(meas.first.c_str(), "11001");
	}

}

int main(int argc, char **argv) {

	xacc::Initialize();
	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	xacc::Finalize();
	return ret;


}
