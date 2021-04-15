/*
 * iterativeLatticeReductionAlgo.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */


#include <gtest/gtest.h>
#include "../src/latticeAlgorithms/iterativeLatticeReduction.h"
#include "../src/executionStatistics.h"
#include "../src/vqaConfig.h"
#include "../src/qaoa/qaoa.h"
#include "../src/indicators/progress_bar.hpp"

using namespace indicators;

TEST(iterativeLatticeReductionAlgo, opt_config_valid){

	//xacc::setOption("quest-verbose", "true");
	//xacc::setOption("quest-debug", "true");

	std::string config_file = "../tests/test_files/config_iterativeLatticeReduction.txt";
	VqaConfig * vqaConfig;
	ASSERT_NO_THROW(vqaConfig = new VqaConfig(config_file));

	MapOptions* mapOptions = new MapOptions();
	mapOptions->verbose = false;

	ExecutionStatistics* execStats = new ExecutionStatistics();

	AcceleratorPartial accelerator = [](std::shared_ptr<xacc::Observable> observable) {
				return xacc::getAccelerator("quest", {std::make_pair("nbQbits", observable->nBits()),
						 // Doesn't require to prepare the same circuit over and over again, but needs to clone statevect.
						 std::make_pair("repeated_measurement_strategy", true)});
	};

	OptimizerPartial optimizer = [](std::vector<double> initialParams, int max_iters) {
		return xacc::getOptimizer("nlopt", xacc::HeterogeneousMap {std::make_pair("initial-parameters", initialParams),
																   std::make_pair("nlopt-maxeval", max_iters)});
	};

	QAOAOptions qaoaOptions;
	qaoaOptions.max_iters = 5;//000;
	qaoaOptions.calcVarAssignment = true;
	qaoaOptions.verbose = true;
	qaoaOptions.parametrizedMode=true;
	qaoaOptions.accelerator = accelerator;
	qaoaOptions.optimizer = optimizer;

	ProgressBar bar{
			option::BarWidth{50},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::ForegroundColor{Color::blue},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
			option::MaxProgress{qaoaOptions.max_iters}
	};

	qaoaOptions.set_default_stats_function(execStats, &bar);

	QOracle quantum_oracle = [execStats, &bar, &qaoaOptions]
							  (xacc::qbit** buffer, std::string hamiltonian, std::string name) {
		run_qaoa(buffer, hamiltonian, name, &bar, execStats, &qaoaOptions);
	};

	for(auto &lattice : vqaConfig->getLattices()){

		IterativeLatticeReduction ilr(&lattice, mapOptions, quantum_oracle, 1);
		ilr.run_test();

	}

	/*auto qpu = xacc::getAccelerator("quest", {{"nbQbits", 5}, {"shots", 5}});

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

	ASSERT_EQ(qubitReg->getMeasurementCounts().size(), 1); // because only one unique

	for(std::pair<std::string, int> meas : qubitReg->getMeasurementCounts()){
		ASSERT_STREQ(meas.first.c_str(), "11001");
	}*/

}

int main(int argc, char **argv) {

	xacc::Initialize();
	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	xacc::Finalize();
	return ret;


}



