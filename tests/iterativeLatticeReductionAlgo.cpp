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
#include "../src/lattice/hmlLattice.hpp"

using namespace indicators;

TEST(iterativeLatticeReductionAlgo, gaussian_heuristics){

	std::string config_file = "../tests/test_files/config_iterativeLatticeReduction.txt";
	VqaConfig * vqaConfig;
	ASSERT_NO_THROW(vqaConfig = new VqaConfig(config_file));

	//value calculated by fpylll
	ASSERT_NEAR(vqaConfig->getLattices()[0].get_orig_gh().get_d(), 22.61827991781903, 0.001);


}

TEST(iterativeLatticeReductionAlgo, opt_config_valid){

	//xacc::setOption("quest-verbose", "true");
	//xacc::setOption("quest-debug", "true");

	const int num_repetitions = 5;
	const int max_iters = 5;

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
	qaoaOptions.max_iters = max_iters;//000;
	qaoaOptions.calcVarAssignment = true;
	qaoaOptions.verbose = true;
	qaoaOptions.extendedParametrizedMode=true;
	qaoaOptions.accelerator = accelerator;
	qaoaOptions.optimizer = optimizer;

	ProgressBar bar;
	bar.disabled=true;

	qaoaOptions.set_default_stats_function(execStats, &bar);

	QOracle quantum_oracle = [execStats, &bar, &qaoaOptions]
							  (xacc::qbit** buffer, std::string hamiltonian, std::string name) {
		Qaoa::run_qaoa(buffer, hamiltonian, name, &bar, execStats, &qaoaOptions);
	};

	for(auto &lattice : vqaConfig->getLattices()){

		IterativeLatticeReduction ilr(&lattice, mapOptions, quantum_oracle, 1);

		for(int i = 0; i < num_repetitions; ++i){

			std::pair<VectorInt, double> x_vect_energy = ilr.run_test();

			VectorInt x_vect = x_vect_energy.first;
			double energy = x_vect_energy.second;

			VectorInt short_lattice_vector = ilr.xVectToShortVect(&x_vect);

			int len_squared = 0;

			for(auto &x:short_lattice_vector){
				len_squared += pow(mpz_class(x).get_si(), 2);
			}

			ASSERT_EQ(len_squared, energy);

		}
	}

}

int main(int argc, char **argv) {

	xacc::Initialize();
	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	xacc::Finalize();
	return ret;


}



