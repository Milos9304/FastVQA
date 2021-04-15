#include <iostream>
#include <iterator>
#include <exception>
#include <chrono>

#include "popl.hpp"
#include "logger.h"
#include "executionStatistics.h"
#include "qaoa/qaoa.h"
#include "vqaConfig.h"
#include "indicators/progress_bar.hpp"
#include "latticeAlgorithms/iterativeLatticeReduction.h"
#include <xacc.hpp>

#include "xacc_observable.hpp" //del
#include "PauliOperator.hpp" //del

using namespace popl;
using namespace indicators;

Color colors[9] = {Color::grey, Color::red, Color::green, Color::yellow, Color::blue, Color::magenta, Color::cyan, Color::white};

int main(int ac, char** av){

	OptionParser op("Allowed options");
	auto help_option = op.add<Switch>("h", "help", "produce help message");
    auto qaoa 		 = op.add<Switch>("", "qaoa", "run qaoa algorithm");
    auto config 	 = op.add<Value<std::string>>("", "config", "set config file location");

    op.parse(ac, av);

    if (help_option->is_set())
    		std::cout << op << "\n";

    else if(qaoa -> is_set()){

    	if(!config->is_set()){
    		loge("No config file specified");
    		return 1;
    	}

    	VqaConfig* vqaConfig;

    	try{
    		vqaConfig = new VqaConfig(config->value());
    	}
    	catch(std::exception &e){
    	    loge(e.what());
    	    return 1;
    	}

		logi("Running QAOA. Configuration loaded from " + config->value());
		int num_lattices = vqaConfig->getLattices().size();

		logi(std::to_string(num_lattices) + " lattice(s) in the config file");

		int i = 0;
		/*for(auto &lattice : vqaConfig->getLattices()){

			logi("Running lattice " + std::to_string(i++) + " / " + std::to_string(num_lattices));
			lattice.toHamiltonianString(Lattice::x_zero_one, true);
			//run_qaoa(lattice.toHamiltonianString(Lattice::x_zero_one), vqaConfig->verbose);

		}*/

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
		qaoaOptions.max_iters = 5000;
		qaoaOptions.verbose = true;
		qaoaOptions.optimizer = optimizer;
		qaoaOptions.accelerator = accelerator;
		qaoaOptions.parametrizedMode = true;

		MapOptions* mapOptions = new MapOptions();
		mapOptions->verbose = false;

		//std::string hamiltonian = vqaConfig->getLattices()[1].toHamiltonianString(options, true);

		ExecutionStatistics* execStats = new ExecutionStatistics();
		xacc::Initialize();

		int counter = 0;

		for(auto &lattice : vqaConfig->getLattices()){

			if(lattice.name != "q_5_1_20_b_0")
				continue;

			ProgressBar bar{
				option::BarWidth{50},
				option::Start{"["},
				option::Fill{"="},
				option::Lead{">"},
				option::Remainder{" "},
				option::End{"]"},
				option::PrefixText{std::to_string(counter+1) + "/" + std::to_string(num_lattices) + " " + lattice.name},
				option::ForegroundColor{colors[counter % 9]},
				option::ShowElapsedTime{true},
				option::ShowRemainingTime{true},
				option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
				option::MaxProgress{qaoaOptions.max_iters}
			};

			//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

				qaoaOptions.set_default_stats_function(execStats, &bar);

				QOracle quantum_oracle = [&bar, execStats, &qaoaOptions]
										  (xacc::qbit** buffer, std::string hamiltonian, std::string name) {
					run_qaoa(buffer, hamiltonian, name, &bar, execStats, &qaoaOptions);
				};

				IterativeLatticeReduction ilr(&lattice, mapOptions, quantum_oracle, 1);
				ilr.run();

			//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

			counter++;
		}

		xacc::Finalize();
		/*auto observable = xacc::quantum::getObservable("pauli", vqaConfig->getLattices()[0].toHamiltonianString(options));
        std::cout << "obs1 " << observable->toString() << "\n";

        std::map<int, std::pair<std::string, std::complex<double>>> operators;

        operators.emplace(0, std::pair<std::string, std::complex<double>>("Z", std::complex<double>(5,0)));
        operators.emplace(1, std::pair<std::string, std::complex<double>>("Z", std::complex<double>(6,0)));

        auto obs2 = xacc::quantum::PauliOperator(operators);
        std::cout << "obs2 " << obs2.toString() << "\n";*/

		return 0;
    }

    else{
    	std::cout << "Invalid argument. Use ./fastVQA -h to see help." << "\n";
    	return 1;
    }

    return 0;
}
