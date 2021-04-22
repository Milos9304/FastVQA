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
#include "lattice/hmlLattice.hpp"
#include <xacc.hpp>

#include "lattiQ.h"

#include "xacc_observable.hpp" //del
#include "PauliOperator.hpp" //del

using namespace popl;

int main(int ac, char** av){

	OptionParser op("Allowed options");
	auto help_option  = op.add<Switch>("h", "help", "produce help message");
    auto qaoa 		  = op.add<Switch>("", "qaoa", "run qaoa algorithm");
    auto config 	  = op.add<Value<std::string>>("", "config", "config file location", "");
    auto lattice_file = op.add<Value<std::string>>("l", "lattice", "lattice file location", "");
    auto niters       = op.add<Value<int>>("i", "iters", "max num of iterations", 0);
    auto save_hml     = op.add<Value<std::string>>("", "savehml", "save hamiltonian to file", "");
    auto load_hml     = op.add<Value<std::string>>("", "loadhml", "save hamiltonian to file", "");
    auto debug        = op.add<Switch>("d", "debug", "print debug messages");

    auto save_interm  = op.add<Value<std::string>>("", "si", "save intermediate results (for specific experiments only)", "");
    auto load_interm  = op.add<Value<std::string>>("", "li", "load intermediate results (for specific experiments only)", "");


    op.parse(ac, av);

    if (help_option->is_set()){
    	std::cout << op << "\n";
    	return 0;
    }

	int num_lattices;
	bool hml_lattice_mode=false;
	HmlLattice* hmlLat;
	Lattice* loadL;
	std::vector<Lattice> lattices_in;
	std::vector<AbstractLatticeInput*> lattices;

    if(qaoa -> is_set()){

    	if(!load_hml->is_set()){

			VqaConfig* vqaConfig;

			if(!config->is_set()){

				if(!lattice_file->is_set()){
					loge("Neither config nor lattice file specified");
					return 1;
				}

				bool success;
				MatrixInt m = VqaConfig::loadLatticeFromFile(lattice_file->value(), &success);
				loadL = new Lattice(m, "");
				lattices.push_back(loadL);
				num_lattices = 1;

				if(success)
					logi(lattice_file->value() + " loaded");
				else
					loge("Problem loading " + lattice_file->value());

			}else{

				//config file load
				try{
					vqaConfig = new VqaConfig(config->value());
				}
				catch(std::exception &e){
					loge(e.what());
					return 1;
				}

				logi("Running QAOA. Configuration loaded from " + config->value());
				lattices_in = vqaConfig->getLattices();
				for(auto &l : lattices_in)
					lattices.push_back(&l);

				num_lattices = lattices.size();

				logi(std::to_string(num_lattices) + " lattice(s) in the config file");

			}

    	}else{ //load from hml file

    		hml_lattice_mode = true;

    		std::ifstream ifs(load_hml->value(), std::ios::binary | std::ios::in);
    		if(!ifs.is_open()){
    			loge("Cannot open " + load_hml->value());
    			return 1;
    		}

    		int n;
    		std::string n_str, hamiltonian;
    		std::getline(ifs, n_str);
    		std::getline(ifs, hamiltonian);
    		std::istringstream iss(n_str);
    		iss >> n;
    		ifs.close();

    		num_lattices = 1;
    		hmlLat = new HmlLattice(n, hamiltonian);
    		hmlLat->name = load_hml->value();
    	}


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
		qaoaOptions.max_iters = niters->is_set() ? (niters->value() == 0 ? 5000 : niters->value()): 5000;
		qaoaOptions.detailed_log_freq = 50;
		qaoaOptions.verbose = debug->is_set();
		qaoaOptions.debug = debug->is_set();
		qaoaOptions.verbose |= true;
		qaoaOptions.optimizer = optimizer;
		qaoaOptions.accelerator = accelerator;
		qaoaOptions.simplifiedSimulation = true;
		qaoaOptions.logEnergies = true;
		qaoaOptions.extendedParametrizedMode = false;//true;
		qaoaOptions.calcVarAssignment = true;
		qaoaOptions.provideHamiltonian = true;
		qaoaOptions.saveIntermediate = save_interm->is_set() ? (save_interm->value() == "" ? false : true) : false;
		qaoaOptions.s_intermediateName = qaoaOptions.saveIntermediate ? save_interm->value() : "";
		qaoaOptions.loadIntermediate = load_interm->is_set() ? (load_interm->value() == "" ? false : true) : false;
		qaoaOptions.l_intermediateName = qaoaOptions.loadIntermediate ? load_interm->value() : "";

		MapOptions* mapOptions = new MapOptions();
		mapOptions->verbose = false;

		ExecutionStatistics* execStats = new ExecutionStatistics();
		xacc::Initialize();
		xacc::setOption("quest-verbose", "true");
		xacc::setOption("quest-debug", "true");

		if(hml_lattice_mode){
			xacc::qbit* buffer;
			ProgressBar bar{bar_opts(0, 1, hmlLat->name)};
			qaoaOptions.set_default_stats_function(execStats, &bar, hmlLat);
    		Qaoa::run_qaoa(&buffer, hmlLat->toHamiltonianString(), load_hml->value(), &bar, execStats, &qaoaOptions);
    		logd("Hml mode");
		}else{
			int counter = 0;
			for(auto &lattice_abs : lattices){

				logi("Running " + lattice_abs->name);

				Lattice *lattice = static_cast<Lattice*>(lattice_abs);

				lattice->toHamiltonianString(mapOptions); //remake
				std::pair<std::vector<double>, std::vector<int>> hamiltonian2 = lattice->getHmlInQuestFormulation();

				ProgressBar bar{bar_opts(counter, num_lattices, lattice->name)};

				qaoaOptions.set_default_stats_function(execStats, &bar, lattice);

				QOracle quantum_oracle = [&bar, execStats, &qaoaOptions, hamiltonian2]
										  (xacc::qbit** buffer, std::string hamiltonian, std::string name) {
					Qaoa::run_qaoa(buffer, hamiltonian, hamiltonian2, name, &bar, execStats, &qaoaOptions);
				};

				IterativeLatticeReduction ilr(lattice, mapOptions, quantum_oracle, 1);

				if(save_hml->is_set() && save_hml->value() != ""){
					std::ofstream ofs(save_hml->value(), std::ios::binary | std::ios::out);
					ofs << lattice->n << "\n";
					ofs << lattice->toHamiltonianString(mapOptions);
					ofs.close();
				}

				ilr.run();

				counter++;
			}
		}

		xacc::Finalize();

		return 0;
    }

    else{
    	std::cout << "Invalid argument. Use ./lattiQ -h to see help." << "\n";
    	return 1;
    }

    return 0;
}
