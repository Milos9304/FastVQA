#include <iostream>
#include <iterator>
#include <exception>

#include "popl.hpp"
#include "logger.h"
#include "qaoa/qaoa.h"
#include "vqaConfig.h"

using namespace popl;

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
		for(auto &lattice : vqaConfig->getLattices()){

			logi("Running lattice " + std::to_string(i++) + " / " + std::to_string(num_lattices));
			lattice.toHamiltonianString(Lattice::x_zero_one, true);break;
			//run_qaoa(lattice.toHamiltonianString(Lattice::x_zero_one), vqaConfig->verbose);

		}

		return 0;
    }

    else{
    	std::cout << "Invalid argument. Use ./fastVQA -h to see help." << "\n";
    	return 1;
    }

    return 0;
}
