#include <iostream>
#include <iterator>
#include <exception>

#include "popl.hpp"
#include "config_io.h"
#include "logger.h"
#include "qaoa/qaoa.h"

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
		int num_hamiltonians = vqaConfig->getHamiltonians().size();

		logi(std::to_string(num_hamiltonians) + " hamiltonian(s) in the config file");

		int i = 0;
		for(auto &hamiltonian : vqaConfig->getHamiltonians()){

			logi("Running hamiltonian " + std::to_string(i++) + " / " + std::to_string(num_hamiltonians));

			vqaConfig->setCurrentHamiltonian(hamiltonian);
			run_qaoa(vqaConfig);

		}

		return 0;
    }

    else{
    	std::cout << "Invalid argument. Use ./fastVQA -h to see help." << "\n";
    	return 1;
    }

    return 0;
}
