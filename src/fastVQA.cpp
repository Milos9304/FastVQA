#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
#include <exception>

#include "vqe/qaoa.h"

#include "config_io.h"
#include "logger.h"


namespace po = boost::program_options;

int main(int ac, char** av){

	std::string config;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
		("qaoa", "run qaoa algorithm")
        ("config,c", po::value<std::string>(&config), "set config file location");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    if (vm.count("qaoa")) {
    	if(config.empty()){
    		loge("No config file specified");
    		return 1;
    	}

    	try{
    		VqeConfig vqeConfig = loadConfigFile(config);
    	}catch(std::exception &e){
    		loge(e.what());
    		return 1;
    	}

    	logi("Running QAOA. Configuration loaded from " + config);
    	run_dummy_qaoa();


    	return 0;
    }

    loge("Invalid argument settings");
    return 1;
}
