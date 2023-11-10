#include "vqas/qaoa.h"
#include "logger.h"

#include <iostream>
//#include <xacc.hpp>

#include "QuEST.h"
#include <random>

#include <fstream>
#include "mpi.h"

namespace FastVQA{

void Qaoa::run_qaoa(ExperimentBuffer* buffer, PauliHamiltonian* hamiltonian, QAOAOptions* options){

	this->num_qubits = hamiltonian->nbQubits;
	this->__initialize(buffer, options);
	options->accelerator->initialize(hamiltonian);

	this->ansatz = getAnsatz("qaoa", this->num_qubits, this->p, 0);
	this->num_params = ansatz.num_params;
	options->accelerator->set_ansatz(&ansatz);

	logd("QAOA starting", log_level);
	__execute(buffer, options->accelerator, options->optimizer);
	logd("QAOA execution done", this->log_level);

}

void Qaoa::__initialize(ExperimentBuffer* buffer, QAOAOptions* options){
	this->instance_name = options->instance_name;
	this->max_iters = options->max_iters;
	this->ftol = options->ftol;
	this->nbSamples_calcVarAssignment = options->nbSamples_calcVarAssignment;
    this->log_level = options->log_level;
    this->p = options->p;
}

void Qaoa::__execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt){

	std::string instance_prefix = "[["+this->instance_name+"]] ";

	int iteration_i = 0;
	OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {
			double ground_state_overlap;
			double expectation = acc->calc_expectation(buffer, x, iteration_i++, &ground_state_overlap);
			return expectation;
		}, num_params);

	std::vector<double> initial_params;
	std::mt19937 gen(1997); //rd() instead of 1997
	std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
	for(int i = 0; i < num_params/2; ++i){
		initial_params.push_back(dis(gen));
		initial_params.push_back(dis(gen));
	}
	std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
	std::vector<double> upperBounds(initial_params.size(), 3.141592654);

	logd("QAOA starting optimization", this->log_level);
	OptResult result = opt->optimize(f, initial_params, this->ftol, this->max_iters, lowerBounds, upperBounds);
	logd("QAOA finishing optimization", this->log_level);

	std::string opt_config;
	acc->finalConfigEvaluator(buffer, result.second, nbSamples_calcVarAssignment);
	if(log_level <= 1){
		logi(instance_prefix + "Final opt-val: " + std::to_string(buffer->opt_val));
		for(auto &solution : buffer->final_solutions){
			logi(instance_prefix + "Final opt-config: " + solution.opt_config);
			logi(instance_prefix + "Final hit-rate: " + std::to_string(solution.hit_rate));
		}
	}

}

/*void Qaoa::run_qaoa_slave_process(){
	//auto acc = xacc::getAccelerator("quest");
	loge("NOT IMPLEMENTED");
	throw;
}*/
}
