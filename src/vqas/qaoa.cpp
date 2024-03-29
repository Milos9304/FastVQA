#include "vqas/qaoa.h"
#include "logger.h"

#include <iostream>
//#include <xacc.hpp>

#include "QuEST.h"
#include <random>

#include <fstream>
//#include "mpi.h"

std::string nlopt_res_to_string(int result){
  switch(result)
  {
    case -1: return "FAILURE";
    case -2: return "INVALID_ARGS";
    case -3: return "OUT_OF_MEMORY";
    case -4: return "ROUNDOFF_LIMITED";
    case -5: return "FORCED_STOP";
    case 1: return "SUCCESS";
    case 2: return "STOPVAL_REACHED";
    case 3: return "FTOL_REACHED";
    case 4: return "XTOL_REACHED";
    case 5: return "MAXEVAL_REACHED";
    case 6: return "MAXTIME_REACHED";
    default: return NULL;
  }
}

namespace FastVQA{

void Qaoa::run_qaoa(ExperimentBuffer* buffer, PauliHamiltonian* hamiltonian, QAOAOptions* options){

	this->num_qubits = hamiltonian->nbQubits;
	this->__initialize(buffer, options);
	options->accelerator->initialize(hamiltonian);

	if(buffer->storeQuregPtr){
		buffer->stateVector = options->accelerator->getQuregPtr();
	}

	this->ansatz = getAnsatz("qaoa", this->num_qubits, this->p, 0);
	this->num_params = ansatz.num_params;
	options->accelerator->set_ansatz(&ansatz);

	logd("QAOA starting", log_level);
	__execute(buffer, options->accelerator, options->optimizer);
	logd("QAOA execution done", this->log_level);
}

void Qaoa::run_cm_qaoa(ExperimentBuffer* buffer, PauliHamiltonian* hamiltonian, QAOAOptions* options, long long int zero_index){

	this->num_qubits = hamiltonian->nbQubits;
	this->__initialize(buffer, options);
	options->accelerator->initialize(hamiltonian);

	if(buffer->storeQuregPtr){
		buffer->stateVector = options->accelerator->getQuregPtr();
	}

	this->ansatz = getAnsatz("cm_qaoa", this->num_qubits, this->p, 0);
	this->num_params = ansatz.num_params;
	options->accelerator->set_ansatz(&ansatz);
	options->accelerator->zero_index = zero_index;

	logd("CM_QAOA starting", log_level);
	__execute(buffer, options->accelerator, options->optimizer);
	logd("CM_QAOA execution done", this->log_level);
}

void Qaoa::run_qaoa_fixed_angles(ExperimentBuffer* buffer, PauliHamiltonian* hamiltonian, QAOAOptions* options, const double *angles){

	this->num_qubits = hamiltonian->nbQubits;
	this->__initialize(buffer, options);
	options->accelerator->initialize(hamiltonian, true);
	//return;
	if(buffer->storeQuregPtr){
		buffer->stateVector = options->accelerator->getQuregPtr();
	}

	this->ansatz = getAnsatz("qaoa", this->num_qubits, this->p, 0);
	this->num_params = ansatz.num_params;
	options->accelerator->set_ansatz(&ansatz);

	logd("QAOA starting", log_level);
	double ground_state_overlap;
	std::vector<double> x;
	x.assign(angles, angles+2*this->p);
	/*double expectation = */options->accelerator->calc_expectation(buffer, x, 0, &ground_state_overlap);
	//std::cerr<<buffer->stateVector->stateVec.real[0];
	logd("QAOA execution done", this->log_level);

}


void Qaoa::__initialize(ExperimentBuffer* buffer, QAOAOptions* options){
	this->instance_name = options->instance_name;
	this->max_iters = options->max_iters;
	this->ftol = options->ftol;
	this->nbSamples_calcVarAssignment = options->nbSamples_calcVarAssignment;
    this->log_level = options->log_level;
    this->p = options->p;
    this->options_ptr = options;
}

void Qaoa::__execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt){

	std::string instance_prefix = "[["+this->instance_name+"]] ";

	int iteration_i = 0;
	std::vector<std::vector<double>> intermediateAngles;
	std::vector<double> intermediateEnergies;

	OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {
			double ground_state_overlap;
			double expectation = acc->calc_expectation(buffer, x, iteration_i++, &ground_state_overlap);
			intermediateAngles.push_back(x);
			intermediateEnergies.push_back(expectation);
			return expectation;
		}, num_params);

	std::vector<double> initial_params;
	if(this->options_ptr->initial_params.size() == 0){
		std::mt19937 gen(0); //rd() instead of 0 - seed
		std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
		for(int i = 0; i < num_params/2; ++i){
			double param1 = dis(gen);
			double param2 = dis(gen);
			initial_params.push_back(param1);
			initial_params.push_back(param2);
			buffer->initial_params.push_back(std::pair<std::string, double>("alpha_"+std::to_string(i),param1));
			buffer->initial_params.push_back(std::pair<std::string, double>("beta_"+std::to_string(i),param2));
		}
	}else{
		for(int i = 0; i < num_params/2; ++i){
			double param1 = this->options_ptr->initial_params[2*i];
			double param2 = this->options_ptr->initial_params[2*i+1];
			initial_params.push_back(param1);
			initial_params.push_back(param2);
			buffer->initial_params.push_back(std::pair<std::string, double>("alpha_"+std::to_string(i),param1));
			buffer->initial_params.push_back(std::pair<std::string, double>("beta_"+std::to_string(i),param2));
		}
	}

	std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
	std::vector<double> upperBounds(initial_params.size(), 3.141592654);

	logd("QAOA starting optimization", this->log_level);
	OptResult result = opt->optimize(f, initial_params, this->ftol, this->max_iters, lowerBounds, upperBounds);
	logd("QAOA finishing optimization", this->log_level);

	std::string opt_config;
	acc->finalConfigEvaluator(buffer, result.first.second, nbSamples_calcVarAssignment);
	if(log_level <= 1){
		logi(instance_prefix + "Final opt-val: " + std::to_string(buffer->opt_val));
		for(auto &solution : buffer->final_solutions){
			logi(instance_prefix + "Final opt-config: " + solution.opt_config);
			logi(instance_prefix + "Final hit-rate: " + std::to_string(solution.hit_rate));
		}
	}

	buffer->num_iters = iteration_i;
	buffer->finalParams = result.first.second;
	buffer->opt_message = nlopt_res_to_string(result.second);
	buffer->intermediateAngles = intermediateAngles;
	buffer->intermediateEnergies = intermediateEnergies;

}

/*void Qaoa::run_qaoa_slave_process(){
	//auto acc = xacc::getAccelerator("quest");
	loge("NOT IMPLEMENTED");
	throw;
}*/
}
