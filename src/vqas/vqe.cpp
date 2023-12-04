#include "vqas/vqe.h"
#include "logger.h"

#include <iostream>
//#include <xacc.hpp>

//#include "QuEST.h"
#include <random>

#include <fstream>


//#include "mpi.h"
namespace FastVQA{
void Vqe::run_vqe(ExperimentBuffer* buffer,
		PauliHamiltonian* hamiltonian,
		VQEOptions* options){

	this->num_qubits = hamiltonian->nbQubits;
	this->__initialize(buffer, options);

	logd("Ansatz generated.", log_level);
	execute(buffer, options->accelerator, options->optimizer, options->zero_reference_states, hamiltonian, options->expectationToStandardOutput, options->keepReferenceToQureg);
	logd("Vqe execution done", this->log_level);
	logi("Min QUBO found: " + std::to_string(buffer->opt_val), this->log_level);

}

void Vqe::run_vqe(ExperimentBuffer* buffer,
		CostFunction cost_function,
		int num_qubits,
		VQEOptions* options){

	this->num_qubits = num_qubits;
	this->__initialize(buffer, options);

	logd("Ansatz generated.", log_level);
	execute(buffer, options->accelerator, options->optimizer, options->zero_reference_states, cost_function, options->expectationToStandardOutput, options->keepReferenceToQureg);
	logd("Vqe execution done", this->log_level);
	logi("Min QUBO found: " + std::to_string(buffer->opt_val), this->log_level);

}

void Vqe::__initialize(ExperimentBuffer* buffer, VQEOptions* options){

	this->max_iters = options->max_iters;
	this->ftol = options->ftol;
	this->instance_name = options->instance_name;
	this->nbSamples_calcVarAssignment = options->nbSamples_calcVarAssignment;
	this->log_level = options->log_level;

	logd("Running VQE with " + std::to_string(this->num_qubits) + " qubits", log_level);

	this->ansatz = getAnsatz(options->ansatz_name, this->num_qubits, 1997);
	this->num_params = ansatz.num_params;

}


void Vqe::execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* optimizer, std::vector<long long unsigned int> zero_reference_states, PauliHamiltonian* PauliHamiltonian, bool logExpecStd, bool keepQureg){
	acc->options.zero_reference_states = zero_reference_states;
	acc->initialize(PauliHamiltonian);
	__execute(buffer, acc, optimizer, logExpecStd, keepQureg);
}

void Vqe::execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* optimizer, std::vector<long long unsigned int> zero_reference_states, CostFunction cost_function, bool logExpecStd, bool keepQureg){
	acc->options.zero_reference_states = zero_reference_states;
	acc->initialize(cost_function, this->num_qubits);
	logd("Num qubits set to " + std::to_string(this->num_qubits));
	__execute(buffer, acc, optimizer, logExpecStd, keepQureg);
}

void Vqe::__execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* optimizer, bool logExpecStd, bool keepQureg){

	std::string instance_prefix = "[["+instance_name+"]] ";

	std::vector<double> intermediateEnergies;
	acc->set_ansatz(&ansatz);

	if(acc->alpha_f(0,1,1,1)==0){
		if(log_level <= 2)
			logw("f: Constant Alpha: " + std::to_string(acc->options.samples_cut_ratio));
	}else{
		if(log_level <= 2)
			logw("f: Linear Init alpha: " + std::to_string(acc->options.samples_cut_ratio) + " Final alpha: " +	std::to_string(acc->options.final_alpha) + " Max iters: " + std::to_string(acc->options.max_alpha_iters));
	}
	int iteration_i = 0;

	OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {
		double ground_state_overlap;
		double expectation = acc->calc_expectation(buffer, x, iteration_i++, &ground_state_overlap);
		buffer->intermediateEnergies.push_back(expectation);
		buffer->intermediateGroundStateOverlaps.push_back(ground_state_overlap);
		if(logExpecStd)
			std::cout << instance_prefix << "e=" << std::to_string(expectation) << std::endl;
		return (expectation);

	}, num_params);

	std::vector<double> initial_params;

	OptResult result;

	if(ansatz.circuit.qaoa_ansatz){

		//std::random_device rd;
		std::mt19937 gen(1997); //rd() instead of 1997
		std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);

		for(int i = 0; i < ansatz.num_params/2; ++i){
			initial_params.push_back(dis(gen));
			initial_params.push_back(dis(gen));
		}

		std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
		std::vector<double> upperBounds(initial_params.size(), 3.141592654);

		result = optimizer->optimize(f, initial_params, this->ftol, this->max_iters, lowerBounds, upperBounds);

	}else{

		for(auto &gate : ansatz.circuit.gates)
			if(gate.param->name != "")
				initial_params.push_back(gate.param->value);

		std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
		std::vector<double> upperBounds(initial_params.size(),  3.141592654);

		//logw("are these really good ones?");
		result = optimizer->optimize(f, initial_params, this->ftol, this->max_iters, lowerBounds, upperBounds);
	}

	//double finalCost = result.first;
	std::string opt_config;

	acc->finalConfigEvaluator(buffer, result.second, nbSamples_calcVarAssignment);
	if(log_level <= 1){
		logi(instance_prefix + "Final opt-val: " + std::to_string(buffer->opt_val));
		for(auto &solution : buffer->final_solutions){
			logi(instance_prefix + "Final opt-config: " + solution.opt_config);
			logi(instance_prefix + "Final hit-rate: " + std::to_string(solution.hit_rate));
		}
	}

	if(keepQureg)
		buffer->stateVector = acc->getQuregPtr();
	else
		acc->finalize();

}
}
