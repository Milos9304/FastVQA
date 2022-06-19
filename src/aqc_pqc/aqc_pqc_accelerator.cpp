#include "AqcPqcAccelerator.h"
#include "logger.h"
#include <iomanip>
#include <fstream>
#include <algorithm>

namespace FastVQA{

AqcPqcAccelerator::AqcPqcAccelerator(AqcPqcAcceleratorOptions options){

	if(options.accelerator_type != "quest"){
		loge("No other accelerator than QuEST implemented");
		throw;
	}

	this->options = options;
}

void AqcPqcAccelerator::initialize(Hamiltonian* h0, Hamiltonian* h1){

	this->nbQubits = h0->nbQubits;
	if(this->nbQubits != h1->nbQubits)
		throw_runtime_error("Dimensions of initial and final Hamiltonian do not match");

	this->h0 = h0;
	this->h1 = h1;

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]), options.log_level);

	this->qureg = createQureg(this->nbQubits, env);

}

Hamiltonian AqcPqcAccelerator::calc_intermediate_hamiltonian(double lambda){

}

double AqcPqcAccelerator::calc_intermediate_expectation(ExperimentBuffer* buffer, double lambda, bool init_zero_state){

	if(lambda == 0){
		return _calc_expectation(h0);
	}else{
		throw_runtime_error("TODO: UNIMPLEMENTED");
	}

}

double AqcPqcAccelerator::_calc_expectation(Hamiltonian *h){
	initZeroState(qureg);
	run_circuit(ansatz.circuit);
	Qureg workspace = createQureg(nbQubits, env);
	PauliHamil hamil;
	h->toPauliHamil(&hamil);
	return calcExpecPauliHamil(qureg, hamil, workspace);
}

}
