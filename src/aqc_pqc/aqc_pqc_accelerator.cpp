#include "aqc_pqc/AqcPqcAccelerator.h"
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
	this->env = createQuESTEnv();

}

void AqcPqcAccelerator::initialize(Hamiltonian* h0, Hamiltonian* h1){

	this->nbQubits = h0->nbQubits;
	if(this->nbQubits != h1->nbQubits)
		throw_runtime_error("Dimensions of initial and final Hamiltonian do not match");

	hamil_int.coeffs0 = h0->coeffs;
	hamil_int.pauliOpts = h0->pauliOpts;

	std::vector<std::string> sequences;

	for(int i = 0; i < hamil_int.coeffs0.size(); ++i){
		std::string s;

		for(int j = 0; j < nbQubits; ++j)
			s+=std::to_string(h0->pauliOpts[nbQubits*i+j]);

		std::vector<std::string>::iterator itr = std::find(sequences.begin(), sequences.end(), s);
		if (itr != sequences.cend())
			throw_runtime_error("Hamiltonian contains multiple equivalent terms");
		sequences.push_back(s);
		hamil_int.coeffs1.push_back(0);
	}

	for(int i = 0; i < h1->coeffs.size(); ++i){
		std::string s;

		for(int j = 0; j < nbQubits; ++j)
			s+=std::to_string(h1->pauliOpts[nbQubits*i+j]);

		std::vector<std::string>::iterator itr = std::find(sequences.begin(), sequences.end(), s);
		if (itr != sequences.cend()){
			int index = std::distance(sequences.begin(), itr);
			hamil_int.coeffs1[index]+=h1->coeffs[i];
		}
		else{
			sequences.push_back(s);
			for(auto &c:s)
				hamil_int.pauliOpts.push_back(c-'0');
			hamil_int.coeffs0.push_back(0);
			hamil_int.coeffs1.push_back(h1->coeffs[i]);
		}
	}
	if(hamil_int.coeffs0.size() != hamil_int.coeffs1.size())
		throw_runtime_error("Problem with mixer hamiltonian calculation");

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]), options.log_level);
	this->qureg = createQureg(this->nbQubits, env);

}

void AqcPqcAccelerator::run(){

	int nbSteps = options.nbSteps;
	logi(std::to_string(nbSteps) + " steps");

	Ansatz ansatz = getAnsatz(options.ansatz_name, hamil_int.nbQubits);

}

Hamiltonian AqcPqcAccelerator::_calc_intermediate_hamiltonian(double lambda){

	std::vector<double> coeffs;

	for(int i = 0; i < hamil_int.coeffs0.size(); ++i){
		coeffs.push_back((1-lambda)*hamil_int.coeffs0[i]+lambda*hamil_int.coeffs1[i]);
	}

	Hamiltonian res(nbQubits, coeffs, hamil_int.pauliOpts);

	return res;
}

double AqcPqcAccelerator::calc_intermediate_expectation(ExperimentBuffer* buffer, double lambda, bool init_zero_state){

	if(lambda < 0 || lambda > 1)
		throw_runtime_error("Invalid lambda value");

	Hamiltonian h = _calc_intermediate_hamiltonian(lambda);
	return _calc_expectation(&h);

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
