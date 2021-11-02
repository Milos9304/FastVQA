#include "accelerator.h"

void Accelerator::run(Circuit circuit, const std::vector<double> &x){

	int i = 0;
	double param_val=0;

	for(auto &gate : circuit.gates){

		if(gate.param->name != "")
			param_val = x[i++];

		apply_gate(gate, param_val);
	}

}


double Accelerator::calc_expectation(ExperimentBuffer* buffer, Ansatz* ansatz, const std::vector<double> &x){

	initZeroState(qureg);
	run(ansatz->circuit, x);
	return calcExpecDiagonalOp(qureg, hamDiag).real;


}

void Accelerator::initialize(Hamiltonian* hamIn){

	int num_qubits = hamIn->nbQubits;

	logd("Initializing " + std::to_string(num_qubits) + " qubits");

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]));

	qureg = createQureg(num_qubits, env);

	hamiltonian = createPauliHamil(num_qubits, hamIn->coeffs.size());
	hamiltonian.termCoeffs = &hamIn->coeffs[0]; //conversion to c array
	hamiltonian.pauliCodes = (enum pauliOpType*)(&hamIn->pauliOpts[0]);

	hamDiag = createDiagonalOp(num_qubits, env, 1);
	initDiagonalOpFromPauliHamil(hamDiag, hamiltonian);

}

void Accelerator::finalize(){
	logd("Destroying qureg");
	destroyQureg(qureg, env);
}

Accelerator::Accelerator(std::string accelerator){

	if(accelerator != "quest"){
		loge("No other accelerator than QuEST implemented");
		throw;
	}

}
