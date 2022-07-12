#include "accelerator_base.h"
#include "logger.h"

namespace FastVQA{

void AcceleratorBase::run_circuit(Circuit circuit, bool init_zero_state){

	if(init_zero_state)
		initZeroState(qureg);

	for(auto &gate : circuit.gates)
		apply_gate(gate);
}

void AcceleratorBase::set_ansatz(Ansatz *ansatz){

	this -> ansatz = *ansatz;

	std::vector<double> gateCodes;
	for(auto &gate : this->ansatz.circuit.gates){
		gateCodes.push_back(gate.code);
		gateCodes.push_back(gate.qubit1);
		gateCodes.push_back(gate.qubit2);
	}
}

void AcceleratorBase::apply_gate(Gate gate){

	double param = gate.param->value;
	apply_gate(gate, param);

}

void AcceleratorBase::apply_gate(Gate gate, double param){

	switch(gate.code){
	case Gate::g_H:
		logd("H " + std::to_string(gate.qubit1));
		hadamard(qureg, gate.qubit1);
		break;
	case Gate::g_Ry:
		logd("Ry " + std::to_string(gate.qubit1) + " @ " + std::to_string(param));
		rotateY(qureg, gate.qubit1, param);
		break;
	case Gate::g_Rz:
		logd("Rz " + std::to_string(gate.qubit1) + " @ " + std::to_string(param));
		rotateZ(qureg, gate.qubit1, param);
		break;
	case Gate::g_CNOT:
		logd("CNOT " + std::to_string(gate.qubit1) + " " + std::to_string(gate.qubit2));
		controlledNot(qureg, gate.qubit1, gate.qubit2);
		break;
	case Gate::g_CZ:
		logd("CZ " + std::to_string(gate.qubit1) + " " + std::to_string(gate.qubit2));
		controlledPhaseFlip(qureg, gate.qubit1, gate.qubit2);
		break;
	default:
		loge("Unknown gate");
		throw;

		break;
	}

}
}
