#include "accelerator.h"
#include "logger.h"

namespace fastVQA{

void Accelerator::apply_gate(Gate gate, double param){

	std::string message;

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
