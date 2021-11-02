#include "accelerator.h"

void Accelerator::apply_gate(Gate gate, double param){

	switch(gate.code){
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
	default:
		loge("Unknown gate");
		throw;

		break;
	}

}
