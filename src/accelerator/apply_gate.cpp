#include "accelerator.h"

inline void logdd(std::string s){
	if(false)
		logd(s);
}

void Accelerator::apply_gate(Gate gate, double param){

	std::string message;

	switch(gate.code){
	case Gate::g_Ry:
		logdd("Ry " + std::to_string(gate.qubit1) + " @ " + std::to_string(param));
		rotateY(qureg, gate.qubit1, param);
		break;
	case Gate::g_Rz:
		logdd("Rz " + std::to_string(gate.qubit1) + " @ " + std::to_string(param));
		rotateZ(qureg, gate.qubit1, param);
		break;
	case Gate::g_CNOT:
		logdd("CNOT " + std::to_string(gate.qubit1) + " " + std::to_string(gate.qubit2));
		controlledNot(qureg, gate.qubit1, gate.qubit2);
		break;
	default:
		loge("Unknown gate");
		throw;

		break;
	}

}
