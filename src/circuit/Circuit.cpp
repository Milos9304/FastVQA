#include "Circuit.h"

void Circuit::addGate(GateCode code, int qubit1){
	gates.push_back(Gate(code, qubit1));
}

void Circuit::addGate(GateCode code, int qubit1, int qubit2){
	gates.push_back(Gate(code, qubit1, qubit2));

}

void Circuit::addParametrizedGate(GateCode code, int qubit1, std::shared_ptr<Parameter> param){
	gates.push_back(Gate(code, qubit1, param));
}

void Circuit::addParametrizedGate(GateCode code, int qubit1, int qubit2, std::shared_ptr<Parameter> param){
	gates.push_back(Gate(code, qubit1, qubit2, param));
}
