/*
 * Circuit.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_CIRCUIT_H_
#define SRC_CIRCUIT_H_

#include "gate.h"
#include <vector>
namespace FastVQA{

typedef std::vector<std::string> Parameters;

class Circuit{

public:

	bool qaoa_ansatz = false;
	bool cm_qaoa_ansatz = false;

	Circuit(){}

	void addGate(GateCode code, int qubit1);
	void addGate(GateCode code, int qubit1, int qubit2);
	void addParametrizedGate(GateCode code, int qubit1, std::shared_ptr<Parameter> param);
	void addParametrizedGate(GateCode code, int qubit1, int qubit2, std::shared_ptr<Parameter> param);

	std::vector<std::shared_ptr<Parameter>> getParamsPtrs();
	std::vector<Gate> gates;

};
}
#endif /* SRC_CIRCUIT_H_ */
