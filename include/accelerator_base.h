/*
 * accelerator_base.h
 *
 *  Created on: June 19, 2022
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_ACCELERATOR_BASE_H_
#define FASTVQA_ACCELERATOR_BASE_H_

#include "QuEST.h"
#include "gate.h"
#include "ansatz.h"

#include <vector>
#include <bitset>
#include <string>

namespace FastVQA{

//typedef std::pair<qreal, long long int> RefEnergy;
struct RefEnergy{
	qreal value; //.first
	long long int index; //.second
	bool isConsideredSolution;

	RefEnergy(qreal value, long long int index, bool isConsideredSolution){
		this->value = value;
		this->index = index;
		this->isConsideredSolution = isConsideredSolution;
	}
};

typedef std::vector<RefEnergy> RefEnergies;

class AcceleratorBase{

public:
	QuESTEnv env;
	Ansatz ansatz{"", 1};

	void set_ansatz(Ansatz *ansatz);
	void run_circuit(Circuit circuit, bool init_zero_state=true);

protected:

	Qureg qureg;

	void apply_gate(Gate gate);
	void apply_gate(Gate gate, qreal param);
};

}
#endif /* FASTVQA_ACCELERATOR_BASE_H_ */
