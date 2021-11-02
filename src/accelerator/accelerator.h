/*
 * accelerator.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_ACCELERATOR_H_
#define SRC_ACCELERATOR_H_

#include <string>
#include "../QuEST/QuEST/include/QuEST.h"
#include "../fastVQA.h"
#include "../hamiltonian/hamiltonian.h"
#include "../ansatz/ansatz.h"

class Accelerator{

public:

	QuESTEnv env;

	Accelerator(std::string accelerator_type);
	void initialize(Hamiltonian* hamiltonian);
	void finalize();

	double calc_expectation(ExperimentBuffer* buffer, Ansatz* ansatz, const std::vector<double> &x);


private:
	Qureg qureg;

	PauliHamil hamiltonian;
    DiagonalOp hamDiag;
    void run(Circuit circuit, const std::vector<double> &x);
    void apply_gate(Gate gate, double param);

};

#endif /* SRC_ACCELERATOR_H_ */
