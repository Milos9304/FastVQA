/*
 * accelerator.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_ACCELERATOR_H_
#define SRC_ACCELERATOR_H_

#include <mpi.h>
#include <string>
#include "../QuEST/QuEST/include/QuEST.h"
#include "../fastVQA.h"
#include "../hamiltonian/hamiltonian.h"
#include "../ansatz/ansatz.h"

class Accelerator{

public:

	QuESTEnv env;
	Ansatz ansatz;

	Accelerator(std::string accelerator_type);
	void initialize(Hamiltonian* hamiltonian);
	void finalize();

	double calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int zero_reference_state);
	void finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples);

	void set_ansatz(Ansatz* ansatz);

	void run_vqe_slave_process();

private:
	Qureg qureg;

	PauliHamil hamiltonian;
    DiagonalOp hamDiag;
    void run(Circuit circuit, const std::vector<double> &x);

    void apply_gate(Gate gate, double param);

    double evaluate_assignment(PauliHamil isingHam, std::string measurement);

    const int control_tag = 42;
    	const int control_val_exit=1;
    	//const int control_reset_all=2;
    	const int control_setting_new_seed=3;

    //const int new_overlap_penalty_tag = 16;
    const int new_qureg_tag = 17;
    const int new_circuit_tag = 18;
    const int new_params_tag = 19;
    const int new_hamiltonian_tag = 20;
    const int calc_exp_val_tag = 21;
    const int measure_with_cache_tag = 22;

    //const int data_tag=35;

};

#endif /* SRC_ACCELERATOR_H_ */
