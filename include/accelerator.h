/*
 * accelerator.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_ACCELERATOR_H_
#define FASTVQA_ACCELERATOR_H_

#include "QuEST.h"
#include "hamiltonian.h"
#include "ansatz.h"
#include "experimentBuffer.h"

#include <vector>
#include <functional>
#include <bitset>
#include <string>

namespace fastVQA{

typedef std::pair<qreal, long long int> RefEnergy;
typedef std::vector<RefEnergy> RefEnergies;

typedef std::function<double(double init_val, double final_val, int iter_i, int max_iters)> AlphaFunction; //initial val, final_val, and num_iterations in which alpha is being increased

struct AcceleratorOptions{

	int log_level;

	std::string accelerator_type;
	std::string alpha_f; //constant, linear

	double samples_cut_ratio=1; //1=100% means no cutting is performed
	std::vector<long long unsigned int> zero_reference_states;

	double final_alpha=1;
	int max_alpha_iters=1000;

};

class Accelerator{

public:

	QuESTEnv env;
	Ansatz ansatz;

	AlphaFunction alpha_f = alpha_constant_f;

	static AlphaFunction alpha_constant_f;
	static AlphaFunction alpha_linear_f;

	AcceleratorOptions options;

	Accelerator(std::shared_ptr<AcceleratorOptions> options);
	void initialize(Hamiltonian* hamiltonian);
	void finalize();

	double calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int iteration_i, double* ground_state_overlap_out);
	void finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples);

	void set_ansatz(Ansatz* ansatz);

	void run_vqe_slave_process();

private:

	Qureg qureg;
	int *qubits_list;
	pauliOpType* all_x_list;

	PauliHamil hamiltonian;
	RefEnergies ref_hamil_energies;

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
}
#endif /* SRC_ACCELERATOR_H_ */
