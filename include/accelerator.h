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
#include "cost_function.h"
#include "experimentBuffer.h"

#include <vector>
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

	bool exclude_zero_state = false;
	bool choose_ground_state_with_smallest_index = false;

	// Set to false if running repetitive experiments with the same number of qubits. Need to create qureg manually in this case.
	// If false, need to specify number of qubits.
	bool createQuregAtEachInilization = true;
	int createQuregAtEachInilization_num_qubits = -1;

};

class Accelerator{

public:

	QuESTEnv env;
	Ansatz ansatz;

	AlphaFunction alpha_f = alpha_constant_f;

	static AlphaFunction alpha_constant_f;
	static AlphaFunction alpha_linear_f;

	AcceleratorOptions options;

	Accelerator(AcceleratorOptions options);

	void initialize(CostFunction cost_function, int num_qubits);
	void initialize(Hamiltonian* hamiltonian);
	void finalize();

	double calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int iteration_i, double* ground_state_overlap_out);
	void finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples);

	void set_ansatz(Ansatz* ansatz);

	void run_vqe_slave_process();

	std::shared_ptr<Qureg> getQuregPtr(){
		return std::make_shared<Qureg>(qureg);
	}

private:

	Qureg qureg;
	int *qubits_list;
	pauliOpType* all_x_list;

	//false - hamiltonian is used
	//true - cost function is used
	bool hamiltonian_specified;

	CostFunction cost_function;

	PauliHamil hamiltonian;
	RefEnergies ref_hamil_energies;

    DiagonalOp hamDiag;

    void __initialize(int num_qubits);

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
