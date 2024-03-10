/*
 * accelerator.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_ACCELERATOR_H_
#define FASTVQA_ACCELERATOR_H_

#include "QuEST.h"
#include "pauliHamiltonian.h"
#include "cost_function.h"
#include "experimentBuffer.h"
#include "accelerator_base.h"

#include <vector>
#include <bitset>
#include <string>

namespace FastVQA{

typedef std::pair<qreal, long long int> RefEnergy;
typedef std::vector<RefEnergy> RefEnergies;

typedef std::function<double(double init_val, double final_val, int iter_i, int max_iters)> AlphaFunction; //initial val, final_val, and num_iterations in which alpha is being increased

struct AcceleratorOptions{

	int log_level=1;

	std::string accelerator_type;

	// Set "linear" for linear CVar. "constant" is default behaviour
	std::string alpha_f = "constant";

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

class Accelerator:public AcceleratorBase{

public:

	AlphaFunction alpha_f = alpha_constant_f;

	static AlphaFunction alpha_constant_f;
	static AlphaFunction alpha_linear_f;

	AcceleratorOptions options;

	Accelerator(AcceleratorOptions options);

	void initialize(CostFunction cost_function, int num_qubits);
	void initialize(PauliHamiltonian* hamIn);
	void initialize(int num_qubits);
	void finalize();

	double calc_expectation(ExperimentBuffer* buffer);
	double calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int iteration_i, double* ground_state_overlap_out);
	void finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples);

	void run_vqe_slave_process();

	std::shared_ptr<Qureg> getQuregPtr(){
		return std::make_shared<Qureg>(qureg);
	}

	RefEnergies getEigenspace();
	RefEnergies getSolutions();

private:

	int log_level;

	int *qubits_list;
	pauliOpType* all_x_list;

	//false - PauliHamiltonian is used
	//true - cost function is used
	bool hamiltonian_specified;

	CostFunction cost_function;

	bool pauliHamilInitialized = false;
	PauliHamil pauliHamiltonian;
	bool ref_hamil_energies_set = false;
	RefEnergies ref_hamil_energies;

    DiagonalOp hamDiag;

    void __initialize(int num_qubits);

    void run_with_new_params(Circuit circuit, const std::vector<double> &x);

    double _energy_evaluation(double* ground_state_overlap_out, int iteration_i);
    double evaluate_assignment(PauliHamil isingHam, std::string measurement);

    const int control_tag = 42;
    	const int control_val_exit=1;
    	//const int control_reset_all=2;
    	const int control_setting_new_seed=3;

    //const int new_overlap_penalty_tag = 16;
    const int new_qureg_tag = 17;
    const int new_circuit_tag = 18;
    const int new_params_tag = 19;
    const int new_PauliHamiltonian_tag = 20;
    const int calc_exp_val_tag = 21;
    const int measure_with_cache_tag = 22;

    //const int data_tag=35;

};
}
#endif /* FASTVQA_ACCELERATOR_H_ */
