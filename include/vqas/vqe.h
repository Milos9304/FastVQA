/*
 * vqe.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef VQAS_VQE_H_
#define VQAS_VQE_H_

#include "ansatz.h"
#include "cost_function.h"
#include "vqeOptions.h"

#include <functional>

namespace fastVQA{

class Vqe{

	public:

		void run_vqe(ExperimentBuffer* buffer, CostFunction cost_f, int num_qubits, VQEOptions* options);
		void run_vqe(ExperimentBuffer* buffer, Hamiltonian* hamiltonian, VQEOptions* options);

	private:

		Ansatz ansatz;

		int num_qubits;
		int num_params;
		int max_iters;
		int nbSamples_calcVarAssignment;

		std::string instance_name;

		int log_level;

		void __initialize(ExperimentBuffer* buffer, VQEOptions* options);

		void execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, std::vector<long long unsigned int> zero_reference_states, CostFunction cost_f, bool logExpecStd=false);
		void execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, std::vector<long long unsigned int> zero_reference_states, Hamiltonian* hamiltonian, bool logExpecStd=false);
		void __execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, bool logExpecStd);

};
}

#endif /* VQAS_VQE_H_ */
