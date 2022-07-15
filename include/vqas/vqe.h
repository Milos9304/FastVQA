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

namespace FastVQA{

class Vqe{

	public:

		void run_vqe(ExperimentBuffer* buffer, CostFunction cost_f, int num_qubits, VQEOptions* options);
		void run_vqe(ExperimentBuffer* buffer, PauliHamiltonian* hamiltonian, VQEOptions* options);

	private:

		Ansatz ansatz;

		int num_qubits;
		int num_params;
		long long int max_iters;
		double ftol;
		int nbSamples_calcVarAssignment;

		std::string instance_name;

		int log_level;

		void __initialize(ExperimentBuffer* buffer, VQEOptions* options);

		void execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, std::vector<long long unsigned int> zero_reference_states, CostFunction cost_f, bool logExpecStd=false, bool keepReferenceToQureg = false);
		void execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, std::vector<long long unsigned int> zero_reference_states, PauliHamiltonian* hamiltonian, bool logExpecStd=false, bool keepReferenceToQureg = false);
		void __execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, bool logExpecStd, bool keepReferenceToQureg);

};
}

#endif /* VQAS_VQE_H_ */
