/*
 * AqcPqcAcelerator.h
 *
 *  Created on: June 19, 2022
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_AQC_PQC_ACCELERATOR_H_
#define FASTVQA_AQC_PQC_ACCELERATOR_H_

#include "QuEST.h"
#include "hamiltonian.h"
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

struct AqcPqcAcceleratorOptions{

	int log_level=1;

	//set to quest
	std::string accelerator_type;

	//nbSteps
	int nbSteps;

	// Name of the ansatz as defined in FastVQA/ansatz.h
	std::string ansatz_name = "Ry_CNOT_all2all_Rz";

};

class AqcPqcAccelerator:public AcceleratorBase{

public:

	AqcPqcAcceleratorOptions options;

	AqcPqcAccelerator(AqcPqcAcceleratorOptions options);

	void initialize(Hamiltonian* h0, Hamiltonian* h1);
	void run();
	void finalize();

	double calc_intermediate_expectation(ExperimentBuffer* buffer, double lambda, bool init_zero_state=true);
	void finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples);

	std::shared_ptr<Qureg> getQuregPtr(){
		return std::make_shared<Qureg>(qureg);
	}

private:

	class GeneralIntermediateHamiltonian{
	public:
		int nbQubits;

		std::vector<double> coeffs0, coeffs1;
		std::vector<int> pauliOpts;

		void toPauliHamil(PauliHamil* hamil);

	}hamil_int;

	int nbQubits;

	double _calc_expectation(Hamiltonian *h);
	Hamiltonian _calc_intermediate_hamiltonian(double lambda);

};
}
#endif /* FASTVQA_AQC_PQC_ACCELERATOR_H_ */
