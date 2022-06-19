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
#include "internal/accelerator_base.h"

#include <vector>
#include <bitset>
#include <string>

namespace FastVQA{

typedef std::pair<qreal, long long int> RefEnergy;
typedef std::vector<RefEnergy> RefEnergies;

typedef std::function<double(double init_val, double final_val, int iter_i, int max_iters)> AlphaFunction; //initial val, final_val, and num_iterations in which alpha is being increased

struct AqcPqcAcceleratorOptions{

	int log_level=1;

	std::string accelerator_type;

};

class AqcPqcAccelerator:public AcceleratorBase{

public:

	AqcPqcAcceleratorOptions options;

	AqcPqcAccelerator(AqcPqcAcceleratorOptions options);

	void initialize(Hamiltonian* h0, Hamiltonian* h1);
	void finalize();

	double calc_intermediate_expectation(ExperimentBuffer* buffer, double lambda, bool init_zero_state=true);
	void finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples);

	std::shared_ptr<Qureg> getQuregPtr(){
		return std::make_shared<Qureg>(qureg);
	}

private:

	Hamiltonian *h0, *h1;
	int nbQubits;

	double _calc_expectation(Hamiltonian *h);
	Hamiltonian calc_intermediate_hamiltonian(double lambda);

};
}
#endif /* FASTVQA_AQC_PQC_ACCELERATOR_H_ */
