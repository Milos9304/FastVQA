/*
 * AqcPqcAcelerator.h
 *
 *  Created on: June 19, 2022
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_AQC_PQC_ACCELERATOR_H_
#define FASTVQA_AQC_PQC_ACCELERATOR_H_

#include "QuEST.h"
#include "pauliHamiltonian.h"
#include "cost_function.h"
#include "experimentBuffer.h"
#include "accelerator_base.h"
#include <Eigen/Dense>

#include <vector>
#include <bitset>
#include <string>
#include <fstream>

namespace FastVQA{

const long double PI =   3.14159265358979323846264338327950288419716939937510L;
const long double PI_2 = 1.57079632679489661923132169163975144209858469968755L;

typedef std::pair<qreal, long long int> RefEnergy;
typedef std::vector<RefEnergy> RefEnergies;

typedef std::function<double(double init_val, double final_val, int iter_i, int max_iters)> AlphaFunction; //initial val, final_val, and num_iterations in which alpha is being increased

enum InitialGroundState {None, PlusState};

struct AqcPqcAcceleratorOptions{

	int log_level=1;

	//set to quest
	std::string accelerator_type;

	//nbSteps
	int nbSteps;

	bool checkHessian = true;

	// Name of the ansatz as defined in FastVQA/ansatz.h
	std::string ansatz_name = "Ry_CNOT_all2all_Rz";

	bool compareWithClassicalEigenSolver = false;

	bool outputLogToFile = false;
	std::string logFileName = "log";

	bool printGroundStateOverlap = false;
	InitialGroundState initialGroundState = None;

	bool checkSolutions = false;
	qreal solutionExpectation;
	std::vector<long long int> solutions;

	//-1 avoids rounding
	int roundDecimalPlaces=-1;

	//print epsilons
	bool printEpsilons = false;

	int optStrategy = 0;

	double xtol = 10e-5;

};

class AqcPqcAccelerator:public AcceleratorBase{

public:

	AqcPqcAcceleratorOptions options;

	AqcPqcAccelerator(AqcPqcAcceleratorOptions options);
	~AqcPqcAccelerator();

	void initialize(PauliHamiltonian* h0, PauliHamiltonian* h1);
	void run();

	qreal calc_intermediate_expectation(ExperimentBuffer* buffer, double lambda, bool init_zero_state=true);
	void finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<qreal> final_params, int nbSamples);

	std::shared_ptr<Qureg> getQuregPtr(){
		return std::make_shared<Qureg>(qureg);
	}

	qreal _calc_expectation(PauliHamiltonian *h);

private:

	bool initialized = false;
	void finalize();

	Eigen::Vector<qreal, Eigen::Dynamic> _optimize_with_rank_reduction(PauliHamiltonian *h, Eigen::Vector<qreal, Eigen::Dynamic> *minus_q, Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> *A, std::vector<std::shared_ptr<Parameter>> *parameters);
	Eigen::Vector<qreal, Eigen::Dynamic> _optimize_trivially(PauliHamiltonian *h, Eigen::Vector<qreal, Eigen::Dynamic> *minus_q, Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> *A, std::vector<std::shared_ptr<Parameter>> *parameters);

	class GeneralIntermediatePauliHamiltonian{
	public:
		int nbQubits;
		PauliHamiltonianType initial_type = PauliHamiltonianType::General;

		std::vector<qreal> coeffs0, coeffs1;
		std::vector<int> pauliOpts;

		PauliHamiltonian h1;

		void toPauliHamil(PauliHamil* hamil);

		Qureg getIntermediateGroundState(double lambda);

	}hamil_int;

	int nbQubits;

	Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> getIntermediateMatrixRepresentation(PauliHamiltonian* h, double* id_coeff, double theta);

	PauliHamiltonian _calc_intermediate_hamiltonian(double lambda);

	Qureg workspace;

	std::ofstream logFile;

	qreal __round(qreal x);

};

}
#endif /* FASTVQA_AQC_PQC_ACCELERATOR_H_ */
