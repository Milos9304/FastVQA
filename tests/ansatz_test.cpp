/*
 * ansatz_test.cpp
 *
 *  Created on: June 18, 2022
 *      Author: Milos Prokop
 */

#include <gtest/gtest.h>
#include <cmath>
#include "fastVQA.h"

#define GTEST_COUT std::cerr << "[          ] [ INFO ]"

TEST(ansatz_test, initOptimalParamsForMinusSigmaXHamiltonian){

	int num_qubits = 5;
	std::vector<std::string> initial_ansatz_names = {"Ry_CNOT_all2all_Rz", "Ry_CNOT_nn_Rz", "Ry_Cz_all2all_Ry", "Ry_Cz_nn_Ry", "Ry_CNOT_nn_Rz_CNOT_Rz"};

	FastVQA::AcceleratorOptions acceleratorOptions;
	acceleratorOptions.accelerator_type = "quest";
	//acceleratorOptions.createQuregAtEachInilization = false;
	//acceleratorOptions.createQuregAtEachInilization_num_qubits = 1;

	FastVQA::Accelerator accelerator(acceleratorOptions);

	accelerator.initialize(num_qubits);

	for(auto &ansatz_name: initial_ansatz_names){
		FastVQA::Ansatz ansatz = FastVQA::getAnsatz(ansatz_name, num_qubits, 0);
		FastVQA::initOptimalParamsForMinusSigmaXHamiltonian(&ansatz);
		GTEST_COUT << "Testing " << ansatz_name << std::endl;
		accelerator.run_circuit(ansatz.circuit);

		std::shared_ptr<Qureg> qureg = accelerator.getQuregPtr();

		//Should return all plus state
		qreal plus_elem = sqrt(1./qureg->numAmpsTotal);
		for(long long i = 0; i < qureg->numAmpsTotal; ++i){
			ASSERT_DOUBLE_EQ(plus_elem, qureg->stateVec.real[i]);
			ASSERT_DOUBLE_EQ(0., qureg->stateVec.imag[i]);
		}
	}
}

TEST(ansatz_test, initializeMinusSigmaXHamiltonian){

	int num_qubits = 5;
	std::vector<std::string> initial_ansatz_names = {"Ry_CNOT_all2all_Rz", "Ry_CNOT_nn_Rz", "Ry_Cz_all2all_Ry", "Ry_Cz_nn_Ry", "Ry_CNOT_nn_Rz_CNOT_Rz"};

	FastVQA::AqcPqcAcceleratorOptions acceleratorOptions;
	acceleratorOptions.accelerator_type = "quest";

	FastVQA::PauliHamiltonian hamiltonian(num_qubits);
	hamiltonian.initializeSumMinusSigmaXHamiltonian();

	FastVQA::AqcPqcAccelerator accelerator(acceleratorOptions);

	accelerator.initialize(&hamiltonian, &hamiltonian);

	for(auto &ansatz_name: initial_ansatz_names){
		FastVQA::Ansatz ansatz = FastVQA::getAnsatz(ansatz_name, num_qubits, 0);
		FastVQA::initOptimalParamsForMinusSigmaXHamiltonian(&ansatz);
		GTEST_COUT << "Testing " << ansatz_name << std::endl;

		FastVQA::ExperimentBuffer buffer;
		accelerator.set_ansatz(&ansatz);
		double expectation = accelerator.calc_intermediate_expectation(&buffer, 0);
		ASSERT_DOUBLE_EQ(expectation, -1*num_qubits);
	}

}

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	return ret;

}
