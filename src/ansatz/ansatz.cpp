#include "ansatz.h"

Ansatz getAnsatz(std::string ansatz_type, int num_qubits, int seed){

	Ansatz ansatz;
    std::mt19937 gen(seed); //rd() instead of seed
	std::uniform_real_distribution<> dis(-2.0, 2.0);

	if(ansatz_type == "EfficientSU2"){

		for(int i = 0; i < num_qubits; ++i){
			double param1 = dis(gen);
			double param2 = dis(gen);
			ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i), param1)));
			ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i), param2)));
		}

		for(int i = 0; i < num_qubits-1; ++i){
			ansatz.circuit.addGate(Gate::g_CNOT, i, i+1);
		}

		for(int i = 0; i < num_qubits; ++i){
			double param1 = dis(gen);
			double param2 = dis(gen);
			ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("c"+std::to_string(i), param1)));
			ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("d"+std::to_string(i), param2)));
		}

		ansatz.num_params = num_qubits * 4;

	}else if(ansatz_type == "Ry_CNOT_all2all_Ry"){

		for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i), param1)));
		}

		for(int i = 0; i < num_qubits-1; ++i){
			for(int j = i+1; j < num_qubits-1; ++j)
				ansatz.circuit.addGate(Gate::g_CZ, i, j);
		}

		for(int i = 0; i < num_qubits; ++i){
			double param1 = dis(gen);
			ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("n"+std::to_string(i), param1)));
		}

		ansatz.num_params = num_qubits * 2;

	}else if(ansatz_type == "qaoa"){

		int p = 1;
		ansatz.circuit.qaoa_ansatz = true;
		ansatz.num_params = p * 2;

	}else{
		std::cerr<<"Unknown ansatz type";
		throw;
	}
	return ansatz;
}
