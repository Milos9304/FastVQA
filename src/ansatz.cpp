#include "ansatz.h"
#include "logger.h"

const double pi = 3.141592653589793238463;

namespace FastVQA{

Ansatz getAnsatz(std::string ansatz_type, int num_qubits, int depth, int seed){

	Ansatz ansatz(ansatz_type, depth);
    std::mt19937 gen(seed); //rd() instead of seed
	std::uniform_real_distribution<> dis(-2.0, 2.0);

	if(depth < 1)
		throw_runtime_error("Incorrect ansatz depth value");

	if(ansatz_type == "EfficientSU2"){

		for(int d = 0; d < depth; ++d){
			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				double param2 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i)+"_"+std::to_string(d), param1)));
				ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i)+"_"+std::to_string(d), param2)));
			}

			for(int i = 0; i < num_qubits-1; ++i){
				ansatz.circuit.addGate(Gate::g_CNOT, i, i+1);
			}
			ansatz.circuit.addGate(Gate::g_CNOT, num_qubits-1, 0);


			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				double param2 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("c"+std::to_string(i)+"_"+std::to_string(d), param1)));
				ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("d"+std::to_string(i)+"_"+std::to_string(d), param2)));
			}
		}

		ansatz.num_params = num_qubits * 4 * depth;

	}else if(ansatz_type == "Ry_CNOT_all2all_Ry"){

		for(int d = 0; d < depth; ++d){

			for(int i = 0; i < num_qubits; ++i){
					double param1 = dis(gen);
					ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

			for(int i = 0; i < num_qubits-1; ++i){
				for(int j = i+1; j < num_qubits; ++j)
					ansatz.circuit.addGate(Gate::g_CNOT, i, j);
			}
			ansatz.circuit.addGate(Gate::g_CNOT, num_qubits-1, 0);

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}
		}

		ansatz.num_params = num_qubits * 2 * depth;

	}


	else if(ansatz_type == "Ry_CNOT_all2all_Rz"){

		for(int d = 0; d < depth; ++d){

			for(int i = 0; i < num_qubits; ++i){
					double param1 = dis(gen);
					ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

			for(int i = 0; i < num_qubits-1; ++i){
				for(int j = i+1; j < num_qubits; ++j)
					ansatz.circuit.addGate(Gate::g_CNOT, i, j);
			}
			ansatz.circuit.addGate(Gate::g_CNOT, num_qubits-1, 0);

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

		}

		ansatz.num_params = num_qubits * 2 * depth;

	}else if(ansatz_type == "Ry_CNOT_nn_Rz"){

		for(int d = 0; d < depth; ++d){

			for(int i = 0; i < num_qubits; ++i){
					double param1 = dis(gen);
					ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

			for(int i = 0; i < num_qubits-1; ++i)
					ansatz.circuit.addGate(Gate::g_CNOT, i, i+1);
			ansatz.circuit.addGate(Gate::g_CNOT, num_qubits-1, 0);

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

		}

		ansatz.num_params = num_qubits * 2 * depth;
	}else if(ansatz_type == "Ry_Cz_all2all_Ry"){

		for(int d = 0; d < depth; ++d){

			for(int i = 0; i < num_qubits; ++i){
					double param1 = dis(gen);
					ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

			for(int i = 0; i < num_qubits-1; ++i){
				for(int j = i+1; j < num_qubits; ++j)
					ansatz.circuit.addGate(Gate::g_CZ, i, j);
			}
			ansatz.circuit.addGate(Gate::g_CZ, num_qubits-1, 0);

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}
		}

		ansatz.num_params = num_qubits * 2 * depth;

	}else if(ansatz_type == "Ry_Cz_nn_Ry"){

		for(int d = 0; d < depth; ++d){

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

			for(int i = 0; i < num_qubits-1; ++i)
				ansatz.circuit.addGate(Gate::g_CZ, i, i+1);
			ansatz.circuit.addGate(Gate::g_CZ, num_qubits-1, 0);

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}
		}

		ansatz.num_params = num_qubits * 2 * depth;
	}else if(ansatz_type == "Ry_CNOT_nn_Rz_CNOT_Rz"){

		for(int d = 0; d < depth; ++d){

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Ry, i, std::shared_ptr<Parameter>(new Parameter("a"+std::to_string(i)+"_"+std::to_string(d), param1)));
			}

			for(int i = 0; i < num_qubits-1; ++i)
				ansatz.circuit.addGate(Gate::g_CNOT, i, i+1);
			ansatz.circuit.addGate(Gate::g_CNOT, num_qubits-1, 0);

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("b"+std::to_string(i), param1)));
			}

			for(int i = 0; i < num_qubits-1; ++i)
				ansatz.circuit.addGate(Gate::g_CNOT, i, i+1);
			ansatz.circuit.addGate(Gate::g_CNOT, num_qubits-1, 0);

			for(int i = 0; i < num_qubits; ++i){
				double param1 = dis(gen);
				ansatz.circuit.addParametrizedGate(Gate::g_Rz, i, std::shared_ptr<Parameter>(new Parameter("c"+std::to_string(i), param1)));
			}
		}

		ansatz.num_params = num_qubits * 3 * depth;
	}

	else if(ansatz_type == "qaoa"){

		int p = depth;
		ansatz.circuit.qaoa_ansatz = true;
		ansatz.num_params = p * 2;

	}else{
		throw_runtime_error("Unknown ansatz type");
	}
	return ansatz;
}
void initOptimalParamsForMinusSigmaXHamiltonian(Ansatz *ansatz){

	std::vector<std::shared_ptr<Parameter>> params = ansatz->circuit.getParamsPtrs();

	if(ansatz->name == "Ry_CNOT_all2all_Rz" || ansatz->name == "Ry_CNOT_nn_Rz"){

		int i = 0;
		for(auto &param_ptr:params){
			if(param_ptr->name == Parameter::blank_param_name)
				continue;

			if(i++ < ansatz->num_params/(ansatz->depth*2)){
				param_ptr->value = pi/2;
			}else
				param_ptr->value = 0;


		}

	}else if(ansatz->name == "Ry_Cz_all2all_Ry" || ansatz->name == "Ry_Cz_nn_Ry"){

		int i = 0;
		for(auto &param_ptr:params){
			if(param_ptr->name == Parameter::blank_param_name)
				continue;

			if(i++ < ansatz->num_params-(ansatz->num_params/(ansatz->depth*2))){
				param_ptr->value = 0;
			}else
				param_ptr->value = pi/2;
		}

	}else if(ansatz->name == "Ry_CNOT_nn_Rz_CNOT_Rz"){

		int i = 0;
		for(auto &param_ptr:params){
			if(param_ptr->name == Parameter::blank_param_name)
				continue;

			if(i++ < ansatz->num_params/(ansatz->depth*3)){
				param_ptr->value = pi/2;
			}else
				param_ptr->value = 0;
		}
	}else{
		throw_runtime_error("Optimal -sum_sigma_x parameter selection not implemented for this ansatz");
	}
}
}
