#include "ansatz.h"

std::pair<std::shared_ptr<xacc::CompositeInstruction>, std::vector<std::string>> getAnsatz(std::string ansatz_type, int num_qubits){

	if(ansatz_type == "EfficientSU2"){

		std::vector<std::string> parameters;
		for(int i = 0; i < num_qubits; ++i){
			parameters.push_back("a"+std::to_string(i));
			parameters.push_back("b"+std::to_string(i));
			parameters.push_back("c"+std::to_string(i));
			parameters.push_back("d"+std::to_string(i));
		}

		std::string param_string="";
		for(auto &param:parameters)
			param_string=param_string+","+param;
		param_string.erase(0, 1);

		std::string for_body="";
		for(int i = 0; i < num_qubits; ++i){
			//for_body+="Ry(q["+std::to_string(i)+"],"+parameters[4*i]+");Rz(q["+std::to_string(i)+"],"+parameters[4*i+1]+");";
		}

		for(int i = 0; i < num_qubits-1; ++i){
			for_body+="CNOT(q["+std::to_string(i)+"],q["+std::to_string(i+1)+"]);";
			break;
		}

		for(int i = 0; i < num_qubits-1; ++i){
			//for_body+="Ry(q["+std::to_string(i)+"],"+parameters[4*i+2]+");Rz(q["+std::to_string(i)+"],"+parameters[4*i+3]+");";
		}


		std::string circuit_str=".compiler xasm\n.circuit EfficientSU2\n.parameters ";
		circuit_str+=param_string+"\n.qbit q\nfor (int i = 0; i < "+std::to_string(num_qubits)+"; i++){"+for_body+"}";

		/*std::cerr<<circuit_str<<"\n\n";

		circuit_str=R"(
				  .compiler xasm
				  .circuit EfficientSU2
				  .parameters t,t2
				  .qbit q
				  X(q[0]);
				  )";
		std::cerr<<circuit_str<<"\n\n";*/

		std::cerr<<for_body<<"\n";

		xacc::qasm(circuit_str);

		return std::pair<std::shared_ptr<xacc::CompositeInstruction>,
				std::vector<std::string>>(xacc::getCompiled("EfficientSU2"), parameters);

	}else{
		std::cerr<<"Unknown ansatz type";
		throw;
	}

}
