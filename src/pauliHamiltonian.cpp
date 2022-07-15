#include "pauliHamiltonian.h"
#include "logger.h"

#include <sstream>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

namespace FastVQA{

std::string PauliHamiltonian::getPauliHamiltonianString(int double_precision){

	std::string res;
	bool start = true;

	if(this->nbQubits < 1 || this->pauliOpts.size() != this->nbQubits * this->coeffs.size() )
    	throw_runtime_error("PauliHamiltonian unitialized");

	for(int i = 0; i < this->coeffs.size(); ++i){

		if(this->coeffs[i] == 0)
			continue;

		if(this->coeffs[i] > 0){//std::cerr<<"jebo"<<this->coeffs[i]<<"   ";
			if(!start)
				res+="+ ";
		}
		start = false;

		res+=to_string_with_precision(this->coeffs[i], double_precision)+" ";
		for(int q = 0; q < this->nbQubits; ++q){

			switch(this->pauliOpts[i*this->nbQubits + q]){
				case 0:
					break;
				case 1:
					res+="X"+std::to_string(q)+" ";
					break;
				case 2:
					res+="Y"+std::to_string(q)+" ";
					break;
				case 3:
					res+="Z"+std::to_string(q)+" ";
					break;
				default:
					throw_runtime_error("Invalid PauliOpts");
					break;
			}
		}
	}

	return res;
}

void PauliHamiltonian::initializeMinusSigmaXHamiltonian(){

	this->coeffs = std::vector<double>{-1};

	std::vector<int> pauliOpts(this->nbQubits);
	for(int i = 0; i < this->nbQubits; ++i)
		pauliOpts[i] = 1;
	this->pauliOpts =pauliOpts;

}

void PauliHamiltonian::toQuestPauliHamil(PauliHamil* hamil){

	*hamil = createPauliHamil(this->nbQubits, this->coeffs.size());

	hamil->termCoeffs = &(this->coeffs)[0]; //conversion to c array
	hamil->pauliCodes = (enum pauliOpType*)(&this->pauliOpts[0]);
}

//WIP
/*Eigen::MatrixXd PauliHamiltonian::getMatrixRepresentation(bool diagonalOp){


	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(this->nbQubits, this->nbQubits);

	if(diagonalOp){
		logw("Not implemented");
		return m;
	}

	if(this->nbQubits < 1 || this->pauliOpts.size() != this->nbQubits * this->coeffs.size() )
	    	throw_runtime_error("PauliHamiltonian unitialized");

	for(int i = 0; i < this->coeffs.size(); ++i){

		if(this->coeffs[i] == 0)
			continue;

		bool all_x=true;
		for(int q = 0; q < this->nbQubits; ++q){

			int opt = this->pauliOpts[i*this->nbQubits + q];
			if(opt != 1){
				all_x = false;

			}

		}

	}
}*/

}
