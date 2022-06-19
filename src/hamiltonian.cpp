#include "hamiltonian.h"

namespace FastVQA{

void Hamiltonian::initializeMinusSigmaXHamiltonian(){

	this->coeffs = std::vector<double>{-1};

	std::vector<int> pauliOpts(this->nbQubits);
	for(int i = 0; i < this->nbQubits; ++i)
		pauliOpts[i] = 1;
	this->pauliOpts =pauliOpts;

}

void Hamiltonian::toPauliHamil(PauliHamil* hamil){

	*hamil = createPauliHamil(this->nbQubits, this->coeffs.size());

	hamil->termCoeffs = &(this->coeffs)[0]; //conversion to c array
	hamil->pauliCodes = (enum pauliOpType*)(&this->pauliOpts[0]);
}

}
