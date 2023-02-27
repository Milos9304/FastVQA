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

// Function to computes the Kronecker Product
// of two matrices
void Kroneckerproduct(Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> A, Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> B, Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> *C)
{
    // i loops till rowa
    for (int i = 0; i < A.rows(); i++) {
        // k loops till rowb
        for (int k = 0; k < B.rows(); k++) {
            // j loops till cola
            for (int j = 0; j < A.cols(); j++) {
                // l loops till colb
                for (int l = 0; l < B.cols(); l++) {
                    // Each element of matrix A is
                    // multiplied by whole Matrix B
                    // resp and stored as Matrix C
                    (*C)(i*B.rows() + l,j*B.cols() + k) = A(i,j) * B(k,l);
                    //cout << C(i + l + 1,j + k + 1) << " ";
                }
            }
        }
    }
}

namespace FastVQA{

std::string PauliHamiltonian::getPauliHamiltonianString(int double_precision){

	std::string res;
	bool start = true;

	if(this->nbQubits < 1 || this->pauliOpts.size() != this->nbQubits * this->coeffs.size() )
    	throw_runtime_error("PauliHamiltonian uninitialized");

	for(unsigned int i = 0; i < this->coeffs.size(); ++i){

		if(this->coeffs[i] == 0)
			continue;

		if(this->coeffs[i] > 0){
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

	this->coeffs = std::vector<qreal>{-1};

	std::vector<int> pauliOpts(this->nbQubits);
	for(int i = 0; i < this->nbQubits; ++i)
		pauliOpts[i] = 1;
	this->pauliOpts =pauliOpts;

	this->type = PauliHamiltonianType::MinusSigmaX;

}

void PauliHamiltonian::initializeSumMinusSigmaXHamiltonian(){

	this->coeffs = std::vector<qreal>(this->nbQubits);
	this->pauliOpts = std::vector<int>(this->nbQubits * this->nbQubits);

	for(int i = 0; i < this->nbQubits; ++i){
		this->coeffs[i] = -1;
		for(int j = i * this->nbQubits; j < (i+1)*this->nbQubits; ++j)
			this->pauliOpts[j] = i == (j-(i * this->nbQubits)) ? 1 : 0;
	}

	this->type = PauliHamiltonianType::SumMinusSigmaX;

}

void PauliHamiltonian::toQuestPauliHamil(PauliHamil* hamil){

	//*hamil = createPauliHamil(this->nbQubits, this->coeffs.size());

	hamil->numQubits = nbQubits;
	hamil->numSumTerms = this->coeffs.size();
	hamil->termCoeffs = &(this->coeffs)[0]; //conversion to c array
	hamil->pauliCodes = (enum pauliOpType*)(&this->pauliOpts[0]);

}

//WIP
Eigen::MatrixXd PauliHamiltonian::getMatrixRepresentation(bool diagonalOp){

	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(1ULL<<(this->nbQubits), 1ULL<<(this->nbQubits));

	if(diagonalOp){
		PauliHamil hamil;
		QuESTEnv env = createQuESTEnv();
		this->toQuestPauliHamil(&hamil);
		DiagonalOp hamDiag = createDiagonalOp(this->nbQubits, env, 1);
		initDiagonalOpFromPauliHamil(hamDiag, hamil);
		destroyQuESTEnv(env);

		for(int i = 0; i < 1ULL<<(this->nbQubits); ++i){
			m(i,i) = hamDiag.real[i];
			//std::cerr<<hamDiag.real[i]<<". ";
		}

		return m;
	}
	throw_runtime_error("Not implemented");
	/*if(this->nbQubits < 1 || this->pauliOpts.size() != this->nbQubits * this->coeffs.size() )
	    	throw_runtime_error("PauliHamiltonian uninitialized");

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

	}*/
}

Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> PauliHamiltonian::getMatrixRepresentation2(bool diagonalOp){

	Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> m = Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>::Zero(1ULL<<(this->nbQubits), 1ULL<<(this->nbQubits));

	if(diagonalOp){
		PauliHamil hamil;
		QuESTEnv env = createQuESTEnv();
		this->toQuestPauliHamil(&hamil);
		DiagonalOp hamDiag = createDiagonalOp(this->nbQubits, env, 1);
		initDiagonalOpFromPauliHamil(hamDiag, hamil);


		for(int i = 0; i < 1ULL<<(this->nbQubits); ++i){
			m(i,i) = hamDiag.real[i];
			//std::cerr<<hamDiag.real[i]<<". ";
		}

		destroyDiagonalOp(hamDiag, env);
		destroyQuESTEnv(env);
		return m;
	}


	if(this->nbQubits < 1 || this->pauliOpts.size() != this->nbQubits * this->coeffs.size() )
	    	throw_runtime_error("PauliHamiltonian uninitialized");

	Eigen::Matrix<qreal, 2, 2> PAULI_I {{1,0},{0,1}};
	Eigen::Matrix<qreal, 2, 2> PAULI_X {{0,1},{1,0}};
	//Eigen::Matrix<qreal, 2, 2> PAULI_Y {{1,0},{0,1}};
	Eigen::Matrix<qreal, 2, 2> PAULI_Z {{1,0},{0,-1}};

	for(int i = 0; i < this->coeffs.size(); ++i){

		if(this->coeffs[i] == 0)
			continue;

		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> m1(1, 1);
		m1(0,0) = 1;

		for(int q = 0; q < this->nbQubits; ++q){
			int opt = this->pauliOpts[i*this->nbQubits + q];

			Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> m2 = m1;
			m1.resize(m2.rows()*2, m2.cols()*2);

			if(opt == 0)
				Kroneckerproduct(PAULI_I, m2, &m1);
			else if(opt == 1)
				Kroneckerproduct(PAULI_X, m2, &m1);
			else if(opt == 2)
				throw_runtime_error("PauliY unimplemented");
			else if(opt == 3)
				Kroneckerproduct(PAULI_Z, m2, &m1);
			else
				throw_runtime_error("Invalid PauliHamiltonian");
		}

		m += this->coeffs[i] * m1;

		/*bool all_x=true;
		for(int q = 0; q < this->nbQubits; ++q){
			int opt = this->pauliOpts[i*this->nbQubits + q];
			if(opt != 1){
				all_x = false;
			}
		}
		if(all_x){
			throw_runtime_error("Not implemented");
		}*/
	}

	return m;
}

}
