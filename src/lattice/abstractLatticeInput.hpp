/*
 * abstractLatticeInput.hpp
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICE_ABSTRACTLATTICEINPUT_HPP_
#define SRC_LATTICE_ABSTRACTLATTICEINPUT_HPP_

#include <string>

//typedef Eigen::Matrix<mpz_class, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixInt;
#include <gmpxx.h>
#include "fplll.h"

typedef ZZ_mat<mpz_t> MatrixInt;
typedef std::vector<mpz_class> VectorInt;

class AbstractLatticeInput{

	public:
		std::string name;
		int n;

		//decode qubo optimal config to x config
		virtual VectorInt quboToXvector(std::string measurement) = 0;

		virtual std::string toHamiltonianString() = 0;
		virtual ~AbstractLatticeInput(){}

};



#endif /* SRC_LATTICE_ABSTRACTLATTICEINPUT_HPP_ */
