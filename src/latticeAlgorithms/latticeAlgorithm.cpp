/*
 * latticeAlgorithm.cpp
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#include "latticeAlgorithm.h"

void LatticeAlgorithm::performLLLonLattice(double delta, double eta, LLLMethod method){

	//lll_reduction(*(lattice->get_orig_lattice()), delta, eta, method, floatType, precision, flags)

	lll_reduction(*(lattice->get_current_lattice()), delta, eta, method, FloatType::FT_DOUBLE);

}

VectorInt LatticeAlgorithm::xVectToShortVect(VectorInt* x_vect){

	VectorInt short_lattice_vector;
	MatrixInt* orig_lattice_t = lattice->get_orig_lattice_transposed();

	for(int i = 0; i < orig_lattice_t->get_cols(); ++i){
		mpz_class val = 0;
		NumVect<Z_NR<mpz_t>> row = orig_lattice_t->matrix[i]; //row because transposed
		for(unsigned int j = 0; j < x_vect->size(); ++j){
			val += (*x_vect)[j] * mpz_class(row[j].get_data());
		}
		short_lattice_vector.push_back(val);
	}

	return short_lattice_vector;

}
