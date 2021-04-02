/*
 * lattice.h
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICE_H_
#define SRC_LATTICE_H_

#include "logger.h"
#include <eigen3/Eigen/Core>
#include "symbolic_manipulation.h"

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixInt;

class Lattice{

	public:

		enum x_init_mode { x_zero_one };
		enum bin_mapping { naive_overapprox };
		enum penalty_mode { penalty_all };

		Lattice(MatrixInt lattice){ // @suppress("Class members should be properly initialized")

			this -> orig_lattice = lattice;

			if(lattice.rows() != lattice.cols()){
				loge("Non-square lattice not supported");
				return;
			}

			this -> n = lattice.rows();
			this -> expression_int = new Expression("expression_int");
		}

		std::string toHamiltonianString(x_init_mode mode);

	private:

		int n;
		MatrixInt orig_lattice;

		bool gram_initialized = false;
		MatrixInt gram_matrix;

		Expression *expression_int, *expression_bin, *expression_penalized, *expression_qubo;

		std::map<std::string, Var*> qubo_to_bin_map;

		bool x_initialized = false;
		void init_x(x_init_mode mode);

		bool bin_initialized = false;
		void init_expr_bin(bin_mapping mapping);

		bool pen_initialized = false;
		void penalize_expr(int penalty, penalty_mode mode);

		bool qubo_generated = false;
		void generate_qubo();

};



#endif /* SRC_LATTICE_H_ */
