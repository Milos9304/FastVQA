/*
 * lattice.h
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICE_H_
#define SRC_LATTICE_H_

#include "abstractLatticeInput.hpp"

#include "../logger.h"
//#include <eigen3/Eigen/Core>
#include "../symbolic_manipulation.h"
#include <vector>

class MapOptions{

	public:

		bool verbose;

		enum x_init_mode { x_symmetric };
		enum bin_mapping { naive_overapprox };
		enum penalty_mode { penalty_all };

		x_init_mode x_mode;
		bin_mapping bin_map;
		penalty_mode pen_mode;

		int penalty;
		int num_qbits_per_x;

		MapOptions(x_init_mode x_mode=x_symmetric,
				bin_mapping bin_map=naive_overapprox,
				penalty_mode pen_mode=penalty_all,
				int penalty_val=1000,
				int num_qbits_per_x=1,
				bool verbose=false){
			this->x_mode = x_mode;
			this->bin_map = bin_map;
			this->pen_mode = pen_mode;
			this->penalty = penalty_val;
			this->num_qbits_per_x = num_qbits_per_x;
			this->verbose = verbose;
		}
};

class Lattice : public AbstractLatticeInput{

	public:

	MatrixInt* lll_transformation;


		Lattice(MatrixInt lattice, std::string name = ""){ // @suppress("Class members should be properly initialized")

			this -> n = lattice.get_rows();

			this -> name = name;
			this -> orig_lattice = lattice;
			this -> orig_lattice_transposed = MatrixInt(lattice);
			this -> orig_lattice_transposed.transpose();
			this -> current_lattice = MatrixInt(lattice);
			this -> orig_gh_sq = calculate_gh_squared(&orig_lattice);

			if(lattice.get_rows()/*.rows()*/ != lattice.get_cols()/*.cols()*/){
				loge("Non-square lattice not supported");
				return;
			}

     		this -> expression_int = new Expression("expression_int");

		}

		//decode qubo optimal config to x config
		VectorInt quboToXvector(std::string measurement);

		mpq_class get_orig_gh(){return orig_gh_sq;}

		MatrixInt* get_orig_lattice(){ return &orig_lattice; }
		MatrixInt* get_orig_lattice_transposed(){ return &orig_lattice_transposed; }
		MatrixInt* get_current_lattice(){ return &current_lattice; }

		std::pair<std::vector<double>, std::vector<int>> getHmlInQuestFormulation();

		std::string toHamiltonianString();
		std::string toHamiltonianString(MapOptions* options);


		//xacc::quantum::PauliOperator getHamiltonian(MapOptions* options);

	private:

		mpq_class orig_gh_sq; //gaussian heuristics
		MatrixInt orig_lattice, orig_lattice_transposed, current_lattice;

		bool gso_current_initialized = false, gso_orig_initialized = false;
		MatGSO<Z_NR<mpz_t>, FP_NR<double>>* gso_current, *gso_orig;

		Expression *expression_int, *expression_bin, *expression_penalized, *expression_qubo;

		//std::map<std::string, Var*> qubo_to_bin_map;
		std::map<int, int> qbit_to_varId_map;

		bool x_initialized = false;
		void init_x(MapOptions::x_init_mode mode, int num_qbits_per_x, bool print=false);

		std::vector<int> x_ids;
		std::map<int, std::map<int, mpq_class>> int_to_bin_map;
		bool bin_initialized = false;
		void init_expr_bin(MapOptions::bin_mapping mapping, bool print=false);

		int z0_id=0, z1_id=0, x1_id;
		bool pen_initialized = false;
		void penalize_expr(int penalty, MapOptions::penalty_mode mode, bool print=false);

		bool qubo_generated = false;
		void generate_qubo(bool print=false);

		//xacc::quantum::PauliOperator hamiltonian;
		void calcHamiltonian(MapOptions* options, bool print);

		mpq_class calculate_gh_squared(MatrixInt* lattice);

};



#endif /* SRC_LATTICE_H_ */
