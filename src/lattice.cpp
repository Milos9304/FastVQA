/*
 * lattice.cpp
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#include "lattice.h"

bool pen_initialized = false;

void Lattice::generate_qubo(bool print){

	expression_qubo = new Expression(*expression_penalized);
	expression_qubo->name = "QUBO";

	int qubit = 0;
	for(auto &var : expression_qubo->getVariables()){

		if(var->id < 0) //id
			continue;

		std::pair<int, std::string> z = expression_qubo -> addZ(qubit++);
		qubo_to_bin_map.emplace(z.second, var);

		std::map<int, double> subs_expr; //id, coeff

		//(1-z)/2
		subs_expr.emplace(-1, 0.5);
		subs_expr.emplace(z.first, -0.5);
		expression_qubo->substitute(var->id, subs_expr);

	}

	if(print)
		expression_qubo -> print();

}

void Lattice::penalize_expr(int penalty, penalty_mode mode, bool print){

	expression_penalized = new Expression(*expression_bin);
	expression_penalized->name = "expression_penalized";

	if(mode == penalty_all){

		int z1_id=0, z2_id=0;
		std::vector<Var*>::iterator x2_it;

		expression_penalized -> addConstant(penalty);
		std::vector<Var*> variables = expression_penalized->getVariables();

		if(variables[0]->id != -1){
			loge("Error! id not the first val");
			return;
		}

		int counter = 0;
		for(std::vector<Var*>::iterator it = variables.begin() + 1;
				it != variables.end(); ++it){

			int z_id = expression_penalized -> addBinaryVar("z_"+(*it)->name);

			if(counter == 0)
				z1_id = z_id;
			else if(counter == 1){
				z2_id = z_id;
				x2_it = it;
			}

			expression_penalized -> addNewTerm((*it)->id, z_id, -penalty);

			for(std::vector<Var*>::iterator it2 = it+1; it2 != variables.end(); it2++){
				expression_penalized -> addNewTerm((*it2)->id, z_id, penalty);
			}

		counter++;
		}

		//add z0=1, z1=x1
		expression_penalized->substituteVarToDouble(z1_id, 1);
		std::map<int, double> subs_expr; //id, coeff
		subs_expr.emplace((*x2_it)->id, 1);
		expression_penalized->substitute(z2_id, subs_expr);

		//std::cout << "subs " << z1_id << " c" << 1 << "\n";
		//std::cout << "subs " << z2_id << " " << (*x2_it)->id << "\n";
	}

	if(print)
		expression_penalized->print();

}


void Lattice::init_x(x_init_mode mode, bool print){

	if(mode == x_zero_one){
		for(int i = 0; i < n; ++i){
			expression_int->addBinaryVar("x"+std::to_string(i));
		}
	}

	for(int i = 0; i < n; ++i){
		expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(i)), gram_matrix.coeff(i, i)); // G_ii*x_i^2
		for(int j = 0; j < i; ++j){
			expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(j)), 2*gram_matrix.coeff(i, j)); //2*G_ij*xi
		}
	}

	if(print)
		expression_int->print();

}

void Lattice::init_expr_bin(bin_mapping mapping, bool print){

	expression_bin = new Expression(*expression_int);
	expression_bin->name = "expression_bin";

	for(auto &var : expression_bin->getVariables()){

		if(var->id < 0) //identity coeff
			continue;

		int lb = var -> lb;
		int ub = var -> ub;

		std::string name = var -> name;

		std::map<int, double> subs_expr; //id, coeff

		if(mapping == naive_overapprox){

			subs_expr.emplace(-1, lb); //set lb to identity coeff
			for(int i = 0; i < ceil(log2(ub-lb+1)); ++i){

				int id = expression_bin->addBinaryVar(name + "_b"+std::to_string(i));
				subs_expr.emplace(id, pow(2, i));
			}

		}

		expression_bin->substitute(var->id, subs_expr);
	}

	if(print)
		expression_bin->print();
}

std::string Lattice::toHamiltonianString(x_init_mode mode, int penalty, bool print){

	if(!gram_initialized){
		gram_matrix = orig_lattice * orig_lattice.transpose();
		gram_initialized = true;
	}

	if(!x_initialized){
		init_x(mode, print);
		x_initialized = true;
	}

	if(!bin_initialized){
		init_expr_bin(naive_overapprox, print);
		bin_initialized = true;
	}

	if(!pen_initialized){
		penalize_expr(penalty, penalty_all, print);
		pen_initialized = true;
	}

	if(!qubo_generated){
		generate_qubo(print);
		qubo_generated = true;
	}

	return expression_qubo->expression_line_print();

}

