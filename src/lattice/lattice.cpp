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

		std::pair<int, std::string> z = expression_qubo -> addZ(qubit);
		//qubo_to_bin_map.emplace(z.second, var);

		std::map<int, mpq_class> subs_expr; //id, coeff

		//(1-z)/2
		subs_expr.emplace(-1, 0.5);
		subs_expr.emplace(z.first, -0.5);
		expression_qubo->substitute(var->id, subs_expr);

		qbit_to_varId_map.emplace(qubit, var->id);

		qubit++;

	}

	if(print)
		expression_qubo -> print();

}

void Lattice::penalize_expr(int penalty, MapOptions::penalty_mode mode, bool print){

	expression_penalized = new Expression(*expression_bin);
	expression_penalized->name = "expression_penalized";

	if(mode == MapOptions::penalty_all){

		std::vector<Var*>::iterator x1_it;

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
				z0_id = z_id;
			else if(counter == 1){
				z1_id = z_id;
				x1_it = it;
			}

			expression_penalized -> addNewTerm((*it)->id, z_id, -penalty);

			for(std::vector<Var*>::iterator it2 = it+1; it2 != variables.end(); it2++){
				expression_penalized -> addNewTerm((*it2)->id, z_id, penalty);
			}

		counter++;
		}

		//add z0=1, z1=x1
		x1_id = (*x1_it)->id;
		expression_penalized->substituteVarToDouble(z0_id, 1);
		std::map<int, mpq_class> subs_expr; //id, coeff
		subs_expr.emplace(x1_id, 1);
		expression_penalized->substitute(z1_id, subs_expr);

		//std::cout << "subs " << z1_id << " c" << 1 << "\n";
		//std::cout << "subs " << z2_id << " " << (*x2_it)->id << "\n";
	}else if(mode == MapOptions::overlap_trick){
		//do nothing as no penalty qubits needed
	}

	if(print)
		expression_penalized->print();

}


void Lattice::init_x(MapOptions::x_init_mode mode, int num_qbits_per_x, bool print){

	Z_NR<mpz_t> coeff;

	if(mode == MapOptions::x_symmetric){

		if(num_qbits_per_x == 1)
			for(int i = 0; i < n_rows; ++i){
				int id = expression_int->addBinaryVar("x"+std::to_string(i));
				x_ids.push_back(id);
			}
		else{

			int lb = -pow(2, num_qbits_per_x)/ 2 + 1;
			int ub = 1-lb;

			for(int i = 0; i < n_rows; ++i){
				int id = expression_int->addIntegerVar("x"+std::to_string(i), lb, ub);
				x_ids.push_back(id);
			}
		}
	}

	for(int i = 0; i < n_rows; ++i){
		gso_current->get_int_gram(coeff, i, i);
		expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(i)), coeff.get_data()/*.coeff(i, i)*/); // G_ii*x_i^2
		for(int j = 0; j < i; ++j){
			mpz_class c(gso_current->get_int_gram(coeff, i, j).get_data());

			//std::cout << "i" << i << " j" << j << " c"<<c << "\n";

			expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(j)), 2*c/*.coeff(i, j)*/); //2*G_ij*xi
		}
	}

	if(print)
		expression_int->print();

}

void Lattice::init_expr_bin(MapOptions::bin_mapping mapping, bool print){

	expression_bin = new Expression(*expression_int);
	expression_bin->name = "expression_bin";

	for(auto &var : expression_bin->getVariables()){

		if(var->id < 0) //identity coeff
			continue;

		int lb = var -> lb;
		int ub = var -> ub;

		std::string name = var -> name;

		std::map<int, mpq_class> subs_expr; //id, coeff

		if(mapping == MapOptions::naive_overapprox){

			subs_expr.emplace(-1, lb); //set lb to identity coeff
			for(int i = 0; i < ceil(log2(ub-lb+1)); ++i){

				int id = expression_bin->addBinaryVar(name + "_b"+std::to_string(i));
				subs_expr.emplace(id, pow(2, i));
			}

		}

		expression_bin->substitute(var->id, subs_expr);
		int_to_bin_map.emplace(var->id, subs_expr);
	}

	if(print)
		expression_bin->print();
}

void Lattice::calcHamiltonian(MapOptions* options, bool print){

		if(!gso_current_initialized){

			logw("calc gso");

			ZZ_mat<mpz_t> blank;

			gso_current = new MatGSO<Z_NR<mpz_t>, FP_NR<double>>(current_lattice, blank, blank, GSO_INT_GRAM);
			gso_current->update_gso();

			gso_current_initialized = true;
		}

		if(!x_initialized){
			init_x(options->x_mode, options->num_qbits_per_x, print);
			x_initialized = true;
		}

		if(!bin_initialized){
			init_expr_bin(options->bin_map, print);
			bin_initialized = true;
		}

		if(!pen_initialized){
			penalize_expr(options->penalty, options->pen_mode, print);
			pen_initialized = true;
		}

		if(!qubo_generated){
			generate_qubo(print);
			qubo_generated = true;
		}
}



/*xacc::quantum::PauliOperator Lattice::getHamiltonian(MapOptions* options){

	calcHamiltonian(options, options->verbose);

	std::map<int, std::pair<std::string, std::complex<double>>> operators;

	for(auto &term : expression_qubo->polynomial){

		std::pair<int, int> vars = term.first;
		if(vars.first == -1){

			if(vars.second == -1){
				//operators.emplace(expression_qubo->getQubit(var), std::pair<std::string, std::complex<double>>("Z", std::complex<double>(term.second.get_d(),0)));
			}else{
				operators.emplace(expression_qubo->getQubit(vars.second), std::pair<std::string, std::complex<double>>("Z", std::complex<double>(term.second.get_d(),0)));
			}

		}else{
			//operators.emplace(expression_qubo->getQubit(vars.second), std::pair<std::string, std::complex<double>>("Z", std::complex<double>(term.second.get_d(),0)));
		}


	}

	operators.emplace(0, std::pair<std::string, std::complex<double>>("Z", std::complex<double>(5,0)));
	operators.emplace(1, std::pair<std::string, std::complex<double>>("Z", std::complex<double>(6,0)));

	hamiltonian = xacc::quantum::PauliOperator(operators);

	return hamiltonian;

}*/

std::string Lattice::toHamiltonianString(){
	if(!qubo_generated){
		loge("Hamiltonian referenced but not yet generated!");
		return "";
	}
	return expression_qubo->expression_line_print();

}

std::string Lattice::toHamiltonianString(MapOptions* options){

	calcHamiltonian(options, options->verbose);
	return expression_qubo->expression_line_print();

}

