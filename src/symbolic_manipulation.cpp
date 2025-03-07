/*
 * symbolic_manipulation.cpp
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#include "symbolic_manipulation.h"
#include <sstream>
#include <iomanip>
#include "logger.h"
#include <iostream>

namespace FastVQA{

void Expression::substitute(int id, std::map<int, mpq_class> subs_expr){

	std::vector<std::pair<int, int>> toDel;

	for(auto &term : polynomial){

		int second_id;

		if(term.first.first == id){
			second_id = term.first.second;
		}else if(term.first.second == id){
			second_id = term.first.first;
		}else{
			continue;
		}

		if(id == second_id){ //^2 terms
			for(auto &subs_var : subs_expr){
				for(auto &subs_var2 : subs_expr){
					this -> addTermCoeff(subs_var.first, subs_var2.first, subs_var.second * subs_var2.second * term.second);
				}
		}

		}else{
			for(auto &subs_var : subs_expr){
				this -> addTermCoeff(second_id, subs_var.first, subs_var.second * term.second);
			}
		}

		toDel.push_back(term.first);
	}

	for(auto &d : toDel){
		polynomial.erase(d);
	}

	delId(id);

}

std::string Expression::expression_line_print(){

	std::stringstream ss;
	ss << std::fixed << std::setprecision(2);

	bool firstTerm = true;
	ss << "	";

	for(auto &term : polynomial){

		if(term.second == 0)
			continue;

		int id1 = term.first.first;
		int id2 = term.first.second;

		if(firstTerm){
			ss << mpf_class(term.second);
			firstTerm = false;
		}
		else{
			term.second < 0 ? ss << " - " << mpf_class(term.second * -1) : ss << " + " << mpf_class(term.second);
		}

		if(id1 == -1 && id2 == -1){} //id
		else if(id1 == -1)
			ss << " " << idMap[id2]->name;
		else if(id2 == -1)
			ss << " " << idMap[id1]->name;
		else
			ss << " " << idMap[id1]->name << " " << idMap[id2]->name;

	}

	return ss.str();
}

mpq_class Expression::evaluate_bin_expr(std::map<Var*, int> *varBoolMap){

	mpq_class result = 0;
	for(auto &term : polynomial){

		if(term.second == 0)
			continue;

		int id1 = term.first.first;
		int id2 = term.first.second;

		if(id1 == -1 && id2 == -1){
			result += term.second;
		} //id
		else if(id1 == -1){
			if(!idMap[id2]->isBinary())
				throw_runtime_error("evaluate_bin_expr only for binary variables");
			result += term.second * (*varBoolMap)[idMap[id2]];
		}else if(id2 == -1)
			result += term.second * (*varBoolMap)[idMap[id1]];
		else
			result += term.second * (*varBoolMap)[idMap[id1]] * (*varBoolMap)[idMap[id2]];

	}

	return result;
}

void Expression::print(){

	std::cout<<"\n  Optimization problem ["+name+"]:\n";
	std::cout << "\n    Expression:\n";
	std::cout << expression_line_print();
	std::cout << "\n\n    Variables (total " << variables.size() - 1 << "):\n";
	for(auto &var: variables){

		if(var->lb == 1 && var->ub == -1) //Z var
			std::cout << "	" << var->name << "\n";

		else if(var->id >= 0) //do not print identity variable
			std::cout << "	" << var->name << " " << "[" << var->lb << ", " << var->ub << "]" << " " << var->extra_information << "\n";

	}
	std::cout << "\n";
}
}

