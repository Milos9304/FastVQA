/*
 * symbolic_manipulation.cpp
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#include "logger.h"
#include "symbolic_manipulation.h"

void Expression::substitute(int id, std::map<int, int> subs_expr){

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

void Expression::print(){

	std::cout<<"\n  Optimization problem ["+name+"]:\n";

	std::cout << "\n    Expression:\n";
	bool firstTerm = true;
	std::cout << "	";
	for(auto &term : polynomial){

		if(term.second == 0)
			continue;

		int id1 = term.first.first;
		int id2 = term.first.second;

		if(firstTerm){
			term.second < 0 ? std::cout << "-" << term.second : std::cout << "" << term.second;
			firstTerm = false;
		}
		else{
			term.second < 0 ? std::cout << " - " << term.second * -1 : std::cout << " + " << term.second;
		}

		if(id1 == id2)
			std::cout << " " << idMap[id1]->name << "^2" ;
		else
			std::cout << " " << idMap[id1]->name << " " << idMap[id2]->name;

	}

	std::cout << "\n\n    Variables:\n";
	for(auto &var: variables){

		if(var->id >= 0) //do not print identity variable
			std::cout << "	" << var->name << " " << "[" << var->lb << ", " << var->ub << "]\n";
	}



	std::cout << "\n";

}
