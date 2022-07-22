/*
 * symbolic_manipulation.cpp
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#include "symbolic_manipulation.h"
#include <sstream>
#include <iomanip>
#include <iostream>

namespace fastVQA{

Expression::Expression(std::string name){

	Var* idVar = new Var(-1, "id", 1, 1);
	variables.push_back(idVar); //identity variable
	idMap.emplace(-1, idVar);
	varMap.emplace("id", -1);

	this -> name = name;
}

int Expression::addIntegerVar(std::string name, int lowerBound, int upperBound, int qubit){
	Var* var = new Var(id_counter, name, lowerBound, upperBound, qubit);
	variables.push_back(var);
	idMap.emplace(id_counter, var);
	varMap.emplace(name, id_counter);
	return id_counter++;
}

std::pair<int, std::string> Expression::addZ(int qubit){
	std::string z_name = "Z"+std::to_string(qubit);
	return std::pair<int, std::string>(this->addIntegerVar(z_name, 1, -1, qubit), z_name);
}

void Expression::addNewTerm(int id_a, int id_b, mpq_class coeff){

	if(id_a == id_b && idMap[id_a]->isBinary())
		polynomial.emplace(std::pair<int, int>(-1, id_a), coeff);
	else if(id_a <= id_b)
		polynomial.emplace(std::pair<int, int>(id_a, id_b), coeff);
	else
		polynomial.emplace(std::pair<int, int>(id_b, id_a), coeff);
}

void Expression::substituteVarByNumeric(int id, mpq_class val){
	std::map<int, mpq_class> subs_expr;
	subs_expr.emplace(-1, val);
	substitute(id, subs_expr);
}

void Expression::delId(int id){
	Var* var = idMap[id];
	varMap.erase(var->name);
	std::vector<Var*>::iterator it = std::find(variables.begin(), variables.end(), var);
	variables.erase(it);
	idMap.erase(id);
}


void Expression::addTermCoeff(int id_a, int id_b, mpq_class coeff){

	std::pair<int, int> to_search;
	if(id_a <= id_b)
		to_search = std::pair<int, int>(id_a, id_b);
	else
		to_search = std::pair<int, int>(id_b, id_a);

	mpq_class coeff2;

	auto search = polynomial.find(to_search);
	if (search != polynomial.end()) {
		coeff2 = search->second;
		polynomial.erase(to_search);
		addNewTerm(id_a, id_b, coeff + coeff2);
	} else {

		//maybe it is a binary variable where x^2=x;
		if(id_a == id_b && idMap[id_b]->isBinary()){

			to_search = std::pair<int, int>(-1, id_b);
			auto search = polynomial.find(to_search);
			if (search != polynomial.end()) {
				coeff2 = search->second;
				polynomial.erase(to_search);
				addNewTerm(id_a, id_b, coeff + coeff2);
				return;
			}

		}else if(id_a == -1 && idMap[id_b]->isBinary()){

			to_search = std::pair<int, int>(id_b, id_b);
			auto search = polynomial.find(to_search);
			if (search != polynomial.end()) {
				coeff2 = search->second;
				polynomial.erase(to_search);
				addNewTerm(id_a, id_b, coeff + coeff2);
				addNewTerm(id_a, id_b, coeff);
				return;
			}
		}
		addNewTerm(id_a, id_b, coeff);
	}
}

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

void Expression::print(){

	std::cout<<"\n  Optimization problem ["+name+"]:\n";
	std::cout << "\n    Expression:\n";
	std::cout << expression_line_print();
	std::cout << "\n\n    Variables (total " << variables.size() - 1 << "):\n";
	for(auto &var: variables){

		if(var->lb == 1 && var->ub == -1) //Z var
			std::cout << "	" << var->name << "\n";

		else if(var->id >= 0) //do not print identity variable
			std::cout << "	" << var->name << " " << "[" << var->lb << ", " << var->ub << "]\n";

	}
	std::cout << "\n";
}
}
