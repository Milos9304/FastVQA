/*
 * symbolic_manipulation.h
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_SYMBOLIC_MANIPULATION_H_
#define SRC_SYMBOLIC_MANIPULATION_H_

#include <map>
#include <vector>
#include <string>
#include <algorithm>

class Var{

	public:

		//lb = -1, ub = 1 => qubo Z var
		Var(int id, std::string name, int lowerBound, int upperBound){
			this -> id = id;
			this -> lb = lowerBound;
			this -> ub = upperBound;
			this -> name = name;
		}

		std::string name;

		int lb, ub;
		int id;

};

class Expression{

	private:

		int id_counter = 0;
		std::vector<Var*> variables;
		std::map<int, Var*> idMap;
		std::map<std::string, int> varMap;

		//id, id, coeff
		std::map<std::pair<int, int>, double> polynomial;

	public:

		std::string name;

		Expression(std::string name){

			Var* idVar = new Var(-1, "id", 1, 1);
			variables.push_back(idVar); //identity variable
			idMap.emplace(-1, idVar);
			varMap.emplace("id", -1);

			this -> name = name;
		}

		int addIntegerVar(std::string name, int lowerBound, int upperBound){
			Var* var = new Var(id_counter, name, lowerBound, upperBound);
			variables.push_back(var);
			idMap.emplace(id_counter, var);
			varMap.emplace(name, id_counter);
			return id_counter++;
		}

		int addBinaryVar(std::string name){
			return this->addIntegerVar(name, 0, 1);
		}

		std::pair<int, std::string> addZ(int qubit){ //when creating qubo formulation
			std::string z_name = "Z"+std::to_string(qubit);
			return std::pair<int, std::string>(this->addIntegerVar(z_name, 1, -1), z_name);
		}

		void addNewTerm(int id_a, int id_b, double coeff){

			if(id_a == id_b)
				polynomial.emplace(std::pair<int, double>(id_a, -1), coeff);
			else if(id_a < id_b)
				polynomial.emplace(std::pair<int, double>(id_a, id_b), coeff);
			else
				polynomial.emplace(std::pair<int, double>(id_b, id_a), coeff);
		}

		void addTermCoeff(int id_a, int id_b, double coeff){

			std::pair<int, int> to_search;
			if(id_a <= id_b)
				to_search = std::pair<int, double>(id_a, id_b);
			else
				to_search = std::pair<int, double>(id_b, id_a);

			auto search = polynomial.find(to_search);
			if (search != polynomial.end()) {
				int coeff2 = search->second;
				polynomial.erase(to_search);
				addNewTerm(id_a, id_b, coeff + coeff2);
			} else {
				addNewTerm(id_a, id_b, coeff);
			}

		}

		void addConstant(double constant){
			polynomial[std::pair<int,int>(-1,-1)] += constant;
		}

		void substituteVarToDouble(int id, double val){
			std::map<int, double> subs_expr;
			subs_expr.emplace(-1, val);
			substitute(id, subs_expr);
		}

		void substitute(int id, std::map<int, double> subs_expr);

		std::vector<Var*> getVariables(){
			return variables;
		}

		int getId(std::string name){
			return varMap[name];
		}

		void delId(int id){
			Var* var = idMap[id];
			varMap.erase(var->name);
			std::vector<Var*>::iterator it = std::find(variables.begin(), variables.end(), var);
			variables.erase(it);
			idMap.erase(id);
		}

		std::string expression_line_print();
		void print();

};
#endif /* SRC_SYMBOLIC_MANIPULATION_H_ */
