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
#include <gmpxx.h>

class Var{

	public:

		//lb = -1, ub = 1 => qubo Z var
		Var(int id, std::string name, int lowerBound, int upperBound, int qubit=-1){
			this -> id = id;
			this -> lb = lowerBound;
			this -> ub = upperBound;
			this -> name = name;
			this -> qubit = qubit;
		}

		std::string name;

		int lb, ub;
		int id;
		int qubit;

		bool isBinary(){
			return lb == 0 && ub == 1;
		}

};

class Expression{

	private:

		int id_counter = 0;
		std::vector<Var*> variables;
		std::map<int, Var*> idMap;
		std::map<std::string, int> varMap;


	public:

		std::string name;

		//id, id, coeff
		std::map<std::pair<int, int>, mpq_class> polynomial;

		Expression(std::string name){

			Var* idVar = new Var(-1, "id", 1, 1);
			variables.push_back(idVar); //identity variable
			idMap.emplace(-1, idVar);
			varMap.emplace("id", -1);

			this -> name = name;
		}

		int getQubit(int id){
			return idMap[id]->qubit;
		}

		int addIntegerVar(std::string name, int lowerBound, int upperBound, int qubit=-1){
			Var* var = new Var(id_counter, name, lowerBound, upperBound, qubit);
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
			return std::pair<int, std::string>(this->addIntegerVar(z_name, 1, -1, qubit), z_name);
		}

		void addNewTerm(int id_a, int id_b, mpz_t coeff){
					addNewTerm(id_a, id_b, mpz_class(coeff));
		}


		void addNewTerm(int id_a, int id_b, mpq_t coeff){
			addNewTerm(id_a, id_b, mpq_class(coeff));
		}

		void addNewTerm(int id_a, int id_b, mpq_class coeff){

			if(id_a == id_b && idMap[id_a]->isBinary())
				polynomial.emplace(std::pair<int, int>(-1, id_a), coeff);
			else if(id_a <= id_b)
				polynomial.emplace(std::pair<int, int>(id_a, id_b), coeff);
			else
				polynomial.emplace(std::pair<int, int>(id_b, id_a), coeff);
		}

		void addTermCoeff(int id_a, int id_b, mpq_class coeff){

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
						//double coeff2 = search->second;
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

		void addConstant(mpq_class constant){
			polynomial[std::pair<int,int>(-1,-1)] += constant;
		}

		void substituteVarToDouble(int id, mpq_class val){
			std::map<int, mpq_class> subs_expr;
			subs_expr.emplace(-1, val);
			substitute(id, subs_expr);
		}

		void substitute(int id, std::map<int, mpq_class> subs_expr);

		std::vector<Var*> getVariables(){
			return variables;
		}

		int getId(std::string name){
			return varMap[name];
		}

		std::string getName(int id){
			return idMap[id]->name;
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
