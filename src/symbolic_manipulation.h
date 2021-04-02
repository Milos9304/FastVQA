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

class Expression{

	private:

		class Var{

			public:
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

		int id_counter = 0;
		std::vector<Var*> variables;
		std::map<int, Var*> idMap;
		std::map<std::string, int> varMap;

		//id, id, coeff
		std::map<std::pair<int, int>, int> polynomial;

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

		void addNewTerm(int id_a, int id_b, int coeff){
			if(id_a <= id_b)
				polynomial.emplace(std::pair<int, int>(id_a, id_b), coeff);
			else
				polynomial.emplace(std::pair<int, int>(id_b, id_a), coeff);
		}

		void addTermCoeff(int id_a, int id_b, int coeff){

			std::pair<int, int> to_search;
			if(id_a <= id_b)
				to_search = std::pair<int, int>(id_a, id_b);
			else
				to_search = std::pair<int, int>(id_b, id_a);

			auto search = polynomial.find(to_search);
			if (search != polynomial.end()) {
				int coeff2 = search->second;
				polynomial.erase(to_search);
				addNewTerm(id_a, id_b, coeff + coeff2);
			} else {
				addNewTerm(id_a, id_b, coeff);
			}

		}

		void substitute(int id, std::map<int, int> subs_expr);

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

		void print();

};
#endif /* SRC_SYMBOLIC_MANIPULATION_H_ */
