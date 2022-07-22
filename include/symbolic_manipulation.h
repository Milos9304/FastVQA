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

namespace fastVQA{

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
		Expression(std::string name);

		//id, id, coeff
		std::map<std::pair<int, int>, mpq_class> polynomial;

		int getIdMapSize(){return idMap.size();}

		std::vector<Var*> getVariables(){return variables;}

		int getId(std::string name){return varMap[name];}

		std::string getName(int id){return idMap[id]->name;}

		void delId(int id);

		int getQubit(int id){return idMap[id]->qubit;}

		int addIntegerVar(std::string name, int lowerBound, int upperBound, int qubit=-1);
		int addBinaryVar(std::string name){return this->addIntegerVar(name, 0, 1);}

		//For qubo formulation
		std::pair<int, std::string> addZ(int qubit);

		//Does not check for pre-existing terms.
		void addNewTerm(int id_a, int id_b, mpz_t coeff){addNewTerm(id_a, id_b, mpz_class(coeff));}
		void addNewTerm(int id_a, int id_b, mpq_t coeff){addNewTerm(id_a, id_b, mpq_class(coeff));}
		void addNewTerm(int id_a, int id_b, mpq_class coeff);

		//Checks if the term already exists. If yes, the coeff gets modified
		void addTermCoeff(int id_a, int id_b, mpq_class coeff);

		void addConstant(mpq_class constant){polynomial[std::pair<int,int>(-1,-1)] += constant;}

		void substituteVarByNumeric(int id, mpq_class val);

		void substitute(int id, std::map<int, mpq_class> subs_expr);

		std::string expression_line_print();
		void print();
};
}
#endif /* SRC_SYMBOLIC_MANIPULATION_H_ */
