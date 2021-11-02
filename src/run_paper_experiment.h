/*
 * run_paper_exp.h
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_RUN_PAPER_EXP_H_
#define SRC_RUN_PAPER_EXP_H_

#include "vqaConfig.h"
#include <fstream>
#include <vector>
#include <utility>
#include "io/logger.h"

class Solution{
public:

	int rank;
	int lattice_id;

	double svLength;
	std::vector<int> coeffs;

};

class SolutionDataset{

public:
	int num_ranks;
	int rank_min;
	int dim;

	std::vector<MatrixInt> matrices;
	std::vector<Solution> dataset;

	SolutionDataset(int num_ranks, int rank_min, int dim){
		this->num_ranks = num_ranks;
		this->rank_min = rank_min;
		this->dim = dim;
	}

	void addDataset(Solution s){
		this->dataset.push_back(s);
	}

	void addLattice(MatrixInt l){
		this->matrices.push_back(l);
	}

	std::pair<std::vector<MatrixInt>, std::vector<Solution>> getMatricexAndDataset(){
		return std::pair<std::vector<MatrixInt>, std::vector<Solution>>(this->matrices, this->dataset);
	}

};

SolutionDataset read_experiment_file(int num_ranks, int rank_min, int dim);
SolutionDataset run_paper_exp(int num_ranks, int rank_min, int dim);

#endif /* SSRC_RUN_PAPER_EXP_H_ */
