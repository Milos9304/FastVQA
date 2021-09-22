/*
 * run_paper_exp.h
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_RUN_PAPER_EXP_H_
#define SRC_RUN_PAPER_EXP_H_

#include "vqaConfig.h"

typedef std::pair<double, std::vector<int>> Solution;

class SolutionDataset{

public:
	int num_ranks;
	int rank_min;
	int dim;

	std::vector<Lattice> matrices;
	std::vector<Solution> dataset;

	SolutionDataset(int num_ranks, int rank_min, int dim){
		this->num_ranks = num_ranks;
		this->rank_min = rank_min;
		this->dim = dim;
	}

	void addDataset(Solution s){
		this->dataset.push_back(s);
	}

	void addLattice(Lattice l){
		this->matrices.push_back(l);
	}

};

void read_experiment_file(int num_ranks, int rank_min, int dim);
void run_paper_exp(int num_ranks, int rank_min, int dim);

#endif /* SSRC_RUN_PAPER_EXP_H_ */
