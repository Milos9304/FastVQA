#include "run_paper_experiment.h"

SolutionDataset read_experiment_file(int num_ranks, int rank_min, int dim){

	std::ifstream expFile("../paper_experiment/out_higher_dims_small_dims.csv");
	std::ifstream matrixFile("../paper_experiment/out_higher_dims_small_dims_matrices.csv");
	if(expFile.fail() || matrixFile.fail()){
		loge("Failed opening experiments file");
	}
	SolutionDataset solutionDataset(num_ranks, rank_min, dim);

	std::string line;
	int inst_counter = 0;
	while(std::getline(expFile, line)){

		std::istringstream s(line);
		std::string field;
		while (getline(s, field,';')){ //iterate over different ranks

			Solution sol;

			//loge(s.str());
			std::string solution = field.substr(1, field.size()-3);

			int i = 0;
			while(solution[i++]!=',');

			double svLength = std::stod(solution.substr(0,i-1));
			std::string rest = solution.substr(i+2, solution.size()-1);

			std::vector<int> vectOfCoeffs;
			std::stringstream ss(rest);
		    for (std::string i; ss >> i;) {
		    	if(i[i.size()-1]==',')
		    		i=i.substr(0,i.size()-1);
		    	vectOfCoeffs.push_back(std::stoi(i));
		        while(ss.peek() == ' ')
		            ss.ignore();
		    }

		    sol.lattice_id = inst_counter;
		    sol.rank = vectOfCoeffs.size();
		    sol.svLength = svLength;
		    sol.coeffs = vectOfCoeffs;
		    solutionDataset.addDataset(sol);

		}
		inst_counter++;
	}

	int matrix_counter = 0;
	while(std::getline(matrixFile, line)){

		MatrixInt lattice;
		int rows = solutionDataset.rank_min+solutionDataset.num_ranks;
		int cols = solutionDataset.dim;
		lattice.resize(rows, cols);

		std::istringstream s(line);
		std::string field;

		int row_index = 0;
		while (getline(s, field,';')){ //iterate over martix rows

			//loge(s.str());
			//logw(field);

			//std::vector<int> vectOfCoeffs;
			std::stringstream ss(field.substr(1,field.size()-2));
			int col_index = 0;
			for (std::string i; ss >> i;) {
				if(i[i.size()-1]==',')
					i=i.substr(0,i.size()-1);
				//vectOfCoeffs.push_back(std::stoi(i));
				//loge(i);
				lattice(row_index, col_index++) = std::stoi(i);
				while(ss.peek() == ' ')
					ss.ignore();
			}
			row_index++;

			//lattice(matrix_counter, col_index++) = n;
			//return;

		}
		solutionDataset.addLattice(lattice);
		matrix_counter++;
	}

	expFile.close();
	matrixFile.close();

	return solutionDataset;
}

SolutionDataset run_paper_exp(int num_ranks, int rank_min, int dim){

	return read_experiment_file(num_ranks, rank_min, dim);

}
