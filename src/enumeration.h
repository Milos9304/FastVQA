/*
 * enumeration.h
 *
 *  Created on: May 3, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_ENUMERATION_H_
#define SRC_ENUMERATION_H_

#include "indicators/progress_bar.hpp"
#include "lattice/lattice.h"

class EnumerationC{

public:

	int n;
	int n_qubits; //<=n, if e.g. penalty not used, less qubits needed
	int rank;
	bool* fixed;
	int qubits_fixed;

	bool ignoreZero = false;

	std::vector<double> coeffs;
	std::vector<int> pauliOpts;

	Lattice* lattice;
	ProgressBar* bar;

	double shortestVectLen=9999999999;

	// Function to print the output
	void updateShortestVect(bool* arr){

		if(rank == 0)
			bar->tick();

		double result = 0;
			//cout << arr[i] << " ";

			//VectorInt x = lattice->quboToXvector(arr, n);
			for(int i = 0; i < coeffs.size(); ++i){

				int val = 1;
				for(int j = i * n; j < i*n+n_qubits; ++j){

					if(pauliOpts[j] == 3){

						if(j - (i * n) < qubits_fixed){
							if(fixed[j - (i * n)]){
								val *= -1;
							}
						}

						else{
							if(arr[j - qubits_fixed - (i * n)]){
								val *= -1;
							}
						}
					}
				}
				result += coeffs[i] * val;
			}

			if(ignoreZero && result == 0)
				return;

			shortestVectLen = result < shortestVectLen ? result : shortestVectLen;

		//cout << endl;
	}

	// Function to generate all binary strings
	void generateAllBinaryStrings(bool* arr, int i){

		if (i == n_qubits-qubits_fixed) {
			updateShortestVect(arr);
			return;
		}

		arr[i] = 0;
		generateAllBinaryStrings(arr, i + 1);

		arr[i] = 1;
		generateAllBinaryStrings(arr, i + 1);
	}
};


#endif /* SRC_ENUMERATION_H_ */
