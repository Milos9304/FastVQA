#include "accelerator.h"

AlphaFunction Accelerator::alpha_constant_f = [](double init_val, double final_val, int iter_i, int max_iters){
	return init_val;
};

AlphaFunction Accelerator::alpha_linear_f = [](double init_val, double final_val, int iter_i, int max_iters){

	int steps=35;

	if(iter_i < max_iters*steps){
		int j = iter_i/steps;
		return init_val+((double)j/max_iters)*(final_val-init_val);
	}
	/*if(iter_i < max_iters)
		return init_val+((double)iter_i/max_iters)*(final_val-init_val);*/;

	return final_val;
};
