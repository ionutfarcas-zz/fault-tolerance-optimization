#include "helper.hpp"

std::string set_aux_var_name(const std::string& var_name, const int& index)
{
	std::stringstream aux_var;
	aux_var << var_name << index;

	return aux_var.str();
}

std::vector<double> gen_rand(const int& size)
{
	double rand_var = 0.0;
	std::vector<double> output;

	for(int i = 0 ; i < size ; ++i)
	{
	 	rand_var = std::rand()/(double)RAND_MAX;
	 	output.push_back(rand_var);
	}

	return output;
}
