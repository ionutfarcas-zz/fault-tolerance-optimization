#include "helper.hpp"

std::string python_code_caller(const std::string& script_name, const int& level_min, const int& level_max)
{
	std::stringstream caller;

	caller << "python " << script_name << " " << level_min << " " << level_max;

	return caller.str();
}

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
	 	rand_var = 1e-2*(std::rand()%10);
	 	output.push_back(rand_var);
	}

	
	return output;
}

void get_python_data(const std::map<std::vector<double>, double >& comb_technique_dict)
{

}
