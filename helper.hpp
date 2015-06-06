#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <cassert>

/* used for calling the python code as python script_name level_min level_max */
std::string python_code_caller(const std::string& script_name, const int& level_min, const int& level_max);
/* used to set the row and column variables for the optimization problem */
std::string set_aux_var_name(const std::string& var_name, const int& index);
/* used to generate random variables for the W matrix in the optimization problem */
std::vector<double> gen_rand(const int& size);
/* used to get data for the GCP when minimizing the interpolation error */
void get_python_data(const std::map<std::vector<double>, double >& comb_technique_dict);

#endif /* HELPER_HPP_ */