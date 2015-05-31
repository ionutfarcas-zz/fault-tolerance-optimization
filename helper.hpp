#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <cassert>

/* used to set the row and column variables for the optimization problem */
std::string set_aux_var_name(const std::string& var_name, const int& index);
/* used to generate random variables for the W matrix in the optimization problem */
std::vector<double> gen_rand(const int& size);

#endif /* HELPER_HPP_ */