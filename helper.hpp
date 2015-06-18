#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <cassert>

typedef std::multimap<std::vector<int>, double> combi_grid_dict;

/* used to convert a string to a number in any format */
template<typename T>
T str_to_number(const std::string& no);
/* used for calling the python code as python script_name level_min level_max */
std::string python_code_caller(
    const std::string& script_name, 
    const int& level_min_x,
    const int& level_min_y, 
    const int& level_max_x,
    const int& level_max_y);
/* used to get data for the GCP when minimizing the interpolation error */
combi_grid_dict get_python_data(const std::string& script_run);
/* used to create the M matrix for the interpolation based problem */
double** M_matrix(const combi_grid_dict& aux_downset);
/* used to set the N matrix = M - the identity */
double** set_N_matrix(const combi_grid_dict& aux_downset);
/* product of two matrices */
double** mat_prod(double** A, double** B, const int& dim);
/* used to compute powers of the matrix N; helpful when computing its inverse */
double** N_pow_k(double** N, const int& size_downset, const int& k);
/* used to compute the sum of the first n-1 powers of N, where n is the dimension of N */
double** sum_pow_N(double** N, const int& size_downset);
/* used to create the inverse of M matrix for the interpolation based problem */
double** M_inv(const combi_grid_dict& aux_downset);
/* used to create the inverse of M matrix for the interpolation based problem */
combi_grid_dict entire_downset_dict(
    const int& level_max_x,
    const int& level_max_y,  
    const int& size_downset, 
    const std::string& script_run);
/* used for setting up the M matrix */
combi_grid_dict aux_dict(const combi_grid_dict& entire_downset);
/* used to set the row and column variables for the optimization problem */
std::string set_aux_var_name(const std::string& var_name, const int& index);
/* used to generate random variables for the W matrix in the optimization problem */
std::vector<double> gen_rand(const int& size); 
/* used to compute the size of the down set for the case level_max_x = level_max_y */
int get_size_downset(const int& level_max_x, const int& level_max_y);
/* test whether j >= i */
bool test_greater(const std::vector<int>& j, const std::vector<int>& i);

#endif /* HELPER_HPP_ */