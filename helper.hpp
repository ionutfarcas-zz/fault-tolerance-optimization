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
#include <algorithm>
#include <cassert>

typedef std::multimap<std::vector<int>, double> combi_grid_dict;
typedef std::vector<std::vector<int>> vec2d;

/* used to convert a string to a number in any format */
template<typename T>
T str_to_number(const std::string& no);
/* used for calling the python code as python script_name level_min level_max */
std::string python_code_caller(const std::string& script_name, const vec2d& levels);
/* used to get data for the GCP when minimizing the interpolation error */
combi_grid_dict get_python_data(const std::string& script_run, const int& dim);
/* used to create the M matrix for the interpolation based problem */
double** M_matrix(const combi_grid_dict& aux_downset);
/* used to set the N matrix = M - the identity */
double** N_matrix(const combi_grid_dict& aux_downset);
/* product of two matrices */
double** mat_prod(double** A, double** B, const int& dim);
/* used to compute powers of the matrix N; helpful when computing its inverse */
double** N_pow_k(double** N, const int& size_downset, const int& k);
/* used to compute the sum of the first n-1 powers of N, where n is the dimension of N */
double** sum_pow_N(double** N, const int& size_downset);
/* used to create the inverse of M matrix for the interpolation based problem */
double** M_inv(const combi_grid_dict& aux_downset);
/* used to create the inverse of M matrix for the interpolation based problem */
combi_grid_dict set_entire_downset_dict(
    const std::vector<int>& level_max, 
    const int& size_downset, 
    const std::string& script_run,
    const int& dim);
/* used to get a vector of the entire downset indices */
vec2d get_downset_indices(const combi_grid_dict& entire_downset, const int& dim);
/* used to filter the input vec2d faults such that only faults from the partial downset
(input from python) are considered */
vec2d filter_faults(
    const vec2d& faults_input, 
    const int& l_max, 
    const std::string& script_run,
    const int& dim);
/* used to create an entire downset dictionary used for setting up the M matrix */
combi_grid_dict create_aux_entire_dict(const combi_grid_dict& entire_downset, const int& dim);
/* used to print the new dictionary after the optimization is performed */
combi_grid_dict create_out_dict(const combi_grid_dict& given_downset, const std::vector<double>& new_c);
/* used to set the row and column variables for the optimization problem */
std::string set_aux_var_name(const std::string& var_name, const int& index);
/* used to generate random variables for the W matrix in the optimization problem */
std::vector<double> gen_rand(const int& size); 
/* used to compute the size of the downset */
int get_size_downset(const std::vector<int>& level_max);
/* used to compute the L1 norm of a vector */
int l1_norm(const std::vector<int>& u);
/* test whether b >= a */
bool test_greater(const std::vector<int>& b, const std::vector<int>& a);
/* used to create a multi-index based on maximum level */
vec2d mindex(const int& level_max);

#endif /* HELPER_HPP_ */