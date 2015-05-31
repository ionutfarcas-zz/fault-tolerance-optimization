#ifndef LPOPT_HPP_
#define LPOPT_HPP_

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cassert>

#include "helper.hpp"
#include "glpk.h"

namespace lp_opt
{
	class LP_OPT
	{
	private:
		/* problem name */
		const std::string prob_name = "error_splitting_optimization";

		/* length of the c vector */
		int i_n;
		/* length of the d vector */
		int i_m;
		/* optimization type: GLP_MIN or GLP_MAX */
		int opt_type;

		/* glp problem; used in every glpk function */
		glp_prob* i_lp_prob;

		/* constraint matrix */
		double* constr_mat;
		/* row index vector for the constraint matrix */
		int* row_index;
		/* colums index vector for the constraint matrix */
		int* col_index;

		/* helper variables, used to compute the dimension of the constraint matrix */
		int rows;
		int cols;
		int total_size;

	public:
		/* default constructor */
		LP_OPT();
		/* additional constructor */
		/* used to initialize the problem parameters n and m, as well as the optimization type (MIN or MAX) */
		/* moreover, size of the constraint matrix is computed and the its memory allocation is performed */
		LP_OPT(const int& _i_n, const int& _i_m, const int& _opt_type);

		/* used for LP optimization problem initialization */
		/* number of auxiliary and structural variables is set */
		/* as well as objective function and the constraints */
		void init_opti_prob();

		/* used to set up the constraint matrix and its row and column indices */
		/*  the constraint matrix has the form */
		/* ----------- */
		/* |  W | -I | */
 		/* ----------- */
 		/* | -W | -I | */
		/* ----------- */
		/* | 0m | 1n | */
		/* ----------- */ 
		void set_constr_matrix();

		/* used to solve the linear programming problem, using the simplex algorithm */
		void solve_opti_problem() const;

		/* used to get the output of the LP problem, i.e. c and d vectors */
		std::vector<double> get_results() const;

		/* destructor; used to deallocate all the allocated memory */
		~LP_OPT();

	};
}
#endif /* LP_HPP_ */
