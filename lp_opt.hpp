#ifndef LPOPT_HPP_
#define LPOPT_HPP_

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <cassert>

#include "glpk.h"
#include "helper.hpp"

namespace lp_opt
{
	class LP_OPT
	{
	protected:
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

	public:
		/* used for LP optimization problem initialization */
		/* number of auxiliary and structural variables is set */
		/* as well as objective function and the constraints */
		virtual void init_opti_prob(const std::string& prob_name) = 0;

		/* used to set up the constraint matrix and its row and column indices when it depends on faults*/ 
		virtual void set_constr_matrix(const vec2d& faults) = 0;

		/* used to set up the constraint matrix and its row and column indices when it depends on matrix W*/ 
		virtual void set_constr_matrix(const std::vector<double>& W) = 0;

		/* used to solve the linear programming problem, using the simplex algorithm */
		virtual void solve_opti_problem() const = 0;

		/* used to get the output of the LP problem, i.e. c and d vectors */
		virtual std::vector<double> get_results() const = 0;

		/* destructor; used to deallocate all the allocated memory */
		virtual ~LP_OPT() {}
	};
}

#endif /* LPOPT_HPP_ */
