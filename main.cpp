#include "opt_combi_technique.hpp"

int main(int argv, char** argc)
{
	std::cout << "Optimization problem started................" << std::endl;
	std::cout << "Please wait................................." << std::endl;
	std::cout << "Computing..................................." << std::endl;
	std::cout << "Results....................................." << std::endl;
	std::cout << std::endl;

	/* specify the name of the optimization problem */
	const std::string prob_name = "interp_based_optimization";

	/* specify levels for which the python code will be called */
	/* the order should be from lower to higher */
	vec2d levels = {{3, 4, 1}, {7, 8, 5}};
	assert(levels.front().size() == levels.back().size());
	
	int dim = levels[0].size();

	/* specify the faults as x and y coordinates; program will check whether the specified constraints are in the problem dictionary */
	vec2d faults = {{1, 2, 1}, {4, 5, 2}};

	/* output (i.e. c vector) of the optimization problem */
	std::vector<double> new_c;

	/* start the optimization problem */
	lp_opt::LP_OPT_INTERP opt_interp(
		levels, 
		dim,
		GLP_MAX,
		faults);
	opt_interp.init_opti_prob(prob_name);
	opt_interp.set_constr_matrix();
	opt_interp.solve_opti_problem();
	new_c = opt_interp.get_results();
	/* end the optimization problem */

	std::cout << std::endl;
	std::cout << "Optimization problem terminated succesfully!" << std::endl;
	return 0;
}