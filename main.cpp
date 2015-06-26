#include "opt_combi_technique.hpp"

int main(int argv, char** argc)
{
	std::cout << "Optimization problem started................" << std::endl;
	std::cout << std::endl;
	std::cout << "Please wait................................." << std::endl;
	std::cout << std::endl;
	std::cout << "Computing..................................." << std::endl;
	std::cout << std::endl;

	/* specify the name of the optimization problem */
	const std::string prob_name = "interp_based_optimization";

	/* specify level min and level max for which the python code will be called */
	/* for the moment, only implemented for the same x and y levels! */
	int level_min_x = 1;
	int level_min_y = 1;
	int level_max_x = 7;
	int level_max_y = 7;

	/* specify the faults as x and y coordinates; program will check whether the specified constraints are in the problem dictionary */
	vec2d faults = {{1, 1}, {3, 4}};

	/* output (i.e. c vector) of the optimization problem */
	std::vector<double> new_c;

	/* start the optimization problem */
	lp_opt::LP_OPT_INTERP opt_interp(
		level_min_x, 
		level_min_y, 
		level_max_x, 
		level_max_y, 
		GLP_MAX,
		faults);
	opt_interp.init_opti_prob(prob_name);
	opt_interp.set_constr_matrix();
	opt_interp.solve_opti_problem();
	new_c = opt_interp.get_results();
	/* end the optimization problem */

	std::cout << "Optimization problem terminated succesfully!" << std::endl;
	return 0;
}