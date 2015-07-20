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
	vec2d levels = {{1, 1}, {3, 3}};
	assert(levels.front().size() == levels.back().size());
	
	int dim = levels[0].size();

	/* specify the faults as x and y coordinates; program will check whether the specified constraints are in the problem dictionary */
	vec2d faults = {{2, 2}};
	assert(faults[0].size() == static_cast<unsigned int>(dim));

	/* output (i.e. c vector) of the optimization problem */
	std::vector<double> new_c;

	/* start the optimization problem */
	auto t1 = std::chrono::high_resolution_clock::now();
	lp_opt::LP_OPT_INTERP opt_interp(
		levels, 
		dim,
		GLP_MAX,
		faults);
	auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Constructor call "
              << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()
              << " seconds\n";
    auto t3 = std::chrono::high_resolution_clock::now();
	opt_interp.init_opti_prob(prob_name);
	auto t4 = std::chrono::high_resolution_clock::now();
    std::cout << "init_opti_prob call "
              << std::chrono::duration_cast<std::chrono::seconds>(t4-t3).count()
              << " seconds\n";
    auto t5 = std::chrono::high_resolution_clock::now();
	opt_interp.set_constr_matrix();
	auto t6 = std::chrono::high_resolution_clock::now();
    std::cout << "set_constr_matrix "
              << std::chrono::duration_cast<std::chrono::seconds>(t6-t5).count()
              << " seconds\n";
    auto t7 = std::chrono::high_resolution_clock::now();
	opt_interp.solve_opti_problem();
	auto t8 = std::chrono::high_resolution_clock::now();
    std::cout << "solve_opti_problem "
              << std::chrono::duration_cast<std::chrono::seconds>(t8-t7).count()
              << " seconds\n";
	new_c = opt_interp.get_results();
	/* end the optimization problem */

	std::cout << std::endl;
	std::cout << "Optimization problem terminated succesfully!" << std::endl;

	std::cout << "Results " << std::endl;
	for(unsigned int i = 0 ; i < new_c.size() ; ++i)
	{
		std::cout << "c[" << i << "] = " << new_c[i] << std::endl;
	}
	return 0;
}