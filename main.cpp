#include "opt_combi_technique.hpp"

int main(int argv, char** argc)
{
	int level_min_x = 1;
	int level_min_y = 1;
	int level_max_x = 5;
	int level_max_y = 5;

	const std::string prob_name = "interp_optimization";

	lp_opt::LP_OPT_INTERP opt_interp(
		level_min_x, 
		level_min_y, 
		level_max_x, 
		level_max_y, 
		GLP_MIN);
	opt_interp.init_opti_prob(prob_name);

	return 0;
}