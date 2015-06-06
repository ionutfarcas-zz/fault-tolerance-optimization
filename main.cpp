#include "lp_opt_err_split.hpp"
#include "lp_opt_interp.hpp"

int main(int argv, char** argc)
{
	const std::string prob_name = "error_splitting_optimization";

	std::vector<double> results;
	int m = 3;
	int n = 4;

	std::vector<double> rand_vars;
	rand_vars = gen_rand(m*n);

	lp_opt::LP_OPT_ERR_SPLIT opt(n, m, GLP_MIN);
	lp_opt::LP_OPT_INTERP opt_interp(2, 7, GLP_MIN);

	opt.init_opti_prob(prob_name);
	opt.set_constr_matrix(rand_vars);
	opt.solve_opti_problem();
	results = opt.get_results();

	opt_interp.init_opti_prob(prob_name);

	for(int i = 0 ; i < m ; ++i)
	{
		std::cout << "d" << i+1 << " = " << results[i] << " ; ";
	}
	std::cout << std::endl;

	for(int i = m ; i < n + m ; ++i)
	{
		std::cout << "c" << i+1 << " = " << results[i] << " ; ";
	}
	std::cout << std::endl;

	return 0;
}