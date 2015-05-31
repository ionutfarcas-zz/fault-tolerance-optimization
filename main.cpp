#include "opt_combi_technique.hpp"

int main(int argv, char** argc)
{
	std::vector<double> results;
	int m = 3;
	int n = 4;

	lp_opt::LP_OPT opt(n, m, GLP_MIN);

	opt.init_opti_prob();
	opt.set_constr_matrix();
	opt.solve_opti_problem();
	results = opt.get_results();

	for(int i = 0 ; i < m ; ++i)
	{
		std::cout << "d" << i+1 << "= " << results[i] << std::endl;
	}
	for(int i = m ; i < n + m ; ++i)
	{
		std::cout << "c" << i+1 << "= " << results[i] << std::endl;
	}

	return 0;
}