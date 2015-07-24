#ifndef LPOPTERRSPLIT_HPP_
#define LPOPTERRSPLIT_HPP_

#include "lp_opt.hpp"

namespace lp_opt
{
	class LP_OPT_ERR_SPLIT : public LP_OPT
	{
	private:
		/* length of the c vector */
		int i_n;
		/* length of the d vector */
		int i_m;
		
		/* helper variables, used to compute the dimension of the constraint matrix */
		int rows;
		int cols;
		int total_size;

	public:
		LP_OPT_ERR_SPLIT();

		LP_OPT_ERR_SPLIT(const int& _i_n, const int& _i_m, const int& _opt_type);

		LP_OPT_ERR_SPLIT(const LP_OPT_ERR_SPLIT& obj);

		LP_OPT_ERR_SPLIT& operator= (const LP_OPT_ERR_SPLIT& rhs);

		virtual void init_opti_prob(const std::string& prob_name);

		virtual void set_constr_matrix();

		/* the constraint matrix has the form */
		/* ----------- */
		/* |  W | -I | */
 		/* ----------- */
 		/* | -W | -I | */
		/* ----------- */
		/* | 0m | 1n | */
		/* ----------- */
		virtual void set_constr_matrix(const std::vector<double>& W);

		virtual void solve_opti_problem() const;

		virtual std::vector<double> get_results() const;

		virtual ~LP_OPT_ERR_SPLIT();
	};
}

#endif /* LPOPTERRSPLIT_HPP_ */