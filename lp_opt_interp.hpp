#ifndef LPOPTINTERP_HPP_
#define LPOPTINTERP_HPP_

#include "lp_opt.hpp"

namespace lp_opt
{
	class LP_OPT_INTERP : public LP_OPT
	{
	private:
		/* maximum level */
		int i_level_min;

		/* maximum level */
		int i_level_max;
	public:
		LP_OPT_INTERP() {}

		LP_OPT_INTERP(const int& _i_level_min, const int& _i_level_max, const int& _opt_type)
		{
			assert(_i_level_min >= 1);
			assert(_i_level_max >= 1);
			assert(_i_level_max >= _i_level_min);
			assert(_opt_type == GLP_MIN || _opt_type == GLP_MAX);

			i_level_min = _i_level_min;
			i_level_max = _i_level_max;
			opt_type = _opt_type;
		}

		virtual void init_opti_prob(const std::string& prob_name)
		{
			std::string name = python_code_caller("../python_code/main.py ", i_level_min, i_level_max);
			std::cout << system(name.c_str());
		}

		virtual void set_constr_matrix()
		{
			
		}

		virtual void set_constr_matrix(const std::vector<double>& W)
		{
			
		}

		virtual void solve_opti_problem() const
		{
			
		}

		virtual std::vector<double> get_results() const
		{
			std::vector<double> v;	

			return v;
		}

		virtual ~LP_OPT_INTERP()
		{

		}
	};
}

#endif /* LPOPTINTERP_HPP_ */
