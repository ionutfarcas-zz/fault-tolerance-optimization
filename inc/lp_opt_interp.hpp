#ifndef LPOPTINTERP_HPP_
#define LPOPTINTERP_HPP_

#include "lp_opt.hpp"

namespace lp_opt
{
	class LP_OPT_INTERP : public LP_OPT
	{
	private:
		/* levels of grid indices */
		vec2d i_levels;
		/* top level of grid indices*/
		std::vector<int> level_max;
		/* dimension of the problem */
		int i_dim;
		/* total size of the optimization problem */
		int total_size;

		/* new dimensionality after the input levels are checked */
		int new_dim;
		/* new levels based on new_dim */
		vec2d new_levels;
		/* if it is the case, the dimension(s) that is/are ignored */
		std::vector<int> ignored_dimensions;
		/* if it is the case, new faults, based on ignored dimensions */
		vec2d new_faults;

		/* no. of constraints */
		int no_faults;
		/* level max sum */
		int l_max;

		//  total size of the optimization problem 
		/* down set size */
		int size_downset;

		/* inverse of M s.t. w = Mc*/
		matrix inv_M;

		/* given dictionary */
		combi_grid_dict given_downset;
		/* modified given downset based on ignored dimensions */
		combi_grid_dict new_given_downset;
		/* entire donwset with corresponding indices */
		combi_grid_dict entire_downset;
		/* auxiliary dictionary, used to create M and inv(M) */
		combi_grid_dict aux_entire_dict;
		/* downset indices as a 2d vector*/
		vec2d downset_indices;
		/* faults input by user */
		vec2d input_faults;
		/* faults that are in the given downset from the set of input faults */
		vec2d valid_input_faults;

	public:
		LP_OPT_INTERP();

		LP_OPT_INTERP(
			const vec2d& _levels,  
			const int& _dim, 
			const int& _opt_type,
			const combi_grid_dict& _given_downset,
			const vec2d& _input_faults);

		LP_OPT_INTERP(const LP_OPT_INTERP& obj);

		LP_OPT_INTERP& operator= (const LP_OPT_INTERP& rhs);

		virtual void init_opti_prob(const std::string& prob_name);
		
		virtual void set_constr_matrix();
	
		virtual void set_constr_matrix(const std::vector<double>& W);
	
		virtual void solve_opti_problem() const;

		virtual std::vector<double> get_results() const;

		virtual ~LP_OPT_INTERP();
	};
}

#endif /* LPOPTINTERP_HPP_ */