#ifndef LPOPTINTERP_HPP_
#define LPOPTINTERP_HPP_

#include "lp_opt.hpp"

namespace lp_opt
{
	class LP_OPT_INTERP : public LP_OPT
	{
	private:
		/* name of the python script */
		const std::string script_name = "../python_code/main.py";
		/* python caller */
		std::string get_dict;

		/* maximum level x direction*/
		int i_level_min_x;
		/* maximum level y direction*/
		int i_level_min_y;
		/* maximum level x direction*/
		int i_level_max_x;
		/* maximum level y direction*/
		int i_level_max_y;

		/* no. of constraints */
		int i_no_faults;

		/* total size of the optimization problem */
		int total_size;
		/* down set size */
		int size_downset;

		/* inverse of M s.t. w = Mc*/
		double** inv_M;

		/* entire donwset with corresponding indices */
		combi_grid_dict entire_downset;
		/* auxiliary dictionary, used to create M and inv(M) */
		combi_grid_dict aux_entire_dict;
		/* downset indices as a 2d vector*/
		vec2d downset_indices;

	public:
		LP_OPT_INTERP() {}

		LP_OPT_INTERP(
			const int& _i_level_min_x,
			const int& _i_level_min_y,  
			const int& _i_level_max_x,
			const int& _i_level_max_y,
			const int& _i_no_faults, 
			const int& _opt_type)
		{
			assert(_i_level_min_x >= 1);
			assert(_i_level_min_y >= 1);
			assert(_i_level_max_x >= 1);
			assert(_i_level_max_y >= 1);
			assert(_i_no_faults >= 0);

			assert(_i_level_max_x >= _i_level_min_x);
			assert(_i_level_max_y >= _i_level_min_y);

			assert(_opt_type == GLP_MIN || _opt_type == GLP_MAX);

			i_level_min_x = _i_level_min_x;
			i_level_min_y = _i_level_min_y;
			i_level_max_x = _i_level_max_x;
			i_level_max_y = _i_level_max_y;

			i_no_faults = _i_no_faults;
			opt_type = _opt_type;

			size_downset = get_size_downset(_i_level_max_x, _i_level_max_y);
			total_size = _i_no_faults*size_downset;

			get_dict = python_code_caller(
				script_name, 
				i_level_min_x, 
				i_level_min_y,  
				i_level_max_x, 
				i_level_max_y);

			entire_downset = set_entire_downset_dict(
				i_level_max_x,
				i_level_max_y,  
				size_downset, 
				get_dict);

			aux_entire_dict = create_aux_entire_dict(entire_downset);
			inv_M = M_inv(aux_entire_dict);

			downset_indices = get_donwset_indices(entire_downset);

			constr_mat = (double*)malloc((1 + total_size)*sizeof(double));
			row_index = (int*)malloc((1 + total_size)*sizeof(int));
			col_index = (int*)malloc((1 + total_size)*sizeof(int));

			assert(constr_mat!= NULL);
			assert(row_index!= NULL);
			assert(col_index!= NULL);
		}

		virtual void init_opti_prob(const std::string& prob_name)
		{
			std::vector<int> index;
			std::string aux_var;
			double neg_norm = 0.0;
			double coeff = 0.0;

			i_lp_prob = glp_create_prob();
			assert(i_lp_prob != NULL);

			glp_set_prob_name(i_lp_prob, prob_name.c_str());
			glp_set_obj_dir(i_lp_prob, opt_type);

			glp_add_rows(i_lp_prob, i_no_faults);
			glp_add_cols(i_lp_prob, size_downset);

			for(int i = 0 ; i < i_no_faults; ++i)
			{
				aux_var = set_aux_var_name("eq_constr_", i + 1);
				glp_set_row_name(i_lp_prob, i + 1, aux_var.c_str());
				glp_set_row_bnds(i_lp_prob, i + 1, GLP_FX, 0.0, 0.0);
			}

			for(int i = 0 ; i < size_downset ; ++i)
   			{
   				index = {downset_indices[i][0], downset_indices[i][1]};
   				neg_norm = -l1_norm(index);
   				coeff = pow(4.0, neg_norm);

				aux_var = set_aux_var_name("w", i + 1);
				glp_set_col_name(i_lp_prob, i + 1, aux_var.c_str());
				glp_set_col_bnds(i_lp_prob, i + 1, GLP_FR, 0.0, 0.0);
				glp_set_obj_coef(i_lp_prob, i + 1, coeff);
			}
		}

		virtual void set_constr_matrix(const vec2d& faults)
		{
			int inv_M_row_index = 0;
			std::vector<int> fault;

			std::cout << "Constraint matrix" << std::endl;
			for(int i = 0 ; i < i_no_faults ; ++i)
			{
				fault = {faults[i][0], faults[i][1]};
				auto it = aux_entire_dict.find(fault);

				if(it != aux_entire_dict.end())
				{
					inv_M_row_index = static_cast<int>(it->second);

					for(int j = 0 ; j < size_downset ; ++j)
					{
						std::cout << inv_M[inv_M_row_index][j] << " ";	
					}
					std::cout << std::endl;
				}
			}
		}

		virtual void set_constr_matrix(const std::vector<double>& W)
		{	
			/* nothing to do here */
		}

		virtual void solve_opti_problem() const
		{
			/* TO - DO :call one lp optimization routine after setting up the constraints matrix and coeficient vectors */
		}

		virtual std::vector<double> get_results() const
		{
			// combi_grid_dict aux_dict = create_aux_entire_dict(entire_downset);
			// inv_M = M_inv(aux_dict);

			// std::cout << "inverse of M " << std::endl;
			// for(int i = 0 ; i < size_downset ; ++i)
			// {
			// 	for(int j = 0 ; j < size_downset ; ++j)
			// 		std::cout << inv_M[i][j] << " ";

			// 	std::cout << std::endl;
			// }

			/* TO -DO : save the new coefficiants after the lp optimization problem is solved */
			std::vector<double> v;	

			return v;
		}

		virtual ~LP_OPT_INTERP()
		{
			free(constr_mat);
			free(row_index);
			free(col_index);

			glp_delete_prob(i_lp_prob);
		}
	};
}

#endif /* LPOPTINTERP_HPP_ */
