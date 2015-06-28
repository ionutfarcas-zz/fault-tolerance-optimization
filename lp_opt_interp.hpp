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

		/* first level of grid indices */
		std::vector<int> i_level_1;
		/* second level of grid indices */
		std::vector<int> i_level_2;
		/* dimension of the problem */
		int i_dim;

		/* no. of constraints */
		int no_faults;
		/* level max sum */
		int l_max;

		/* total size of the optimization problem */
		int total_size;
		/* down set size */
		int size_downset;

		/* inverse of M s.t. w = Mc*/
		double** inv_M;

		/* given downset from python code */
		combi_grid_dict given_downset;
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
		LP_OPT_INTERP() {}

		LP_OPT_INTERP(
			std::vector<int> _level_1,
			std::vector<int> _level_2,  
			const int& _dim, 
			const int& _opt_type,
			const vec2d& _input_faults)
		{
			assert(_opt_type == GLP_MIN || _opt_type == GLP_MAX);

			i_level_1 = _level_1;
			i_level_2 = _level_2;
			opt_type = _opt_type;
			input_faults = _input_faults;

			size_downset = get_size_downset(_level_2);
			get_dict = python_code_caller(script_name, _level_1, _level_2);

			l_max = _level_1[0] + _level_2[1];

			valid_input_faults = filter_faults(input_faults, l_max, get_dict);
			no_faults = valid_input_faults.size();

			if(no_faults == 0)
			{
				std::cout << "Please introduce valid faults!" << std::endl;
				exit(0);
			}

			total_size = no_faults*size_downset;

			given_downset = get_python_data(get_dict);
			entire_downset = set_entire_downset_dict(_level_2, size_downset, get_dict);
		
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

			glp_add_rows(i_lp_prob, no_faults);
			glp_add_cols(i_lp_prob, size_downset);

			for(int i = 0 ; i < no_faults; ++i)
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
				glp_set_col_bnds(i_lp_prob, i + 1, GLP_DB, -1.0, 1.0);
				glp_set_obj_coef(i_lp_prob, i + 1, coeff);
			}
		}

		virtual void set_constr_matrix()
		{
			int inv_M_row_index = 0;
			std::vector<int> fault;

			for(int i = 0 ; i < no_faults ; ++i)
			{
				fault = {valid_input_faults[i][0], valid_input_faults[i][1]};
				auto it = aux_entire_dict.find(fault);

				if(it != aux_entire_dict.end())
				{
					inv_M_row_index = static_cast<int>(it->second);

					for(int j = 0 ; j < size_downset ; ++j)
					{
						constr_mat[j + i*size_downset + 1] = inv_M[inv_M_row_index][j];
					}
				}
			}

			for(int i = 0 ; i < no_faults ; ++i)
			{
				for(int j = 0 ; j < size_downset ; ++j)
				{
					row_index[j + i*size_downset + 1] = i + 1;
					col_index[j + i*size_downset + 1] = j + 1;
				}
			}
		}

		virtual void set_constr_matrix(const std::vector<double>& W)
		{	
			/* TO DO: nothing here */
		}

		virtual void solve_opti_problem() const
		{
			glp_load_matrix(i_lp_prob, total_size, row_index, col_index, constr_mat);
			glp_simplex(i_lp_prob, NULL);	
		}

		virtual std::vector<double> get_results() const
		{
			double w_i = 0.0;
			double c_i = 0.0;
			std::vector<double> w;		
			std::vector<double> c;

			combi_grid_dict input, output;

			for(int i = 0 ; i < size_downset ; ++i)
			{
				w_i = glp_get_col_prim(i_lp_prob, i + 1);
				w.push_back(w_i);
			}

			for(int i = 0 ; i < size_downset ; ++i)
			{
				c_i = 0.0;
				for(int j = 0 ; j < size_downset ; ++j)
				{
					c_i += inv_M[i][j]*w[j];
				}
				c.push_back(c_i);
			}

			input = get_python_data(get_dict);
			output = create_out_dict(given_downset, c);

			std::cout << std::endl;
			std::cout<< "Dictionary before optimization: " << std::endl;
			for(auto it = input.begin(); it != input.end(); ++it)
			{
				std::cout << "{(" << it->first[0] << ", " << it->first[1] << "), " << it->second << "} ";
			}
			std::cout << std::endl;

			std::cout << std::endl;
			std::cout << "Input faults" << std::endl;
			for(unsigned int i = 0 ; i < input_faults.size() ; ++i)
			{
				std::cout << "{" << input_faults[i][0] << ", " << input_faults[i][1] << "} ";
			}
			std::cout << std::endl;
			std::cout << "Valid input faults" << std::endl;
			for(int i = 0 ; i < no_faults ; ++i)
			{
				std::cout << "{" << valid_input_faults[i][0] << ", " << valid_input_faults[i][1] << "} ";
			}
			std::cout << std::endl;

			std::cout << std::endl;
			std::cout<< "Dictionary after optimization: " << std::endl;
			for(auto it = output.begin(); it != output.end(); ++it)
			{
				std::cout << "{(" << it->first[0] << ", " << it->first[1] << "), " << it->second << "} ";
			}
			std::cout << std::endl;
			
			return c;
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
