#include "lp_opt.hpp"

namespace lp_opt
{
	LP_OPT::LP_OPT()
	{

	}

	LP_OPT::LP_OPT(const int& _i_n, const int& _i_m, const int& _opt_type)
	{
		assert(_i_n >= 1);
		assert(_i_m >= 1);
		assert(_opt_type == GLP_MIN || _opt_type == GLP_MAX);

		i_n = _i_n;
		i_m = _i_m;
		opt_type = _opt_type;

		rows = 2*_i_m + 1;
		cols = _i_n + _i_m;
		total_size = rows*cols;

		constr_mat = (double*)malloc((1 + total_size)*sizeof(double));
		row_index = (int*)malloc((1 + total_size)*sizeof(int));
		col_index = (int*)malloc((1 + total_size)*sizeof(int));

		assert(constr_mat!= NULL);
		assert(row_index!= NULL);
		assert(col_index!= NULL);
	}

	void LP_OPT::init_opti_prob()
	{
		std::string aux_var;

		i_lp_prob = glp_create_prob();
		assert(i_lp_prob != NULL);

		glp_set_prob_name(i_lp_prob, prob_name.c_str());
		glp_set_obj_dir(i_lp_prob, opt_type);

		glp_add_rows(i_lp_prob, rows);
		glp_add_cols(i_lp_prob, cols);

		for(int i = 0 ; i < rows - 1; ++i)
		{
			aux_var = set_aux_var_name("ineq_constr_", i + 1);
			glp_set_row_name(i_lp_prob, i + 1, aux_var.c_str());
			glp_set_row_bnds(i_lp_prob, i + 1, GLP_UP, 0.0, 0.0);
		}

		aux_var = set_aux_var_name("eq_constr_", rows);
		glp_set_row_name(i_lp_prob, rows, aux_var.c_str());
		glp_set_row_bnds(i_lp_prob, rows, GLP_FX, 1.0, 1.0);

		for(int i = 0 ; i < i_m ; ++i)
		{
			aux_var = set_aux_var_name("d", i + 1);
			glp_set_col_name(i_lp_prob, i + 1, aux_var.c_str());
			glp_set_col_bnds(i_lp_prob, i + 1, GLP_FR, 0.0, 0.0);
			glp_set_obj_coef(i_lp_prob, i + 1, 1.0);

		}
		for(int i = i_m ; i < cols ; ++i)
		{
			aux_var = set_aux_var_name("c", i + 1);
			glp_set_col_name(i_lp_prob, i + 1, aux_var.c_str());
			glp_set_col_bnds(i_lp_prob, i + 1, GLP_FR, 0.0, 0.0);
			glp_set_obj_coef(i_lp_prob, i + 1, 0.0);
		}
	}

	void LP_OPT::set_constr_matrix()
	{
		std::vector<double> rand_vars;
		rand_vars = gen_rand(i_m*i_n);

		for(int i = 0 ; i < i_m ; ++i)
		{
			for(int j = 0 ; j < i_n ; ++j)
			{
				constr_mat[j + i*cols + 1] = -rand_vars[j + i*i_n + 1];
			}

			for(int j = i_n ; j < cols ; ++j)
			{
				if(j - i_n == i)
				{
					constr_mat[j + i*cols + 1] = -1.0;
				}
			}	
		}

		for(int i = i_m ; i < rows - 1; ++i)
		{
			for(int j = 0 ; j < i_n ; ++j)
			{
				constr_mat[j + i*cols + 1] = rand_vars[j + (i - i_m)*i_n + 1];
			}

			for(int j = i_n ; j < cols ; ++j)
			{
				if(j - i_n == i - i_m)
				{
					constr_mat[j + i*cols + 1] = -1.0;
				}
			}	
		}

		for(int j = i_m ; j < cols ; ++j)
		{
			constr_mat[j + cols*(rows - 1) + 1] = 1.0;
		}

		for(int i = 0 ; i < rows ; ++i)
		{
			for(int j = 0 ; j < cols ; ++j)
			{
				row_index[j + i*cols + 1] = i + 1;
				col_index[j + i*cols + 1] = j + 1;
			}
		}

		std::cout << "Constraint matrix "<< std::endl;
		std::cout << "-----------" << std::endl;
		std::cout << "|  W | -I |" << std::endl;
		std::cout << "-----------" << std::endl;
		std::cout << "| -W | -I |" << std::endl;
		std::cout << "-----------" << std::endl;
		std::cout << "| 0m | 1n |" << std::endl;
		std::cout << "-----------" << std::endl;
		for(int i = 0 ; i < rows ; ++i)
		{
			for(int j = 0 ; j < cols ; ++j)
				std::cout << constr_mat[j + i*cols + 1] << " ";

			std::cout << std::endl;
		}
	}

	void LP_OPT::solve_opti_problem() const
	{
		glp_load_matrix(i_lp_prob, total_size, row_index, col_index, constr_mat);
		glp_simplex(i_lp_prob, NULL);	
	}

	std::vector<double> LP_OPT::get_results() const
	{
		double x = 0.0;		
		std::vector<double> result;

		for(int i = 0 ; i < cols ; ++i)
		{
			x = glp_get_col_prim(i_lp_prob, i + 1);
			result.push_back(x);
		}
		
		return result;
	}

	LP_OPT::~LP_OPT()
	{
		free(constr_mat);
		free(row_index);
		free(col_index);

		glp_delete_prob(i_lp_prob);
	}
}
