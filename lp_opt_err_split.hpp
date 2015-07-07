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
		LP_OPT_ERR_SPLIT() 
		{
			i_n = 0;
			i_m = 0;
			rows = 0;
			cols = 0;
			total_size = 0;

			constr_mat = NULL;
			row_index = NULL;
			col_index = NULL;
		}

		LP_OPT_ERR_SPLIT(const int& _i_n, const int& _i_m, const int& _opt_type)
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

			i_lp_prob = glp_create_prob();
			assert(i_lp_prob != NULL);
		}

		LP_OPT_ERR_SPLIT(const LP_OPT_ERR_SPLIT& obj)
		{	
			i_n = obj.i_n;
			i_m = obj.i_m;
			opt_type = obj.opt_type;

			rows = obj.rows;
			cols = obj.cols;
			total_size = obj.total_size;

			constr_mat = (double*)malloc((1 + total_size)*sizeof(double));
			assert(constr_mat!= NULL);
			std::memcpy(constr_mat, obj.constr_mat, 1*sizeof(constr_mat));

			row_index = (int*)malloc((1 + total_size)*sizeof(int));
			assert(row_index!= NULL);
			std::memcpy(row_index, obj.row_index, 1*sizeof(row_index));

			col_index = (int*)malloc((1 + total_size)*sizeof(int));
			assert(col_index!= NULL);
			std::memcpy(col_index, obj.col_index, 1*sizeof(col_index));

			i_lp_prob = glp_create_prob();
			assert(i_lp_prob != NULL);
			std::memcpy(i_lp_prob, obj.i_lp_prob, 1*sizeof(i_lp_prob));
		}

		LP_OPT_ERR_SPLIT& operator= (const LP_OPT_ERR_SPLIT& rhs)
		{
			if(&rhs == this)
			{
				return *this;
			}

			i_n = rhs.i_n;
			i_m = rhs.i_m;
			opt_type = rhs.opt_type;

			rows = rhs.rows;
			cols = rhs.cols;
			total_size = rhs.total_size;

			constr_mat = (double*)malloc((1 + total_size)*sizeof(double));
			assert(constr_mat!= NULL);
			std::memcpy(constr_mat, rhs.constr_mat, 1*sizeof(constr_mat));

			row_index = (int*)malloc((1 + total_size)*sizeof(int));
			assert(row_index!= NULL);
			std::memcpy(row_index, rhs.row_index, 1*sizeof(row_index));

			col_index = (int*)malloc((1 + total_size)*sizeof(int));
			assert(col_index!= NULL);
			std::memcpy(col_index, rhs.col_index, 1*sizeof(col_index));

			i_lp_prob = glp_create_prob();
			assert(i_lp_prob != NULL);
			std::memcpy(i_lp_prob, rhs.i_lp_prob, 1*sizeof(i_lp_prob));

			return *this;
		}

		virtual void init_opti_prob(const std::string& prob_name)
		{
			std::string aux_var;

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

		virtual void set_constr_matrix()
		{
			/* TO DO: nothing here */
		}

		/* the constraint matrix has the form */
		/* ----------- */
		/* |  W | -I | */
 		/* ----------- */
 		/* | -W | -I | */
		/* ----------- */
		/* | 0m | 1n | */
		/* ----------- */
		virtual void set_constr_matrix(const std::vector<double>& W)
		{
			for(int i = 0 ; i < i_m ; ++i)
			{
				for(int j = 0 ; j < i_n ; ++j)
				{
					constr_mat[j + i*cols + 1] = -W[j + i*i_n + 1];
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
					constr_mat[j + i*cols + 1] = W[j + (i - i_m)*i_n + 1];
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

		virtual void solve_opti_problem() const
		{
			glp_load_matrix(i_lp_prob, total_size, row_index, col_index, constr_mat);
			glp_simplex(i_lp_prob, NULL);	
		}

		virtual std::vector<double> get_results() const
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

		virtual ~LP_OPT_ERR_SPLIT()
		{
			free(constr_mat);
			free(row_index);
			free(col_index);

			glp_delete_prob(i_lp_prob);
		}
	};
}

#endif /* LPOPTERRSPLIT_HPP_ */
