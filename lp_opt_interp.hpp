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
		std::string python_caller;

		/* maximum level x direction*/
		int i_level_min_x;
		/* maximum level y direction*/
		int i_level_min_y;
		/* maximum level x direction*/
		int i_level_max_x;
		/* maximum level y direction*/
		int i_level_max_y;

		/* total size of the optimization problem */
		int total_size;
		/* down set size */
		int size_downset;

		/* M matrix s.t. w = Mc */
		double** M;
		/* inverse of M */
		double** M_inv;

		/* entire donwset with corresponding indices */
		combi_grid_dict entire_downset;

	public:
		LP_OPT_INTERP() {}

		LP_OPT_INTERP(
			const int& _i_level_min_x,
			const int& _i_level_min_y,  
			const int& _i_level_max_x,
			const int& _i_level_max_y, 
			const int& _opt_type)
		{
			assert(_i_level_min_x >= 1);
			assert(_i_level_min_y >= 1);
			assert(_i_level_max_x >= 1);
			assert(_i_level_max_y >= 1);

			assert(_i_level_max_x >= _i_level_min_x);
			assert(_i_level_max_y >= _i_level_min_y);

			assert(_opt_type == GLP_MIN || _opt_type == GLP_MAX);

			i_level_min_x = _i_level_min_x;
			i_level_min_y = _i_level_min_y;
			i_level_max_x = _i_level_max_x;
			i_level_max_y = _i_level_max_y;

			opt_type = _opt_type;

			size_downset = get_size_downset(_i_level_max_x, _i_level_max_y);
			total_size = size_downset*size_downset;

			python_caller = python_code_caller(
				script_name, 
				i_level_min_x, 
				i_level_min_y,  
				i_level_max_x, 
				i_level_max_y);

			entire_downset = entire_downset_dict(
				i_level_max_x,
				i_level_max_y,  
				size_downset, 
				python_caller);

			M = (double**)calloc(size_downset, sizeof(double*));
			for(int i = 0 ; i < size_downset ; ++i)
			{
				M[i] = (double*)calloc(size_downset, sizeof(double));
			}

			M_inv = (double**)calloc(size_downset, sizeof(double*));
			for(int i = 0 ; i < size_downset ; ++i)
			{
				M_inv[i] = (double*)calloc(size_downset, sizeof(double));
			}

			constr_mat = (double*)malloc((1 + total_size)*sizeof(double));
			row_index = (int*)malloc((1 + total_size)*sizeof(int));
			col_index = (int*)malloc((1 + total_size)*sizeof(int));

			assert(constr_mat!= NULL);
			assert(row_index!= NULL);
			assert(col_index!= NULL);
		}

		virtual void init_opti_prob(const std::string& prob_name)
		{
			i_lp_prob = glp_create_prob();
			assert(i_lp_prob != NULL);

			glp_set_prob_name(i_lp_prob, prob_name.c_str());
			glp_set_obj_dir(i_lp_prob, opt_type);

			std::cout << "Levels(x, y) and coefficients corresponding to the entire downset " << std::endl;
			for(auto ii = entire_downset.begin(); ii != entire_downset.end(); ++ii)
			{
				std::cout << ii->first[0] << " " << ii->first[1] << " " << ii->second << std::endl;
			}

			combi_grid_dict aux = aux_dict(entire_downset);
			
			set_M_matrix(M, aux);
			std::cout << "M matrix: " << std::endl;
			for(int i = 0 ; i < size_downset ; ++i)
			{
				for(int j = 0 ; j < size_downset ; ++j)
					std::cout << M[i][j] << " ";

				std::cout << std::endl;
			}
		}

		virtual void set_constr_matrix()
		{
			/* TO -DO : inverse of the M matrix, obtained from the mapping of levels and coefficients */
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
			/* TO -DO : save the new coefficiants after the lp optimization problem is solved */
			std::vector<double> v;	

			return v;
		}

		virtual ~LP_OPT_INTERP()
		{
			for(int i = 0 ; i < size_downset; ++i)
    		{
    			free(M[i]);
    		}
    		free(M);
    		for(int i = 0 ; i < size_downset; ++i)
    		{
    			free(M_inv[i]);
    		}
    		free(M_inv);

			free(constr_mat);
			free(row_index);
			free(col_index);

			glp_delete_prob(i_lp_prob);
		}
	};
}

#endif /* LPOPTINTERP_HPP_ */
