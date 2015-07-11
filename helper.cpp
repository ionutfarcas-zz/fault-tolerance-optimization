#include "helper.hpp"

template<typename T>
T str_to_number(const std::string& no)
{
	T value;
	std::stringstream stream(no);
	stream >> value;

	if (stream.fail()) 
	{
		std::runtime_error e(no);
        std::cout << "Error in the conversion of " << no << "!" << std::endl;
        throw e;
    }

    return value;
}

std::string python_code_caller(const std::string& script_name, const vec2d& levels)
{
    int levels_no = 0;
    int level_size = 0;
    int one_level = 0;
    std::stringstream caller;

    levels_no = levels.size();
    level_size = levels[0].size();

    caller << "python " << script_name << " ";

    for(int i = 0 ; i < levels_no ; ++i)
    {
        for(int j = 0 ; j < level_size ; ++j)
        {
            one_level = levels[i][j];
            caller << one_level << " ";
        }
    }

    return caller.str();
}

combi_grid_dict get_python_data(const std::string& script_run, const int& dim)
{
	FILE* stream;
    char buffer[256];
    std::string level_x_str;
    std::string level_y_str;
    std::string coeff_str;

    double coeff = 0.0;

    combi_grid_dict dict;   

    stream = popen(script_run.c_str(), "r");

    if(stream) 
    {
        while (!feof(stream))
        {
            if (fgets(buffer, sizeof(buffer), stream) != NULL)
            {
                std::string one_level_str;
                int one_level = 0;
                std::vector<int> levels;
                std::stringstream temp(buffer);

                for(int i = 0 ; i < dim ; ++i)
                {
                    temp >> one_level_str;
                    one_level = str_to_number<int>(one_level_str);
                    levels.push_back(one_level);
                }

                temp >> coeff_str;
                coeff = str_to_number<double>(coeff_str);

                dict.insert(std::make_pair(levels, coeff));
            }
        }

        pclose(stream);
    }
    else
    {
        throw "Error reading script output!"; 
    }

    return dict;
}

double** M_matrix(const combi_grid_dict& aux_downset, const int& dim)
{
    int size_downset = aux_downset.size();
    int i = 0;
    int j = 0;

    double** M = (double**)calloc(size_downset, sizeof(double*));
    for(int i = 0 ; i < size_downset ; ++i)
    {
        M[i] = (double*)calloc(size_downset, sizeof(double));
    }

    assert(M!= NULL);

    for(auto ii = aux_downset.begin(); ii != aux_downset.end(); ++ii)
    {
        i = static_cast<int>(ii->second);
        j = 0;

        std::vector<int> w;
        for(int it = 0 ; it < dim ; ++it)
        {
            w.push_back(ii->first[it]);
        }

        for(auto jj = aux_downset.begin(); jj != aux_downset.end(); ++jj)
        {
            std::vector<int> c;
            for(int it = 0 ; it < dim ; ++it)
            {
                c.push_back(jj->first[it]);
            }
            j = static_cast<int>(jj->second);

            if(test_greater(c, w))
            {
                M[i][j] = 1.0;
            }
            else
            {
                M[i][j] = 0.0;   
            }
        }
    }

    return M;
}

double** N_matrix(const combi_grid_dict& aux_downset, const int& dim)
{
    int size_downset = aux_downset.size();
    int i = 0;
    int j = 0;

    double** N = (double**)calloc(size_downset, sizeof(double*));
    for(int i = 0 ; i < size_downset ; ++i)
    {
        N[i] = (double*)calloc(size_downset, sizeof(double));
    }

    assert(N!= NULL);

    for(auto ii = aux_downset.begin(); ii != aux_downset.end(); ++ii)
    {
        i = static_cast<int>(ii->second);
        j = 0;

        std::vector<int> w;
        for(int it = 0 ; it < dim ; ++it)
        {
            w.push_back(ii->first[it]);
        }

        for(auto jj = aux_downset.begin(); jj != aux_downset.end(); ++jj)
        {
            std::vector<int> c;
            for(int it = 0 ; it < dim ; ++it)
            {
                c.push_back(jj->first[it]);
            }
            j = static_cast<int>(jj->second);

            if(i == j)
            {
                N[i][j] = 0.0;
            }
            else if(test_greater(c, w))
            {
                N[i][j] = 1.0;
            }
            else
            {
                N[i][j] = 0.0;   
            }
        }
    }

    return N;
}

double** mat_prod(double** A, double** B, const int& dim)
{
    double** result = (double**)calloc(dim, sizeof(double*));
    for(int i = 0 ; i < dim ; ++i)
    {
        result[i] = (double*)calloc(dim, sizeof(double));
    }

    for(int i = 0 ; i < dim ; ++i)
    {
        for(int j = 0 ; j < dim ; ++j)
        {
            for(int k = 0 ; k < dim ; ++k)
            {
                result[i][j] += A[i][k]*B[k][j];
            }
        }
    }

    return result;
}

double** N_pow_k(double** N, const int& size_downset, const int& k)
{
    if(k == 1)
    {
        return N;
    }
    else
    {
        return (mat_prod(N, N_pow_k(N, size_downset, k-1), size_downset));
    }
}

double** sum_pow_N(double** N, const int& size_downset)
{
    int sign = 0;

    double** N_to_k = (double**)calloc(size_downset, sizeof(double*));
    for(int i = 0 ; i < size_downset ; ++i)
    {
        N_to_k[i] = (double*)calloc(size_downset, sizeof(double));
    }

    double** result = (double**)calloc(size_downset, sizeof(double*));
    for(int i = 0 ; i < size_downset ; ++i)
    {
        result[i] = (double*)calloc(size_downset, sizeof(double));
    }

    assert(N_to_k!= NULL);
    assert(result!= NULL);

    for (int k = 1 ; k <= size_downset - 1 ; ++k)
    {
        sign = pow(-1, k);
        N_to_k = N_pow_k(N, size_downset, k);

        for(int i = 0 ; i < size_downset ; ++i)
        {
            for(int j = i ; j < size_downset ; ++j)
            {
                result[i][j] += sign*N_to_k[i][j];
            }
        }
    }

    return result;
}

double** M_inv(const combi_grid_dict& aux_downset, const int& dim)
{
    int size_downset = aux_downset.size();

    double** N = (double**)calloc(size_downset, sizeof(double*));
    for(int i = 0 ; i < size_downset ; ++i)
    {
        N[i] = (double*)calloc(size_downset, sizeof(double));
    }

    double** sum_N = (double**)calloc(size_downset, sizeof(double*));
    for(int i = 0 ; i < size_downset ; ++i)
    {
        sum_N[i] = (double*)calloc(size_downset, sizeof(double));
    }

    double** M_inv = (double**)calloc(size_downset, sizeof(double*));
    for(int i = 0 ; i < size_downset ; ++i)
    {
        M_inv[i] = (double*)calloc(size_downset, sizeof(double));
    }

    assert(N!= NULL);
    assert(sum_N!= NULL);
    assert(M_inv!= NULL);

    N = N_matrix(aux_downset, dim);
    sum_N = sum_pow_N(N, size_downset);

    for(int i = 0 ; i < size_downset ; ++i)
    {
        for(int j = i ; j < size_downset ; ++j)
        {
            if(i == j)
            {
                M_inv[i][j] = 1.0;
            }
            else
            {
                M_inv[i][j] = sum_N[i][j];
            }
        }
    }

    return M_inv;
}

combi_grid_dict set_entire_downset_dict(
    const std::vector<int>& level_max, 
    const int& size_downset, 
    const std::string& script_run,
    const int& dim)
{
    int size = 1;
    int in_dict_size = 0;
    int in_out_diff = 0;
    double key = 0.0;
    std::vector<int> level;
    vec2d levels;
    combi_grid_dict input, output, result;

    input = get_python_data(script_run, dim);
    in_dict_size = input.size();

    in_out_diff = size_downset - in_dict_size;
    auto max_level_max = std::max_element(level_max.begin(), level_max.end());

    levels = mindex(dim, *max_level_max);
    size = levels.size();

    if(in_out_diff != 0)
    {
        for(int i = 0 ; i < size ; ++i)
        {
            level = levels[i];
            auto ii = input.find(level);

            if(ii != input.end())
            {
                key = ii->second;
                output.insert(std::make_pair(level, key));
            }
            else
            {
                key = 0.0;
                output.insert(std::make_pair(level, key));
            }
            
        }
        result = output;
    }
    else
    {
        result = input;
    }

    return result;
}

combi_grid_dict create_aux_entire_dict(const combi_grid_dict& entire_downset, const int& dim)
{
    double key = 0;
    int i = 0;

    combi_grid_dict aux_dict;

    for(auto ii = entire_downset.begin(); ii != entire_downset.end(); ++ii)
    {
        std::vector<int> levels;
        key = static_cast<double>(i);

        for(int i = 0 ; i < dim ; ++i)
        {
            levels.push_back(ii->first[i]);
        }

        ++i;
        aux_dict.insert(std::make_pair(levels, key));
    }

    return aux_dict;
}

vec2d get_downset_indices(const combi_grid_dict& entire_downset, const int& dim)
{
    vec2d indices;
    
    for(auto ii = entire_downset.begin(); ii != entire_downset.end(); ++ii)
    {
        std::vector<int> index;

        for(int i = 0 ; i < dim ; ++i)
        {
            index.push_back(ii->first[i]);
        }

        indices.push_back(index);
    }

    return indices;
}

vec2d filter_faults(
    const vec2d& faults_input, 
    const int& l_max, 
    const std::string& script_run,
    const int& dim)
{
    int no_faults = 0;
    int level_fault = 0;
    combi_grid_dict received_dict;
    
    vec2d faults_output;

    no_faults = faults_input.size();
    received_dict = get_python_data(script_run, dim);
    
    for(int i = 0 ; i < no_faults ; ++i)
    {
        auto it = received_dict.find(faults_input[i]);

        if(it != received_dict.end())
        {
            level_fault = std::accumulate(faults_input[i].begin(), faults_input[i].end(), 0);

            if((level_fault == l_max) || (level_fault == (l_max - 1)))
            {
                faults_output.push_back(faults_input[i]);
            }
        }
    }

    return faults_output;
}

combi_grid_dict create_out_dict(const combi_grid_dict& given_downset, const std::vector<double>& new_c, const int& dim)
{
    double key = 0;
    int i = 0;

    combi_grid_dict out_dict;

    for(auto ii = given_downset.begin(); ii != given_downset.end(); ++ii)
    {
        std::vector<int> levels;
        key = new_c[i];

        for(int i = 0 ; i < dim ; ++i)
        {
            levels.push_back(ii->first[i]);
        }

        ++i;
        out_dict.insert(std::make_pair(levels, key));
    }

    return out_dict;
}

std::string set_aux_var_name(const std::string& var_name, const int& index)
{
	std::stringstream aux_var;
	aux_var << var_name << index;

	return aux_var.str();
}

std::vector<double> gen_rand(const int& size)
{
	double rand_var = 0.0;
	std::vector<double> output;

	for(int i = 0 ; i < size ; ++i)
	{
     rand_var = 1e-2*(std::rand()%10);
     output.push_back(rand_var);
 }

 return output;
}

int get_size_downset(const std::vector<int>& level_max, const int& dim)
{
    int size = 1;
    int min_level_max = *std::min_element(level_max.begin(), level_max.end());

    for(int i = 0 ; i < dim ; ++i)
    {
        size *= (min_level_max + i);
    }

    size = static_cast<int>(size/(factorial(dim)));

    return size;
}

int l1_norm(const std::vector<int>& u)
{
    int norm = 0;

    for (int elem : u)
    { 
        norm += abs(elem); 
    }

    return norm;
}

int factorial(const int& dim)
{
    int fact = 0;

    if(dim == 0 || dim == 1)
    {
        fact = 1;
    }
    else
    {
        fact = dim*factorial(dim - 1);
    }

    return fact;
}

bool test_greater(const std::vector<int>& b, const std::vector<int>& a)
{
    int dim = a.size();
    bool test = true;


    for(int i = 0 ; i < dim ; ++i)
    {
        test *= (b[i] >= a[i])?true:false;
    }
    
    return test;
}

vec2d mindex(const int& dimension, const int& upper_limit)
{
    int j = 0;

    std::vector<int> temp(dimension, 1);
    std::vector<int> ones;
    vec2d mindex_result;

    for(int i = 0 ; i < dimension ; ++i)
    {
        ones.push_back(1);
    }

    while(true)
    {
        mindex_result.push_back(temp);

        for(j = dimension - 1 ; j >= 0 ; --j)
        {
            if(++temp[j] <= upper_limit)
                break;
            else
                temp[j] = 1;
        }

        if( j < 0)
            break;
    }

    return mindex_result;
}

void check_input_levels(const vec2d& levels)
{
    std::vector<int> l_min = levels[0];
    std::vector<int> l_max = levels[1];
    std::vector<int> c;

    for(unsigned int i = 0 ; i < l_min.size() ; ++i)
    {
        c.push_back(l_max[i] - l_min[i]);
    }

    if(std::adjacent_find(c.begin(), c.end(), std::not_equal_to<int>() ) == c.end())
    {
        std::cout << "Correct input levels!" << std::endl;
    }
    else
    {
        std::cout << "Input levels are incorrect!" << std::endl;
        std::cout << "Please input them of the form: l_max = l_min + c*ones(dim), c>=1, integer" << std::endl;
        exit(0);
    }

}