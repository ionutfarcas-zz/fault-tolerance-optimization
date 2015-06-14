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

std::string python_code_caller(
    const std::string& script_name, 
    const int& level_min_x,
    const int& level_min_y, 
    const int& level_max_x,
    const int& level_max_y)
{
	std::stringstream caller;

	caller << "python " << script_name << " " << level_min_x << " " << level_min_y << " " << level_max_x << " " << level_max_y;

	return caller.str();
}

combi_grid_dict get_python_data(const std::string& script_run)
{
	FILE* stream;
    char buffer[256];
    std::string level_x_str;
    std::string level_y_str;
    std::string coeff_str;

    double level_x = 0.0;
    double level_y = 0.0;
    double coeff = 0.0;

    combi_grid_dict dict;   

    stream = popen(script_run.c_str(), "r");

    if(stream) 
    {
        while (!feof(stream))
        {
            if (fgets(buffer, sizeof(buffer), stream) != NULL)
            {
            	std::vector<int> levels;
                std::stringstream temp(buffer);
                temp >> level_x_str >> level_y_str >> coeff_str;

                level_x = str_to_number<int>(level_x_str);
                level_y = str_to_number<int>(level_y_str);
                coeff   = str_to_number<double>(coeff_str);
                levels.push_back(level_x);
                levels.push_back(level_y);

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

void set_M_matrix(double** M, const combi_grid_dict& aux_downset)
{
    int i = 0;
    int j = 0;

    std::vector<int> w;
    std::vector<int> c;

    for(auto ii = aux_downset.begin(); ii != aux_downset.end(); ++ii)
    {
        i = static_cast<int>(ii->second);
        j = 0;

        w = {ii->first[0], ii->first[1]};

        for(auto jj = aux_downset.begin(); jj != aux_downset.end(); ++jj)
        {
            j = static_cast<int>(jj->second);
            c = {jj->first[0], jj->first[1]};

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
}

void set_M_inv(double** M_inv, const int& size_downset)
{
    /* To be implemented */
}

combi_grid_dict entire_downset_dict(
    const int& level_max_x,
    const int& level_max_y,  
    const int& size_downset, 
    const std::string& script_run)
{
    int in_dict_size = 0;
    int in_out_diff = 0;
    double key = 0.0;
    std::vector<int> levels;
    combi_grid_dict input, output, result;

    input = get_python_data(script_run);
    in_dict_size = input.size();

    in_out_diff = size_downset - in_dict_size;

    if(in_out_diff != 0)
    {
        for(int i = 0 ; i < level_max_x ; ++i)
        {
            for(int j = 0 ; j < level_max_y - i ; ++j)
            {
                levels = {i + 1, j + 1};
                auto ii = input.find(levels);

                if(ii != input.end())
                {
                    key = ii->second;
                    output.insert(std::make_pair(levels, key));
                }
                else
                {
                    key = 0.0;
                    output.insert(std::make_pair(levels, key));
                }
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

combi_grid_dict aux_dict(const combi_grid_dict& entire_downset)
{
    double key = 0;
    int i = 0;
    std::vector<int> levels;

    combi_grid_dict aux_dict_out;

    for(auto ii = entire_downset.begin(); ii != entire_downset.end(); ++ii)
    {
        key = static_cast<double>(i);
        levels = {ii ->first[0], ii ->first[1]};
        ++i;

        aux_dict_out.insert(std::make_pair(levels, key));
    }

    return aux_dict_out;
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

int get_size_downset(const int& level_max_x, const int& level_max_y)
{
    assert(level_max_x == level_max_y);

    return static_cast<int>((level_max_x)*(level_max_x + 1)*0.5);
}

bool test_greater(const std::vector<int>& j, const std::vector<int>& i)
{
    bool test = false;

    int i_x = i[0];
    int i_y = i[1];

    int j_x = j[0];
    int j_y = j[1];

    test = ((j_x >= i_x) && (j_y >= i_y))?true:false;

    return test;
}