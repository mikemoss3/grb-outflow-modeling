#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <tuple>
#include "fitsio.h"
#include <list>

using namespace std;

// Import Custom Libraries
// #include "SynthGRB.hpp"
// #include "SynthGRBLibrary.hpp"
// #include "ObsGRB.hpp"
// #include "Response.hpp"
// #include "LightCurve.hpp"
// #include "Spectrum.hpp"
// #include "cosmology.hpp"
// #include "utilfuncs.hpp"
// #include "DataAnalysis.hpp"
// #include "TTEs.hpp"
// #include "ModelParams.hpp"
// #include "ShellDist.hpp"

double BPL(float energy, double * param_list);

int main(int argc, char const *argv[])
{

	std::vector<std::vector<std::vector<double>>> model_param_bounds;
	model_param_bounds.push_back({{-2.,2.,5}});
	model_param_bounds.push_back({{1.,1e5,5},{-2.,2.,5},{-1.,-4.,5}});
	model_param_bounds.push_back({{1.,1e5,5},{-2.,2.,5},{-1.,-4.,5}});


	std::vector<std::vector<double>> param_combo_list;


	// Number of models
	int num_models = model_param_bounds.size();

	// Number of total parameters
	int num_params = 0;
	for(int i=0; i < num_models; ++i)
	{
		num_params += model_param_bounds.at(i).size();
	}

	std::vector<std::vector<double>> param_space_grid;
	// std::vector<std::vector<double>> param_space_grid[num_params];
	// int np = 0;

	// For each model, make a parameter space for each parameter
	for(int i=0; i < num_models; ++i)
	{
		for(size_t j = 0; j < model_param_bounds.at(i).size(); ++j)
		{
			
			std::vector<double> tmp_arr;
			float step_size = ( model_param_bounds.at(i).at(j).at(1) - model_param_bounds.at(i).at(j).at(0) ) / model_param_bounds.at(i).at(j).at(2);

			for(int k = 0; k <= model_param_bounds.at(i).at(j).at(2); ++k)
			{
				tmp_arr.push_back(model_param_bounds.at(i).at(j).at(0) + (k*step_size));
				cout << (model_param_bounds.at(i).at(j).at(0) + (k*step_size)) << " ";
			}
			cout <<endl;
			
			param_space_grid.push_back(tmp_arr);
			// ++np;
		}
	}	

	// Find the the number of parameter combinations
	int num_combos=1;
	for(size_t i=0; i < param_space_grid.size(); ++i)
	{
		num_combos *= param_space_grid.at(i).size();
	}
	// Resize the parameter combination list
	param_combo_list.resize(num_combos);

	cout << "num_combos = " << num_combos << endl;

	// To keep track of next element in each of the parameter
	int* indices = new int[num_params];
	// Initialize each to zero
	for (int i = 0; i < num_params; ++i)
	{
		indices[i] = 0;
	}
	

	// Find all parameter combinations
	// while(1)
	for(int k=0; k<num_combos; ++k)
	{
		// Current combination
		for (int i = 0; i < num_params; ++i)
		{
			param_combo_list.at(k).push_back( param_space_grid.at(i).at(indices[i]) );
			cout << param_space_grid.at(i).at(indices[i]) << " ";
		}
		cout << endl;

		// Find the "rightmost" parameter that has more elements left to use after this current element 
		int next = num_params - 1;
		while (next >= 0 && (indices[next] + 1 >= param_space_grid.at(next).size()))
		{
			--next;
		}

		// No such array is found so no more combinations left
		// if (next < 0)
		// 	return;

		// If found move to next element in that array
		++indices[next];

		// For all arrays to the right of this array current index again points to first element
		for (int i = next + 1; i < num_params; ++i)
		{
			indices[i] = 0;
		}
	}


	return 0;
}

