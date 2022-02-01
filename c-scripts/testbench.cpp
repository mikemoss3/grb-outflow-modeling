#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <tuple>
#include "fitsio.h"

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

int main(int argc, char const *argv[])
{

	// Array to store params 
	string inputs[13];
	int i = 0;

	// Load in the file
	string file_name = "input-files/jet-params.txt"; 
	ifstream file_jet_params(file_name);
	string line_jet_params;
	if ( file_jet_params.is_open() ) 
	{ 
		while(getline(file_jet_params, line_jet_params))
		{
			// If it is a commented line, just pass by
			if(line_jet_params[0] == '#'){ continue; }
			// Else: 
			inputs[i] = line_jet_params;
			++i;
		}	
	
		// Close files and free memory 
		file_jet_params.close();
	}
	else{ std::cout << "Unable to open file."; }


	for(int j=0; j<13; ++j)
	{
		cout << inputs[j] << endl;
	}


	return 0;
}
