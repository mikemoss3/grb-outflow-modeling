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
#include "Spectrum.hpp"
// #include "cosmology.hpp"
#include "utilfuncs.hpp"
// #include "DataAnalysis.hpp"
// #include "TTEs.hpp"
// #include "ModelParams.hpp"
// #include "ShellDist.hpp"
#include "GammaFunction.hpp"

int main(int argc, char const *argv[])
{


	float emin = 1e-8;
	float emax = 1e3;
	float num_en_bins = 200;

	std::vector<float> energy_axis(num_en_bins);
	
	// Move to log space to define the energy interval with equally spaced points (in log space)
	float log_emin = log10(emin);
	float log_emax = log10(emax);
	float log_de = (log_emax - log_emin) / num_en_bins;

	// For each bin, calculate the energy axis
	for(int i=0;i<num_en_bins;++i)
	{
		energy_axis.at(i) = pow( 10., log_emin + (i*log_de));
	}


	double nu_1 = 5e+17;
	double nu_2 = 1e+15;
	float p = 2.5;
	double B = 1000.;
	float Gamma = 200.;


	ofstream spec_file; // Construct file 

	double param_list[6] = {nu_1, nu_2, p, B, Gamma, 1.}; // slow
	spec_file.open("test_synch_spec_slow.txt"); // Open text file with this name

	// double param_list[6] = {nu_2, nu_1, p, B, Gamma, 1.}; // fast
	// spec_file.open("test_synch_spec_fast.txt"); // Open text file with this name

	for(int i = 0; i < num_en_bins; ++i )
	{
		spec_file << energy_axis.at(i);
		spec_file << "\t";
		spec_file << Synchrotron(energy_axis.at(i),  param_list);
		spec_file << endl;		
	}

	
	return 0;
}

