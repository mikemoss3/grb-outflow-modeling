#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <tuple>
#include <cfitsio/fitsio.h>
#include <list>

using namespace std;

// Import Custom Libraries
#include "SynthGRB.hpp"
// #include "SynthGRBLibrary.hpp"
// #include "ObsGRB.hpp"
// #include "Response.hpp"
// #include "LightCurve.hpp"
#include "Spectrum.hpp"
#include "cosmology.hpp"
#include "utilfuncs.hpp"
// #include "DataAnalysis.hpp"
// #include "TTEs.hpp"
// #include "ModelParams.hpp"
// #include "ShellDist.hpp"
// #include "GammaFunction.hpp"

using namespace std;

int main(int argc, char const *argv[])
{

	// float en = 150.;
	float en_c = 200./h_planck_kev;
	float en_m = 10./h_planck_kev;
	float p = 2.2;
	float B = 1e3;
	float gamma = 100;
	float Y = 1.;
	// float norm = 1.;
	
	
	double e_lo = 1e-3;
	double e_hi = 1e6;

	int norm_num_bin= 200.*log10(e_hi/e_lo); // Number of energy bins to use to calculate normalization
	std::vector<float> norm_energy_axis(norm_num_bin+1); // Make vector to store energy axis (initialized to zero's) 

	// Move to log space to define the energy interval with equally spaced points (in log space)
	float log_emin = log10(e_lo);
	float log_emax = log10(e_hi);
	float log_de = (log_emax - log_emin) / norm_num_bin;

	// For each bin, calculate the energy axis
	for(int i=0;i<norm_num_bin;++i)
	{
		norm_energy_axis.at(i) = pow( 10., log_emin + (i*log_de));
	}


	float en_curr=0.; // Current energy to evaluate the addition to the normalization
	double norm = 0.; // Set normalization to zero

	for(int i=0;i < norm_num_bin; ++i)
	{
		en_curr = (norm_energy_axis.at(i+1) + norm_energy_axis.at(i) ) / 2.; // Set the current energy
		// Calculate the contribution according to a Left-Riemann-Sum
		double param_list[7] = {en_c, en_m, p, B, gamma, Y, 1.};
		norm += (norm_energy_axis.at(i+1)-norm_energy_axis.at(i)) * Synchrotron(en_curr, param_list);
	}

	std::cout << 1./norm << endl;

	return 0;
}