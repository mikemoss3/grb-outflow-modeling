#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "fitsio.h"

using namespace std;

#include "spectrum_funcs.hpp"
#include "Spectrum.hpp"
#include "response_funcs.hpp"
#include "Response.hpp"
// #include "light_curve_funcs.hpp"
// #include "cosmology.hpp"

int main(int argc, char const *argv[])
{

	char file_name[] = "/Users/mjmoss/Research/instr-rsp-investigate/example-instr-files/sw00883832000b_1speak.rsp";
	Response test_rsp = Response();
	test_rsp.load_rsp_from_file(file_name);

	// Make spectrum class
	Spectrum source_spectrum = Spectrum();
	source_spectrum.E_min = test_rsp.ENERG_LO.at(0);
	source_spectrum.E_max = test_rsp.ENERG_HI.back();
	source_spectrum.num_E_bins = test_rsp.num_phot_bins;


	// Fill in the energy vector
	for(int i=0; i < test_rsp.num_phot_bins; i++)
	{
		source_spectrum.ENERG_LO.push_back(test_rsp.ENERG_LO.at(i));
		source_spectrum.ENERG_MID.push_back(test_rsp.ENERG_MID.at(i));
		source_spectrum.ENERG_HI.push_back(test_rsp.ENERG_HI.at(i));		
	}

	// Zero out the spectrum (this also serves to set the length of the vector)
	source_spectrum.zero_spectrum();

	for(int i=0; i< source_spectrum.num_E_bins;i++)
	{
		source_spectrum.spectrum_dE.at(i) = 1000*(source_spectrum.ENERG_HI.at(i)-source_spectrum.ENERG_LO.at(i))*pow(source_spectrum.ENERG_MID.at(i),-1);
		source_spectrum.spectrum_sum += 1000*(source_spectrum.ENERG_HI.at(i)-source_spectrum.ENERG_LO.at(i))*pow(source_spectrum.ENERG_MID.at(i),-1);
	}


	Spectrum folded_spectrum = Spectrum();

	make_folded_spectrum(&folded_spectrum,test_rsp,source_spectrum);

	
	ofstream src_spec_dE_file;
	src_spec_dE_file.open("./test_src_spectrum_points.txt");
	int i=0;
	// For each energy bin, write the energy bin value and the spectrum rate to file.
	while ( i < source_spectrum.num_E_bins)
	{
		src_spec_dE_file << source_spectrum.ENERG_MID.at(i);
		src_spec_dE_file << " ";
		src_spec_dE_file << source_spectrum.spectrum_dE.at(i);
		src_spec_dE_file << "\n";		
		i++;
	}
	src_spec_dE_file.close(); // Close file

	ofstream folded_spec_dE_file;
	folded_spec_dE_file.open("./test_folded_spectrum_points.txt");

	i=0;
	// For each energy bin, write the energy bin value and the spectrum rate to file.
	while ( i < folded_spectrum.num_E_bins)
	{
		folded_spec_dE_file << folded_spectrum.ENERG_MID.at(i);
		folded_spec_dE_file << " ";
		folded_spec_dE_file << folded_spectrum.spectrum_dE.at(i);
		folded_spec_dE_file << "\n";		
		i++;
	}

	folded_spec_dE_file.close(); // Close file


	return 0;
}