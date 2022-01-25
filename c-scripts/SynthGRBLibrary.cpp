/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-30

SynthGRBLibrary class which will contain all higher level functions to interact with a library of SynthGRB objects/files. 
*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <filesystem>

// Import Custom Libraries
#include "ModelParams.hpp"
#include "SynthGRB.hpp"
#include "SynthGRBLibrary.hpp"

using namespace std;

// Class constructor
SynthGRBLibrary::SynthGRBLibrary(std::string path_name)
{
	set_SynthGRB_lib_dir(path_name);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the path name to the directory containing the SynthGRB files.
void SynthGRBLibrary::set_SynthGRB_lib_dir(std::string path_name)
{
	lib_dir_name = path_name;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make synthetic GRBs and save them into the specified library
void SynthGRBLibrary::MakeLibrary(float tw, float dte, 
	std::vector<float> eps_e_vec, std::vector<float> eps_b_vec, 
	std::vector<float> zeta_vec, std::vector<double> E_dot_iso_vec, 
	std::vector<float> theta_vec, std::vector<float> r_open_vec, 
	std::vector<float> eps_th_vec, std::vector<float> sigma_vec, std::vector<float> p_vec,
	std::string LorentzDist, std::string ShellDistParamsFile,
	float energ_min, float energ_max, int num_energ_bins, float tmin, float tmax)
{
	// Input Parameters: 
	// tw = 10; // sec, Duration of the wind
	// dte = 0.002; // sec, Time between successive shell launches
	// eps_e_vec = Vector which stores the fraction of dissipated energy that goes into the electrons 
	// eps_b_vec = Vector which stores the fraction of dissipated energy that goes into the magnetic field 
	// zeta_vec = Vector which stores the fraction of electrons which are accelerated 
	// E_dot_iso_vec = erg/s, Vector which stores the  injected isotropic energy rate
	// theta_vec = radians, Vector which stores the half-opening angle of the jet
	// r_open_vec = cm, Vector which stores the opening radius of the jet
	// eps_th_vec = Vector which stores the fraction of energy in the outflow in the form of thermal energy 
	// sigma_vec = Vector which stores the magnetization of the outflow 
	// LorentzDist = Distribution of the jet shells
	// ShellDistParamsFile = "Default";

	// energ_min = Minimum energy of the spectra
	// energ_max = Maximum energy of the spectra
	// num_energ_bins = number of energy bins in the spectra
	// tmin = Minimum time of the emission to calculate the spectra across
	// tmax = Maximum time of the emission to calculate the spectra across	

	// Calculate total number of SynthGRBs in the library
	int num_synthGRBs = eps_e_vec.size() * eps_b_vec.size() * zeta_vec.size() * E_dot_iso_vec.size() * theta_vec.size() * r_open_vec.size() * eps_th_vec.size() * sigma_vec.size();

	// Set size of library
	synth_grb_list.resize(num_synthGRBs);
	int iter = 0;

	// Using the range of input parameters, make and save Synthetic GRBs for the entire parameter space
	for(size_t i1=0; i1< eps_e_vec.size(); ++i1)
	{
		for(size_t i3=0; i3< eps_b_vec.size(); ++i3)
		{
			for(size_t i4=0; i4< zeta_vec.size(); ++i4)
			{
				for(size_t i6=0; i6< E_dot_iso_vec.size(); ++i6)
				{
					for(size_t i7=0; i7< theta_vec.size(); ++i7)
					{
						for(size_t i8=0; i8< r_open_vec.size(); ++i8)
						{
							for(size_t i10=0; i10< eps_th_vec.size(); ++i10)
							{
								for(size_t i11=0; i11< sigma_vec.size(); ++i11)
								{
									for(size_t i12=0; i12 < p_vec.size(); ++i12)
									{
										/*
										SynthGRB * template_grb = new SynthGRB();
										(*template_grb).tw = tw;
										(*template_grb).dte = dte;
										(*template_grb).eps_e = eps_e_vec.at(i1);
										(*template_grb).eps_b = eps_b_vec.at(i3);
										(*template_grb).zeta = zeta_vec.at(i4);
										(*template_grb).E_dot_iso = E_dot_iso_vec.at(i6);
										(*template_grb).theta = theta_vec.at(i7);
										(*template_grb).r_open = r_open_vec.at(i8);
										(*template_grb).eps_th = eps_th_vec.at(i10);
										(*template_grb).sigma = sigma_vec.at(i11);
										(*template_grb).LorentzDist = LorentzDist;
										(*template_grb).ShellDistParamsFile = ShellDistParamsFile;
										*/
										SynthGRB * template_grb = new SynthGRB(tw,
											dte,
											eps_e_vec.at(i1),
											eps_b_vec.at(i3),
											zeta_vec.at(i4),
											E_dot_iso_vec.at(i6),
											theta_vec.at(i7),
											r_open_vec.at(i8),
											eps_th_vec.at(i10),
											sigma_vec.at(i11),
											p_vec.at(i12),
											LorentzDist,
											ShellDistParamsFile);
										
										(*template_grb).SimulateJetDynamics();
										(*template_grb).make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax);
										synth_grb_list.at(iter) = (*template_grb);
										iter+=1;
									}

								}
							}
						}
					}
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load SynthGRBs located in the directory specified by path_name
void SynthGRBLibrary::LoadLibrary(std::string path_name)
{
	for (const auto & entry : std::filesystem::directory_iterator(path_name))
	{
		LoadSynthGRB(entry.path());
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load the SynthGRB contained in file designated by file_name 
// (file must be within the currently set library directory)
void SynthGRBLibrary::LoadSynthGRB(std::string file_name)
{
	SynthGRB another_SynthGRB;
	another_SynthGRB.LoadFromFile(file_name);

	// Add SynthGRB to library list
	synth_grb_list.push_back(another_SynthGRB);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write library to FITS file in accordance to XPSEC template model FITS format
int SynthGRBLibrary::WriteLibraryToFITS(std::string out_file_name)
{
	// Declare all variable types 
	fitsfile *fptr; // pointer to FITS file object 
	int status = 0; // Error status indicator 

	// Declaring character array
	std::string s = "!"+out_file_name;
    int n = s.length();
    char char_array[n + 1];
 
    // Copying the contents of the string to char array
    strcpy(char_array, s.c_str());

	// Open FITS file
	// fits_create_file(fitsfile **fptr, char *filename, > int *status)
	fits_create_file(&fptr, char_array, &status);

	// Define number of HDUs
	// For each HDU: 
	// Define Keywords
	// Define Data Table  

	










	/*
	// Get number of HDUs 
	fits_get_num_hdus(fptr, &num_hdus, &status);
	// Skip hdu=1, this is the primary HDU
	for (ihdu=2; ihdu<=num_hdus; ihdu++)
	{
		// Move to current HDU 
		fits_movabs_hdu(fptr, ihdu, &hdutype, &status);
		// We only want Binary Table data units
		if (hdutype == BINARY_TBL) 
		{
			// Extract the extension name and the array name.
			fits_read_key(fptr, TSTRING, "EXTNAME", extname, comment, &status);
			
			// Check for errors 
			if (status){fits_report_error(stderr, status); return status;}

			// Keep track of which HDU corresponds to which extension 
			if (strcmp(extname, "SPECTRUM")==0){spectrum_ihdu = ihdu;}
			if (strcmp(extname, "EBOUNDS")==0){ebounds_ihdu = ihdu;}
		}
	}
	*/

	fits_close_file(fptr, &status);
	return status;
}