/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-30

Header file for GRB.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef SYNTHGRBLIBCLASS_H
#define SYNTHGRBLIBCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
// #include "Spectrum.hpp"
// #include "LightCurve.hpp"
#include "ModelParams.hpp"
#include "SynthGRB.hpp"

using namespace std;

// Declare class
class SynthGRBLibrary
{
public:
	// Class constructor
	SynthGRBLibrary(std::string path_name);
	
	// Class member variables
	std::vector<SynthGRB> synth_grb_list; // Vector of all loaded GRBs in the library directory
	std::string lib_dir_name; // Path to the library directory

	// Class member functions
	// Set the path name to the library directory 
	void set_SynthGRB_lib_dir(std::string path_name);
	
	// Make synthetic GRBs and save them into the specified library
	void MakeLibrary(float tw, float dte, 
		std::vector<float> eps_e_vec, std::vector<float> eps_b_vec, 
		std::vector<float> zeta_vec, std::vector<double> E_dot_iso_vec, 
		std::vector<float> theta_vec, std::vector<float> r_open_vec, 
		std::vector<float> eps_th_vec, std::vector<float> sigma_vec, std::vector<float> p_vec,
		std::string LorentzDist, std::string ShellDistParamsFile,
		float energ_min, float energ_max, int num_energ_bins, float tmin, float tmax);

	// Load SynthGRBs located in the directory specified by path_name
	void LoadLibrary(std::string path_name);
	// Load the SynthGRB contained in file designated by file_name 
	// (file must be within the currently set library directory)
	void LoadSynthGRB(std::string file_name);

	// Write library to FITS file in accordance to XPSEC template model FITS 
	int WriteLibraryToFITS(std::string out_file_name);
};

#endif 