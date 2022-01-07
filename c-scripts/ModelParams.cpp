/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-10-18

Model parameter class which contains all the necessary methods to interface with a set of fit parameters and statistics
*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "ModelParams.hpp"

ModelParams::ModelParams( float tw, float dte, float eps_e, float eps_b, float zeta, double E_dot_iso, float theta, float r_open, float eps_th, float sigma, std::string LorentzDist, std::string ShellDistParamsFile)
{
	// If no arguments are given, default values are assumed

	// Assign values for the jet parameters
	this->tw = tw; // sec, Duration of the wind
	this->dte = dte; // sec, Time between successive shell launches
	this->eps_e = eps_e; // Fraction of dissipated energy that goes into the electrons 
	this->eps_b = eps_b; // Fraction of dissipated energy that goes into the magnetic field 
	this->zeta = zeta; // Fraction of electrons which are accelerated 
	this->E_dot_iso = E_dot_iso; // erg/s, Injected isotropic energy rate
	this->theta = theta; // radians, Half-opening angle of the jet
	this->r_open = r_open; // cm, Opening radius of the jet
	this->eps_th = eps_th; // Fraction of energy in the outflow in the form of thermal energy 
	this->sigma = sigma; // Magnetization of the outflow 
	this->LorentzDist = LorentzDist; // Distribution of the jet shells
	this->ShellDistParamsFile = ShellDistParamsFile; // File that contains the parameters to create the distribution of jet shells. If Default is used, default params are passed

	this->numshells = tw/dte; // Number of shells in the jet

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void ModelParams::set_param_space(
		std::vector<float> eps_e_vec, 
		std::vector<float> eps_b_vec, 
		std::vector<float> zeta_vec, 
		std::vector<double> E_dot_iso_vec, 
		std::vector<float> theta_vec, 
		std::vector<float> r_open_vec, 
		std::vector<float> eps_th_vec, 
		std::vector<float> sigma_vec)
{
	this->eps_e_vec = eps_e_vec;
	this->eps_b_vec = eps_b_vec;
	this->zeta_vec = zeta_vec;
	this->E_dot_iso_vec = E_dot_iso_vec;
	this->theta_vec = theta_vec;
	this->r_open_vec = r_open_vec;
	this->eps_th_vec = eps_th_vec;
	this->sigma_vec = sigma_vec;
	this->LorentzDist = LorentzDist;
	this->ShellDistParamsFile = ShellDistParamsFile;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Copy a ModelParams object
void ModelParams::copy(ModelParams * p_model_params_in)
{
	this->tw = (*p_model_params_in).tw;
	this->dte = (*p_model_params_in).dte;
	this->numshells = (*p_model_params_in).numshells;
	this->eps_e = (*p_model_params_in).eps_e;
	this->eps_b = (*p_model_params_in).eps_b;
	this->zeta = (*p_model_params_in).zeta;
	this->E_dot_iso = (*p_model_params_in).E_dot_iso;
	this->theta = (*p_model_params_in).theta;
	this->r_open = (*p_model_params_in).r_open;
	this->eps_th = (*p_model_params_in).eps_th;
	this->sigma = (*p_model_params_in).sigma;
	this->LorentzDist = (*p_model_params_in).LorentzDist;
	this->ShellDistParamsFile = (*p_model_params_in).ShellDistParamsFile;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load ModelParams from a text file
void ModelParams::LoadFromTXT(std::string filename)
{
	// Load in the file
	ifstream file_model_param(filename);
	string line_model_params;
	if(file_model_param.is_open())
	{
		getline( file_model_param, line_model_params);
		this->tw = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->dte = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->eps_e = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->eps_b = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->zeta = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->E_dot_iso = stod(line_model_params);
		getline( file_model_param, line_model_params);
		this->theta = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->r_open = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->eps_th = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->sigma = stof(line_model_params);
		getline( file_model_param, line_model_params);
		this->LorentzDist = line_model_params;
		getline( file_model_param, line_model_params);
		this->ShellDistParamsFile = line_model_params;
		// Calculate number of jet shells
		this->numshells = tw/dte;

		// Close files and free memory 
		file_model_param.close();
	}
	else std::cout << "Unable to open file.";

}
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write ModelParams to a text file
void ModelParams::WriteToTXT(std::string filename)
{
	/*
	Write all current model parameter values in the ModelParams object to a text file with the file name specified by the user 
	*/

	// Construct file object
	ofstream out_file; 
	// Open text file with this name
	out_file.open(filename); 

	// Write parameters to file
	out_file << "tw = " << this->tw << "\n";
	out_file << "dte = " << this->dte << "\n";
	out_file << "numshells = " << this->numshells << "\n";
	out_file << "eps_e = " << this->eps_e << "\n";
	out_file << "eps_b = " << this->eps_b << "\n";
	out_file << "zeta = " << this->zeta << "\n";
	out_file << "E_dot_iso = " << this->E_dot_iso << "\n";
	out_file << "theta = " << this->theta << "\n";
	out_file << "r_open = " << this->r_open << "\n";
	out_file << "eps_th = " << this->eps_th << "\n";
	out_file << "sigma = " << this->sigma << "\n";
	out_file << "LorentzDist = " << this->LorentzDist << "\n";
	out_file << "ShellDistParamsFile = " << this->ShellDistParamsFile << "\n";

	out_file.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void ModelParams::PrintAllParams()
{
	std::cout << "tw = " << this->tw << "\n";
	std::cout << "dte = " << this->dte << "\n";
	std::cout << "numshells = " << this->numshells << "\n";
	std::cout << "eps_e = " << this->eps_e << "\n";
	std::cout << "eps_b = " << this->eps_b << "\n";
	std::cout << "zeta = " << this->zeta << "\n";
	std::cout << "E_dot_iso = " << this->E_dot_iso << "\n";
	std::cout << "theta = " << this->theta << "\n";
	std::cout << "r_open = " << this->r_open << "\n";
	std::cout << "eps_th = " << this->eps_th << "\n";
	std::cout << "sigma = " << this->sigma << "\n";
	std::cout << "LorentzDist = " << this->LorentzDist << "\n";
	std::cout << "ShellDistParamsFile = " << this->ShellDistParamsFile << "\n\n";
}
