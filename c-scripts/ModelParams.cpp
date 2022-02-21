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

ModelParams::ModelParams(float tw, float dte, 
		double E_dot_iso, float theta, float r_open, float eps_th, float sigma,
		float eps_e_int, float eps_b_int, float zeta_int, float p_int,
		float eps_e_ext, float eps_b_ext, float zeta_ext, float p_ext, 
		float k_med, double rho_not,
		std::string LorentzDist, std::string ShellDistParamsFile)
{
	// If no arguments are given, default values are assumed

	// Assign values for the jet parameters
	this->tw = tw; // sec, Duration of the wind
	this->dte = dte; // sec, Time between successive shell launches
	this->E_dot_iso = E_dot_iso; // erg/s, Injected isotropic energy rate
	this->theta = theta; // radians, Half-opening angle of the jet
	this->r_open = r_open; // cm, Opening radius of the jet
	this->eps_th = eps_th; // Fraction of energy in the outflow in the form of thermal energy 
	this->sigma = sigma; // Magnetization of the outflow 

	this->eps_e_int = eps_e_int; // Fraction of dissipated energy that goes into the electrons (in internal shocks)
	this->eps_b_int = eps_b_int; // Fraction of dissipated energy that goes into the magnetic field (in internal shocks)
	this->zeta_int = zeta_int; // Fraction of electrons which are accelerated in the shock (in internal shocks)
	this->p_int = p_int; // Power law index of the electron population (in internal shocks)

	this->eps_e_ext = eps_e_ext; // Fraction of dissipated energy that goes into the electrons  (in the forward shock)
	this->eps_b_ext = eps_b_ext; // Fraction of dissipated energy that goes into the magnetic field (in the forward shock) 
	this->zeta_ext = zeta_ext; // Fraction of electrons which are accelerated in the shock (in the forward shock)
	this->p_ext = p_ext; // Power law index of the electron population  (in the forward shock)

	this->k_med = k_med; // Indication of a wind or constant medium
	this->rho_not = rho_not; // Normalization of the external medium 

	this->LorentzDist = LorentzDist; // Distribution of the jet shells
	this->ShellDistParamsFile = ShellDistParamsFile; // File that contains the parameters to create the distribution of jet shells. If Default is used, default params are passed

	this->numshells = tw/dte; // Number of shells in the jet

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void ModelParams::set_param_space(
	std::vector<double> E_dot_iso_vec,
	std::vector<float> theta_vec,
	std::vector<float> r_open_vec,
	std::vector<float> eps_th_vec,
	std::vector<float> sigma_vec,
	std::vector<float> eps_e_int_vec,
	std::vector<float> eps_b_int_vec,
	std::vector<float> zeta_int_vec,
	std::vector<float> p_int_vec,
	std::vector<float> eps_e_ext_vec,
	std::vector<float> eps_b_ext_vec,
	std::vector<float> zeta_ext_vec,
	std::vector<float> p_ext_vec,
	std::vector<float> k_med_vec,
	std::vector<double> rho_not_vec)
{
	this->E_dot_iso_vec = E_dot_iso_vec;
	this->theta_vec = theta_vec;
	this->r_open_vec = r_open_vec;
	this->eps_th_vec = eps_th_vec;
	this->sigma_vec = sigma_vec;
	
	this->eps_e_int_vec = eps_e_int_vec;
	this->eps_b_int_vec = eps_b_int_vec;
	this->zeta_int_vec = zeta_int_vec;
	this->p_int_vec = p_int_vec;

	this->eps_e_ext_vec = eps_e_ext_vec;
	this->eps_b_ext_vec = eps_b_ext_vec;
	this->zeta_ext_vec = zeta_ext_vec;
	this->p_ext_vec = p_ext_vec;

	this->k_med_vec = k_med_vec;
	this->rho_not_vec = rho_not_vec;

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
	this->E_dot_iso = (*p_model_params_in).E_dot_iso;
	this->theta = (*p_model_params_in).theta;
	this->r_open = (*p_model_params_in).r_open;
	this->eps_th = (*p_model_params_in).eps_th;
	this->sigma = (*p_model_params_in).sigma;

	this->eps_e_int = (*p_model_params_in).eps_e_int;
	this->eps_b_int = (*p_model_params_in).eps_b_int;
	this->zeta_int = (*p_model_params_in).zeta_int;
	this->p_int = (*p_model_params_in).p_int;

	this->eps_e_ext = (*p_model_params_in).eps_e_ext;
	this->eps_b_ext = (*p_model_params_in).eps_b_ext;
	this->zeta_ext = (*p_model_params_in).zeta_ext;
	this->p_ext = (*p_model_params_in).p_ext;
	
	this->k_med = (*p_model_params_in).k_med;
	this->rho_not = (*p_model_params_in).rho_not;

	this->LorentzDist = (*p_model_params_in).LorentzDist;
	this->ShellDistParamsFile = (*p_model_params_in).ShellDistParamsFile;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load ModelParams from a text file
void ModelParams::LoadFromTXT(std::string filename)
{
	// Array to store params 
	string inputs[19];
	int i = 0;

	// Load in the file
	ifstream file_jet_params(filename);
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

	this->tw = stof(inputs[0]);
	this->dte = stof(inputs[1]);

	this->E_dot_iso = stod(inputs[2]);
	this->theta = stof(inputs[3]);
	this->r_open = stof(inputs[4]);
	this->eps_th = stof(inputs[5]);
	this->sigma = stof(inputs[6]);

	this->eps_e_int = stof(inputs[7]);
	this->eps_b_int = stof(inputs[8]);
	this->zeta_int = stof(inputs[9]);
	this->p_int = stof(inputs[10]);

	this->eps_e_ext = stof(inputs[11]);
	this->eps_b_ext = stof(inputs[12]);
	this->zeta_ext = stof(inputs[13]);
	this->p_ext = stof(inputs[14]);
	
	this->k_med = stof(inputs[15]);
	this->rho_not = stod(inputs[16]);

	this->LorentzDist = inputs[17];
	this->ShellDistParamsFile = inputs[18];
	this->numshells = tw/dte;

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
	
	out_file << "E_dot_iso = " << this->E_dot_iso << "\n";
	out_file << "theta = " << this->theta << "\n";
	out_file << "r_open = " << this->r_open << "\n";
	out_file << "eps_th = " << this->eps_th << "\n";
	out_file << "sigma = " << this->sigma << "\n";

	out_file << "eps_e_int = " << this->eps_e_int << "\n";
	out_file << "eps_b_int = " << this->eps_b_int << "\n";
	out_file << "zeta_int = " << this->zeta_int << "\n";
	out_file << "p_int = " << this->p_int << "\n";

	out_file << "eps_e_ext = " << this->eps_e_ext << "\n";
	out_file << "eps_b_ext = " << this->eps_b_ext << "\n";
	out_file << "zeta_ext = " << this->zeta_ext << "\n";
	out_file << "p_ext = " << this->p_ext << "\n";

	out_file << "k_med = " << this->k_med << "\n";
	out_file << "rho_not = " << this->rho_not << "\n";

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
	std::cout << "E_dot_iso = " << this->E_dot_iso << "\n";
	std::cout << "theta = " << this->theta << "\n";
	std::cout << "r_open = " << this->r_open << "\n";
	std::cout << "eps_th = " << this->eps_th << "\n";
	std::cout << "sigma = " << this->sigma << "\n";
	std::cout << "eps_e_int = " << this->eps_e_int << "\n";
	std::cout << "eps_b_int = " << this->eps_b_int << "\n";
	std::cout << "zeta_int = " << this->zeta_int << "\n";
	std::cout << "p_int = " << this->p_int << "\n";
	std::cout << "eps_e_ext = " << this->eps_e_ext << "\n";
	std::cout << "eps_b_ext = " << this->eps_b_ext << "\n";
	std::cout << "zeta_ext = " << this->zeta_ext << "\n";
	std::cout << "p_ext = " << this->p_ext << "\n";
	std::cout << "k_med = " << this->k_med << "\n";
	std::cout << "rho_not = " << this->rho_not << "\n";
	std::cout << "LorentzDist = " << this->LorentzDist << "\n";
	std::cout << "ShellDistParamsFile = " << this->ShellDistParamsFile << "\n";

}
