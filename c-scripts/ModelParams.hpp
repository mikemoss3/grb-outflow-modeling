/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-31

Header file for ModelParameter.cpp class

*/

// Ensure that C++ doesn't have any problems if this header is accessed by multiple scripts simultaneously
#ifndef MODELPARAMSCLASS_H
#define MODELPARAMSCLASS_H

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries

using namespace std;

// Declare light curve class
class ModelParams
{
public:
	// ModelParams constructor
	ModelParams(float tw = 10.0, float dte = 0.002, float alpha_e = 1./3, float alpha_b = 1./3, float zeta_int = 1e-3, float zeta_ext = 1., double E_dot_iso = 1e53, float theta = 0.1, float r_open = 1e6, float eps_th = 0.03, float sigma = 0.1, float p = 2.5, std::string LorentzDist = "step", std::string ShellDistParamsFile = "Default");

	// ModelParams member variables
	// jet parameters:
	float tw; // sec, Duration of the wind
	float dte; // sec, Time between successive shell launches
	float eps_e; //  Fraction of dissipated energy that goes into the electrons 
	float eps_b; // Fraction of dissipated energy that goes into the magnetic field 
	float zeta_int; // Fraction of electrons which are accelerated in internal and reverse 
	float zeta_ext; // Fraction of electrons which are accelerated in forward shock
	double E_dot_iso; // erg/s, Injected isotropic energy rate
	float theta; // radians, Half-opening angle of the jet
	float r_open; // cm, Opening radius of the jet
	float eps_th; // Fraction of energy in the outflow in the form of thermal energy 
	float sigma; // Magnetization of the outflow 
	float p; // Power law index of the electron population 
	std::string LorentzDist; // Distribution of the jet shells
	std::string ShellDistParamsFile; // File that contains the parameters to create the distribution of jet shells

	// Calculated as tw/dte
	int numshells; // Number of shells in the jet


	// These vectors can be used to define the parameters space of the variable
	std::vector<float> eps_e_vec;
	std::vector<float> eps_b_vec;
	std::vector<float> zeta_int_vec;
	std::vector<float> zeta_ext_vec;
	std::vector<double> E_dot_iso_vec;
	std::vector<float> theta_vec;
	std::vector<float> r_open_vec;
	std::vector<float> eps_th_vec;
	std::vector<float> sigma_vec;
	std::vector<float> p_vec;

	// Set parameter space
	void set_param_space(
		std::vector<float> eps_e_vec, 
		std::vector<float> eps_b_vec, 
		std::vector<float> zeta_int_vec, 
		std::vector<float> zeta_ext_vec, 
		std::vector<double> E_dot_iso_vec, 
		std::vector<float> theta_vec, 
		std::vector<float> r_open_vec, 
		std::vector<float> eps_th_vec, 
		std::vector<float> sigma_vec,
		std::vector<float> p_vec);
	
	// Copy a ModelParams object
	void copy(ModelParams * p_model_params_in);

	// Load ModelParams from a text file
	void LoadFromTXT(std::string filename);
	
	// Write ModelParams to a text file
	void WriteToTXT(std::string filename);

	// ModelParams member methods
	void PrintAllParams();

// private:

};

#endif