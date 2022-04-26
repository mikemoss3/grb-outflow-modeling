/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-24

Class that stores the shells which comprise the prompt jet. 
The shells are distributed according to the selected distribution.
The shells will be used to simulated the prompt emission.
*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "ShellDist.hpp"

using namespace std;

// SynthGRB constructor
ShellDist::ShellDist(int numshells, double E_dot)
{
	this->numshells = numshells; // the total number of shells launched
	this->E_dot = E_dot; // total energy injection rate of the central engine 
	
	// Set the size of the shell vectors
	shell_radius.resize(numshells); // Records shell radius 
	shell_gamma.resize(numshells); // Records shell Lorentz factor
	shell_mass.resize(numshells); // Records shell mass 
	shell_te.resize(numshells); // Records launch time of the shells
	shell_status.resize(numshells); // Status indicator: 0 = deactived, 1 = active and launched, 2 = not launched

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Distribute the shells according to a step function with initial Lorentz factor G1 and secondary Lorentz factor G2
// The step occurs when the cumulative mass of the shells reaches a fraction of the total mass specified by mfrac
void ShellDist::step(float dte, float g1, float g2, float mfrac, bool fluctuations)
{
	/*
	Distribute the Lorentz factors of the shells into a step function. 
	Params: 
	dte = time between shell launches, this can be specific by a single float to apply a constant time step through out the jet evolution or can be a array of the shell emission times
	g1 = Lorentz factor of the group of shells launched earlier
	g2 = Lorentz factor of the group of shells launched later
	mfrac = the mass fraction of the group of shells launched earlier, e.g., mfrac = M_1 / M_total
	*/

	// Number of shells with Lorentz factor g1: 
	int n1 = numshells / ( ((1.-mfrac)*g2/(mfrac*g1)) + 1. );
	// Number of shells with Lorentz factor g2: 
	// int n2 = numshells - n1;

	// Set the Lorentz factors for each section of the step distribution
	double gamma_bar = 0.0;// Average Lorentz factor in the outflow
	for(int i = 0; i<n1; ++i)
	{
		shell_gamma.at(i) = g1;
		gamma_bar += shell_gamma.at(i);
	}
	for(int i = n1; i<numshells; ++i)
	{
		shell_gamma.at(i) = g2;
		gamma_bar += shell_gamma.at(i);
	}
	gamma_bar /= numshells;
	
	for(float i=0; i<numshells; ++i)
	{
		// Set the Mass for each shell 
		// Define the mass as M/M_ave, where M_ave is the average mass per shell (M_ave = M_dot * dt = E_dot *dte /gamma_ave/c^2)
		shell_mass.at(i) = gamma_bar / shell_gamma.at(i);
		
		// Calculate the launch time of each shell since the start of the launch 
		shell_te.at(i) = i*dte;

		// Calculate the initial shell position based on when the shell will be launched
		// Notice this is actually R/c 
		shell_radius.at(i) = - beta(shell_gamma.at(i)) * shell_te.at(i);

		// Deactivate all shells
		shell_status.at(i) = 1;
	}

	// Eliminate possible divide by zero error (still insignificantly small).
	shell_radius.at(0) = 1./c_cm;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Distribute the shells according to a step function with initial Lorentz factor G1 and secondary Lorentz factor G2
// The smooth step occurs when the cumulative mass of the shells reaches a fraction of the total time specified by tfrac
void ShellDist::smoothstep(float dte, float g1, float g2, float tfrac, bool fluctuations)
{
	/*
	Distribute the Lorentz factors of the shells into a step function. 
	Params: 
	dte = time between shell launches, this can be specific by a single float to apply a constant time step through out the jet evolution or can be a array of the shell emission times
	g1 = Lorentz factor of the group of shells launched earlier
	g2 = Lorentz factor of the group of shells launched later
	tfrac = describes when the smooth decay begins 
	*/

	float tw = dte*numshells; // Duration of the wind

	// Set the Lorentz factors for each section of the step distribution
	double gamma_bar = 0.0;// Average Lorentz factor in the outflow
	for(double i = 0; i<numshells; ++i)
	{
		if( dte*i >= tfrac*tw)
		{
			shell_gamma.at(i) = g2;	
		}
		else
		{
			shell_gamma.at(i) = ((g2+g1)/2) - (((g2-g1)/2)*cos(M_PI*(dte*i)/(tfrac*tw)));
		}
		gamma_bar += shell_gamma.at(i);
	}
	gamma_bar /= numshells;


	for(float i=0; i<numshells; ++i)
	{
		// Set the Mass for each shell 
		// Define the mass as M/M_ave, where M_ave is the average mass per shell (M_ave = M_dot * dt = E_dot *dte /gamma_ave/c^2)
		shell_mass.at(i) = gamma_bar / shell_gamma.at(i);
		
		// Calculate the launch time of each shell since the start of the launch 
		shell_te.at(i) = i*dte;

		// Calculate the initial shell position based on when the shell will be launched
		// Notice this is actually R/c 
		shell_radius.at(i) = - beta(shell_gamma.at(i)) * shell_te.at(i);

		// Deactivate all shells
		shell_status.at(i) = 1;
	}

	// Eliminate possible divide by zero error (still insignificantly small).
	shell_radius.at(0) = 1./c_cm;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Distribute the shell according to an oscillatory function described in Hascoet et al. 2013. 
void ShellDist::oscillatory(float dte, float median, float amp, float freq, float decay, bool fluctuations)
{
	/*
	Distribute the Lorentz factors of the shells into a step function. 
	Params: 
	dte = time between shell launches, this can be specific by a single float to apply a constant time step through out the jet evolution or can be a array of the shell emission times
	median = median Lorentz factor of the undamped oscillating distribution 
	amp = amplitude of the undamped oscillating distribution 
	freq = frequency of the oscillations in the distribution
	decay = used to modulate the decay rate, number multiplied by the half life time scale
	*/

	// Set the Lorentz factors for each section of the step distribution
	double gamma_bar = 0.0;// Average Lorentz factor in the outflow
	for(double i = 0; i<numshells; ++i)
	{
		shell_gamma.at(i) = median * ( 1. + amp * cos( freq * M_PI * (1. - i/numshells) ) )*exp(-decay*i/numshells);
		gamma_bar += shell_gamma.at(i);
	}
	gamma_bar /= numshells;
	
	for(float i=0; i<numshells; ++i)
	{
		// Set the Mass for each shell 
		// Define the mass as M/M_ave, where M_ave is the average mass per shell (M_ave = M_dot * dt = E_dot *dte /gamma_ave/c^2)
		shell_mass.at(i) = gamma_bar / shell_gamma.at(i);
		
		// Calculate the launch time of each shell since the start of the launch 
		shell_te.at(i) = i*dte;

		// Calculate the initial shell position based on when the shell will be launched
		// Notice this is actually R/c 
		shell_radius.at(i) = - beta(shell_gamma.at(i)) * shell_te.at(i);

		// Activate all shells
		shell_status.at(i) = 1;
	}

	// Eliminate possible divide by zero error (still insignificantly small).
	shell_radius.at(0) = 1./c_cm;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Distribute the shells by injecting Gaussians at specific places. 
void ShellDist::gauss_inject(float dte, float gamma_ave, float decay, int num_gauss, std::vector<float> means, std::vector<float> amps, std::vector<float> sigmas, std::vector<float> sigma_count, bool fluctuations)
{
	/*
	Distribute the Lorentz factors of the shells into a step function. 
	Params: 
	dte = time between shell launches, this can be specific by a single float to apply a constant time step through out the jet evolution or can be a array of the shell emission times

	means = 
	sigmas = 
	amps = 

	decay = used to modulate the decay rate, number multiplied by the half life time scale
	*/

	// First, check if each input array is the same size
	if( (means.size() != (size_t)num_gauss) | (sigmas.size() != (size_t)num_gauss) | (amps.size() != (size_t)num_gauss) )
	{
		cout << "Input Gaussian parameter arrays must have a length equal to the specified Number of Gaussians.";
		return;
	}

	// Initially, we set all Lorentz factors to the input average value
	for(double i = 0; i<numshells; ++i)
	{
		shell_gamma.at(i) = gamma_ave * exp(-decay*i/numshells);
	}

	// Set the Lorentz factors for each section of the step distribution

	// Now, add a Gaussian at each value in the "means" array	
	float curr_t_start = 0.; // Used to find the start time of each Gaussian
	float curr_t_stop = 0.; // Used to find the stop time of each Gaussian
	for(int k=0; k < num_gauss; ++k)
	{
		curr_t_start = means.at(k) - (sigma_count.at(k) * sigmas.at(k) );
		if(curr_t_start < 0){curr_t_start = 0;} // Set the start time limit to zero if it is was found to be less than zero

		curr_t_stop = means.at(k) + (sigma_count.at(k) * sigmas.at(k) );
		if(curr_t_stop > (dte*numshells)){curr_t_stop = (dte*numshells);} // Set the stop time limit to the end duration of the wind if it is was found to be larger than the wind duration

		for(double i = (curr_t_start/dte); i < (curr_t_stop/dte); ++i)
		{
			shell_gamma.at(i) += amps.at(k) * exp( - pow( (i*dte) - means.at(k),2.) / (2.*pow(sigmas.at(k),2.)) );
		}
	}

	// For all shells that still have a Lorentz factor of 0, deactivate them,
	for(float i=0; i<numshells; ++i)
	{
		if(shell_gamma.at(i) <= 1.)
		{
			// De-activate all shells
			shell_status.at(i) = 0;
		}
		else
		{
			// Activate all shells
			shell_status.at(i) = 1;
		}
	}

	// Find the average Lorentz factor of the active shells in the outflow
	double gamma_bar = 0.0;
	for(double i = 0; i<numshells; ++i)
	{
		if(shell_status.at(i) == 1)
		{
			gamma_bar += shell_gamma.at(i);
		}
	}
	gamma_bar /= numshells;


	// For the active shells, find the mass and radius
	for(float i=0; i<numshells; ++i)
	{
		if(shell_status.at(i) == 1)
		{			
			// Set the Mass for each shell 
			// Define the mass as M/M_ave, where M_ave is the average mass per shell (M_ave = M_dot * dt = E_dot *dte /gamma_ave/c^2)
			shell_mass.at(i) = gamma_bar / shell_gamma.at(i);
			
			// Calculate the launch time of each shell since the start of the launch 
			shell_te.at(i) = i*dte;

			// Calculate the initial shell position based on when the shell will be launched
			// Notice this is actually R/c 
			shell_radius.at(i) = - beta(shell_gamma.at(i)) * shell_te.at(i);
		}
	}

	// Eliminate possible divide by zero error (still insignificantly small).
	shell_radius.at(0) = 1./c_cm;

	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Distribute the shells by injecting uniform distributions at specific places. 
void ShellDist::square_inject(float dte, float gamma_ave, float decay, int num_squares, std::vector<float> starts, std::vector<float> durations, std::vector<float> amps, bool fluctuations)
{
	/*
	Distribute the Lorentz factors of the shells into a step function. 
	Params: 
	dte = time between shell launches, this can be specific by a single float to apply a constant time step through out the jet evolution or can be a array of the shell emission times

	means = 
	sigmas = 
	amps = 

	decay = used to modulate the decay rate, number multiplied by the half life time scale
	*/


	// Initially, we set all Lorentz factors to the input average value
	for(double i = 0; i<numshells; ++i)
	{
		shell_gamma.at(i) = gamma_ave * exp(-decay*i/numshells);
	}

	// Set the Lorentz factors for each section of the step distribution

	// Now, add a square starting at each value in the "starts" array	
	float curr_t_start = 0.; // Used to find the start time of each Gaussian
	float curr_t_stop = 0.; // Used to find the stop time of each Gaussian
	for(int k=0; k < num_squares; ++k)
	{
		curr_t_start = starts.at(k);
		if(curr_t_start < 0){curr_t_start = 0;} // Set the start time limit to zero if it is was found to be less than zero

		curr_t_stop = starts.at(k) +  durations.at(k) ;
		if(curr_t_stop > (dte*numshells)){curr_t_stop = (dte*numshells);} // Set the stop time limit to the end duration of the wind if it is was found to be larger than the wind duration

		for(double i = (curr_t_start/dte); i < (curr_t_stop/dte); ++i)
		{
			shell_gamma.at(i) +=  amps[k];
		}
	}

	// For all shells that still have a Lorentz factor of 0, deactivate them,
	for(float i=0; i<numshells; ++i)
	{
		if(shell_gamma.at(i) <= 1.)
		{
			// De-activate all shells
			shell_status.at(i) = 0;
		}
		else
		{
			// Activate all shells
			shell_status.at(i) = 1;
		}
	}

	// Find the average Lorentz factor in the outflow
	double gamma_bar = 0.0;
	for(double i = 0; i<numshells; ++i)
	{
		if(shell_status.at(i) == 1)
		{
			gamma_bar += shell_gamma.at(i);
		}
	}
	gamma_bar /= numshells;


	for(float i=0; i<numshells; ++i)
	{

		if(shell_status.at(i) == 1)
		{
			// Set the Mass for each shell 
			// Define the mass as M/M_ave, where M_ave is the average mass per shell (M_ave = M_dot * dt = E_dot *dte /gamma_ave/c^2)
			shell_mass.at(i) = gamma_bar / shell_gamma.at(i);
			
			// Calculate the launch time of each shell since the start of the launch 
			shell_te.at(i) = i*dte;

			// Calculate the initial shell position based on when the shell will be launched
			// Notice this is actually R/c 
			shell_radius.at(i) = - beta(shell_gamma.at(i)) * shell_te.at(i);
		}
	}

	// Eliminate possible divide by zero error (still insignificantly small).
	shell_radius.at(0) = 1./c_cm;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make a linear distribution
void ShellDist::linear(float dte, float g1, float g2, bool fluctuations)
{
	/*
	Distribute the Lorentz factors of the shells as a linear function
	Params: 
	dte = time between shell launches, this can be specific by a single float to apply a constant time step through out the jet evolution or can be a array of the shell emission times
	g1 = Lorentz factor at the front of the ejecta
	g2 = Lorentz factor at the end of the ejecta

	*/

	// Calculate slope
	float xstart = 0.; // Start time of the ejecta
	float xend = this->numshells * dte; // end time of the ejecta
	float m = (g2-g1)/(xend - xstart);

	// Assign Lorentz factors to each shell
	for(double i = 0; i<numshells; ++i)
	{
		shell_gamma.at(i) = m*i*dte + g1;
	}

	// Find the average Lorentz factor in the outflow
	double gamma_bar = 0.0;
	for(double i = 0; i<numshells; ++i)
	{
		gamma_bar += shell_gamma.at(i);
	}
	gamma_bar /= numshells;


	for(float i=0; i<numshells; ++i)
	{

		// Set the Mass for each shell 
		// Define the mass as M/M_ave, where M_ave is the average mass per shell (M_ave = M_dot * dt = E_dot *dte /gamma_ave/c^2)
		shell_mass.at(i) = gamma_bar / shell_gamma.at(i);
		
		// Calculate the launch time of each shell since the start of the launch 
		shell_te.at(i) = i*dte;

		// Calculate the initial shell position based on when the shell will be launched
		// Notice this is actually R/c 
		shell_radius.at(i) = - beta(shell_gamma.at(i)) * shell_te.at(i);

		// Activate all shells
		shell_status.at(i) = 1;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write shell distribution to a text file
void ShellDist::WriteToTXT(string filename, double time, bool append)
{

	// write out thermal params
	ofstream shell_dist_file; // Construct file 
	if(append == false)
	{
		shell_dist_file.open(filename); // Open text file with this name
	}
	else
	{
		shell_dist_file.open(filename, std::ios_base::app); // Append to text file with this name 
	}
	size_t i=0;
	shell_dist_file << "// Next step" << "\n";
	shell_dist_file << time << "\n";
	// For each shell, write the shell parameters to file.
	while ( i < shell_radius.size())
	{
		shell_dist_file << shell_radius.at(i);
		shell_dist_file << " ";
		shell_dist_file << shell_gamma.at(i);
		shell_dist_file << " ";
		shell_dist_file << shell_mass.at(i);
		shell_dist_file << " ";
		shell_dist_file << shell_te.at(i);
		shell_dist_file << " ";
		shell_dist_file << shell_status.at(i);
		shell_dist_file << "\n";
		++i;
	}
	shell_dist_file.close(); // Close file
}

