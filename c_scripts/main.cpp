/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Main function to create and use Spectra, Light Curves, and Response functions

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "Response.hpp"
#include "LightCurve.hpp"
#include "Spectrum.hpp"
#include "cosmology.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
	Response instrument_response = Response();
	std::string resp_file_name = "/Users/mjmoss/Research/instr-rsp-investigate/example-instr-files/sw00883832000b_1speak.rsp";
	instrument_response.LoadRespFromFile(resp_file_name);

	float tmin = 0;
	float tmax = 20;
	float dt = 0.1;
	float z = 1;

	float Emin = 8;
	float Emax = 40000;
	float num_E_bins = 100;

	// float Emin = instrument_response.phot_energ_min;
	// float Emax = instrument_response.phot_energ_max;
	// float num_E_bins = instrument_response.num_phot_bins;	

	Spectrum source_spectrum_tot = Spectrum(Emin, Emax, num_E_bins, tmin, tmax, z);
	Spectrum source_spectrum_therm = Spectrum(Emin, Emax, num_E_bins, tmin, tmax, z);
	Spectrum source_spectrum_synch = Spectrum(Emin, Emax, num_E_bins, tmin, tmax, z);
	
	std::string therm_emission_file_name = "../sim_results/ordlor_spectrum_therm.txt";
	std::string synch_emission_file_name = "../sim_results/ordlor_spectrum_synch.txt";
	
	// source_spectrum_tot.AddThermComp(therm_emission_file_name, tmin, tmax);
	// source_spectrum_tot.AddSynchComp(synch_emission_file_name, tmin, tmax);
	// source_spectrum_therm.AddThermComp(therm_emission_file_name, tmin, tmax);
	// source_spectrum_synch.AddSynchComp(synch_emission_file_name, tmin, tmax);

	// source_spectrum_tot.WriteToTXT("source_spectrum_tot.txt");
	// source_spectrum_therm.WriteToTXT("source_spectrum_therm.txt");
	// source_spectrum_synch.WriteToTXT("source_spectrum_synch.txt");
	
	LightCurve lightcurve_tot = LightCurve(Emin, Emax, num_E_bins, tmin , tmax, dt, z);
	LightCurve lightcurve_therm = LightCurve(Emin, Emax, num_E_bins, tmin , tmax, dt, z);
	LightCurve lightcurve_synch = LightCurve(Emin, Emax, num_E_bins, tmin , tmax, dt, z);
	
	// lightcurve_tot.set_instrument_response(&instrument_response);
	// lightcurve_therm.set_instrument_response(&instrument_response);
	// lightcurve_synch.set_instrument_response(&instrument_response);
	
	lightcurve_tot.AddThermalLightCurve(therm_emission_file_name);
	lightcurve_tot.AddSynchLightCurve(synch_emission_file_name);
	lightcurve_therm.AddThermalLightCurve(therm_emission_file_name);
	lightcurve_synch.AddSynchLightCurve(synch_emission_file_name);
	
	lightcurve_tot.WriteToTXT("test_tot.txt");
	lightcurve_therm.WriteToTXT("test_therm.txt");
	lightcurve_synch.WriteToTXT("test_synch.txt");
	
	return 0;
}
