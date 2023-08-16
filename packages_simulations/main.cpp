/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

Main function to create and use Spectra, Light Curves, and Response functions

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
#include "SynthGRB.hpp"
#include "ObsGRB.hpp"
#include "Response.hpp"
// #include "LightCurve.hpp"
// #include "Spectrum.hpp"
#include "cosmology.hpp"
#include "utilfuncs.hpp"
#include "DataAnalysis.hpp"
#include "TTEs.hpp"
#include "ModelParams.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
	// Read in command line variables
	std::string save_dir = argv[1]; // indicates which directory to save everything (default is ../files-data/synthetic-data)
	std::string param_file_name = argv[2]; // indicates the jet parameter file (default is jet-params.txt)
	bool write_lor_prof; // indicates whether to or not write the Lorentz distribution to a file
	isstringstream(argv[3]) >> write_lor_prof;
	bool track_lor_prof_evo; // indicates whether or not to write snapshots of the Lorentz distribution at different times
	isstringstream(argv[4]) >> track_lor_prof_evo;
	bool write_jet_params; // indicates whether or not to write out the jet parameters calculated during the simulation to a file
	isstringstream(argv[5]) >> write_jet_params;
	bool light_curve; // indicates whether to calculate a light curve
	isstringstream(argv[6]) >> light_curve;
	float light_curve_time_start = atof(argv[8]); // time interval start to calculate the light curve over
	float light_curve_time_stop = atof(argv[9]); // time interval stop to calculate the light curve over
	float light_curve_energy_min = atof(argv[10]); // energy interval min to calculate the light curve over
	float light_curve_energy_max = atof(argv[11]); // energy interval max to calculate the light curve over
	bool time_log_scale; // indicates whether the time axis should be generated in linear space (False) or log space (True)
	isstringstream(argv[12]) >> time_log_scale;
	float time_step = atof(argv[13]); // the time step along the time axis
	bool spectra; // indicates whether to calculate spectra
	isstringstream(argv[14]) >> spectra;
	float spectra_energy_min = atof(argv[15]); // energy interval min to calculate the spectra over
	float spectra_energy_max = atof(argv[16]); // energy interval max to calculate the spectra over
	int num_energy_bins = atoi(argv[17]); // number of energy bins along the energy axis
	int num_spectra = atoi(argv[18]); // number of spectra to be created
	bool inc_comps; // indicates whether to record the separate emission components that compose the generated spectra or light curves
	isstringstream(argv[19]) >> inc_comps;

	float spectra_time_start = 0.; // time interval start to calculate the spectra over
	float spectra_time_stop = 0.; // time interval stop to calculate the spectra over


	// Create SynthGRB instance
	SynthGRB synth_grb = SynthGRB();

	// Load jet parameters from the supplied file name 
	synth_grb.LoadJetParamsFromTXT(param_file_name);
	
	// Write Lorentz distribution to file if desired 
	if(write_lor_prof){(*synth_grb.p_jet_shells).WriteToTXT(save_dir+"synthGRB_shell_dist.txt");}
	synth_grb.anim_lor_dist = track_lor_prof_evo; // Write evolution of Lorentz distribution to file

	// Simulate jet Dynamics
	synth_grb.SimulateJetDynamics();

	// Write jet parameters to file if desired
	if(write_jet_params){synth_grb.write_out_jet_params(save_dir);}

	// Make total light curve
	synth_grb.make_source_light_curve(light_curve_energy_min, light_curve_energy_max, light_curve_time_start, light_curve_time_stop, time_step, time_log_scale);
	synth_grb.WriteLightCurveToTXT(save_dir+"synthGRB_light_curve.txt");

	// Make total spectra (for each supplied time interval)
	for(int i = 0; i < num_spectra; ++i)
	{
		spectra_time_start = atof(argv[20+(2*i)]);
		spectra_time_stop = atof(argv[21+(2*i)]);

		// Total spectrum
		synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi);
		synth_grb.WriteSpectrumToTXT(save_dir+"synthGRB_spec_TOT.txt");
	}

	// Make component light curves and spectra if desired
	if(inc_comps)
	{
		const char* comp[4]= { "TH", "IS", "FS", "RS" };
		
		// Component light curves
		for(int i = 0; i < 4; ++i)
		{
			synth_grb.make_source_light_curve(light_curve_energy_min, light_curve_energy_max, light_curve_time_start, light_curve_time_stop, time_step, comp[i], time_log_scale);
			synth_grb.WriteLightCurveToTXT(save_dir+"synthGRB_light_curve_"+comp[i]+".txt");
		}
	
		// Component spectra
		for(int i = 0; i < num_spectra; ++i)
		{
			spectra_time_start = atof(argv[20+(2*i)]);
			spectra_time_stop = atof(argv[20+(2*i)]);

			for(int j = 0; j < 4; ++j)
			{
				synth_grb.make_source_spectrum(spectra_energy_min, spectra_energy_max, num_energy_bins, spectra_time_start, spectra_time_stop, comp[j]);
				synth_grb.WriteSpectrumToTXT(save_dir+"synthGRB_spec_"+comp[j]+".txt");
			}
		}
	}

	return 0;
}


// // Directory path to data storage directory
// std::string dir="../files-data/synthetic-data/";

// /* Testing SynthGRB default Light Curve and Spectrum making */	
// float energ_min = 8;
// float energ_max = 5.e4;
// float energ_min_lc = energ_min; 
// float energ_max_lc = energ_max;
// float num_energ_bins = 400;

// // 15 - 5000 keV = Konus-Wind energy band

// // Light curve time interval
// float tmin = 0.;
// float tmax = 10.;
// float dt = 0.05;

// SynthGRB test_grb = SynthGRB();
// test_grb.LoadJetParamsFromTXT("files-input/jet-params.txt");

// test_grb.anim_lor_dist = true;
// (*test_grb.p_jet_shells).WriteToTXT(dir+"synthGRB_shell_dist.txt");

// test_grb.SimulateJetDynamics();
// test_grb.write_out_jet_params(dir);

// // // Spectrum time interval
// float tlo = 0.5; // The T90? essentially.
// float thi = 1.;

// // // Total spectrum
// test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi);
// test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_TOT.txt");

// // // Component spectrum
// test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "TH");
// test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_TH.txt");
// test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "IS");
// test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_IS.txt");
// test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "FS");
// test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_FS.txt");
// test_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tlo, thi, "RS");
// test_grb.WriteSpectrumToTXT(dir+"synthGRB_spec_RS.txt");

// // Total light curve
// test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt);
// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve.txt");

// // Component light curves
// test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "TH", false);
// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_TH.txt");
// test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "IS", false);
// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_IS.txt");
// test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "FS", false);
// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_FS.txt");
// test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, tmin, tmax, dt, "RS", false);
// test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_RS.txt");

// // Looking at after glows 
// // test_grb.make_source_light_curve(energ_min_lc, energ_max_lc, 15, 4e5, dt, "FS", true);
// // test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_afterglow_gbm.txt");
// // test_grb.make_source_light_curve(0.2, 10, 1e2, 1e6, dt, "FS", true);
// // test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_afterglow_xrt.txt");
// // test_grb.make_source_light_curve(1e-3, 5e-3, 1e2, 1e6, dt, "FS", true);
// // test_grb.WriteLightCurveToTXT(dir+"synthGRB_light_curve_afterglow_opt.txt");

// return 0;
