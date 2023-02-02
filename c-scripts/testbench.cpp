#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <tuple>
#include "fitsio.h"
#include <list>

using namespace std;

// Import Custom Libraries
#include "SynthGRB.hpp"
// #include "SynthGRBLibrary.hpp"
// #include "ObsGRB.hpp"
// #include "Response.hpp"
#include "LightCurve.hpp"
#include "Spectrum.hpp"
#include "cosmology.hpp"
#include "utilfuncs.hpp"
// #include "DataAnalysis.hpp"
// #include "TTEs.hpp"
#include "ModelParams.hpp"
#include "ShellDist.hpp"
#include "GammaFunction.hpp"

using namespace std;


int main(int argc, char const *argv[])
{

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* Testing SynthGRB spectrum making */

	/*
	// Initialize all jet parameters to be given to SynthGRB
	// float tw =  50;
	// float dte =  1;
	// // int numshells =  50;
	// float eps_e =  1./3.;
	// float eps_b =  1./3.;
	// float zeta_int =  1e-3;
	// float zeta_ext =  1.;
	// double E_dot_iso = 1e53; 
	// float theta =  0.1;
	// float r_open = 1e6; 
	// float eps_th = 0.03;
	// float sigma = 0.1;
	// float p = 2.5;
	// std::string LorentzDist = "step";
	// std::string ShellDistParamsFile = "Default";
	// std::string ShellDistParamsFile = "./input-files/jet-shells-step.txt";

	// SynthGRB * test_grb_time = new SynthGRB(tw, dte, eps_e, eps_b, zeta_int, zeta_ext, E_dot_iso, theta, r_open, eps_th, sigma, p, LorentzDist, ShellDistParamsFile);

	// Or use a text file with specified jet parameters
	std::string JetParamFile = "./input-files/jet-params.txt";

	SynthGRB * test_grb_time = new SynthGRB();
	(*test_grb_time).LoadJetParamsFromTXT(JetParamFile);

	// Simulate Jet Dynamics
	(*test_grb_time).SimulateJetDynamics();
	
	// Write out jet parameters
	(*test_grb_time).write_out_jet_params("./data-file-dir/");

	// Initialize variables for spectrum creation
	float energ_min = 0.1;
	float energ_max = 1e5;
	int num_energ_bins = 600;
	float tmin = 0.;
	float tmax = 30.;

	// Make spectrum
	(*test_grb_time).make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax);
	
	(*test_grb_time).WriteSpectrumToTXT("data-file-dir/test_spec.txt");

	(*test_grb_time).make_source_light_curve(50, 300, tmin, tmax, 0.05);
	(*test_grb_time).WriteLightCurveToTXT("data-file-dir/test_light_curve.txt");
	*/



	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* Testing DataAnalyis Fitting Package */
	
	
	/*
	// Define parameters for the synthetic source 
	// float tw = 50; 
	// float dte = 0.1;
	// // int numshells = 500;
	// float eps_e = 0.33;
	// float eps_b = 0.33;
	// float zeta_int = 1e-3;
	// float zeta_ext = 1.;
	// double E_dot_iso = 1e53;
	// float theta = 0.1;
	// float r_open = 1e6;
	// float eps_th = 0.03;
	// float sigma = 0.1;
	// float p = 2.5; 
	// std::string LorentzDist = "step";
	// std::string ShellDistParamsFile = "Default";

	// Define parameter space to explore
	std::vector<float> eps_e_vec = {0.33, 0.5};
	std::vector<float> eps_b_vec = {0.33, 0.5};
	std::vector<float> zeta_int_vec = {1e-3, 1.};
	std::vector<float> zeta_ext_vec = {1.};
	std::vector<double> E_dot_iso_vec = {1e52, 1e53};
	std::vector<float> theta_vec = {0.1, 0.2};
	std::vector<float> r_open_vec = {1e6};
	std::vector<float> eps_th_vec = {0.03};
	std::vector<float> sigma_vec = {0.1};
	std::vector<float> p_vec = {2.2,2.5};



	// Make instrument response object
	Response instrument_response = Response();
	// Load instrument response 
	// instrument_response.LoadRespFromFile("~/Research/instr-rsp-investigate/example-instr-files/sw00883832000b_1speak.rsp");
	instrument_response.LoadRespFromFile("~/Research/instr-rsp-investigate/example-instr-files/glg_cspec_n6_bn160625945_v04.rsp");
	// instrument_response.LoadRespFromFile("~/Research/instr-rsp-investigate/example-instr-files/glg_cspec_b0_bn160625945_v04.rsp");

	
	// Make "observed" GRB data from a synthetic GRB and the given instrument response matrix 
	// Initialize synthetic GRB pointer
	// SynthGRB * synth_obs_grb = new SynthGRB(tw, dte, eps_e, eps_b, zeta_int, zeta_ext, E_dot_iso, theta, r_open, eps_th, sigma, p, LorentzDist, ShellDistParamsFile);
	SynthGRB * synth_obs_grb = new SynthGRB();
	(*synth_obs_grb).LoadJetParamsFromTXT("input-files/jet-params.txt");

	// Simulate the Jet Dynamics
	(*synth_obs_grb).SimulateJetDynamics();
	// Make spectrum from jet simulation
	(*synth_obs_grb).make_source_spectrum(instrument_response.phot_energ_min, instrument_response.phot_energ_max,instrument_response.num_phot_bins);
	// Write out synthetic spectrum
	(*(*synth_obs_grb).p_source_spectrum).WriteToTXT("data-file-dir/spec_source.txt");
	// Make an a more realistic spectrum by adding fluctuations
	(*(*synth_obs_grb).p_source_spectrum).add_fluctuations();
	// Add uncertainties, although I don't think it makes sense to add uncertainties to a spectrum before convolving with an instrument response matrix
	(*(*synth_obs_grb).p_source_spectrum).add_unc(0.1);
	// Write out fluctuated spectrum 
	(*(*synth_obs_grb).p_source_spectrum).WriteToTXT("data-file-dir/spec_source_fluc.txt");


	// Initialize a spectrum pointer that will point to a convolved spectrum (i.e., the synthetic spectrum folded with the instrument response matrix) 
	// Also knows as the synthetic observed spectrum
	Spectrum * p_folded_spectrum = new Spectrum(instrument_response.chan_energ_min, instrument_response.chan_energ_max,instrument_response.num_chans);
	// Convolve synthetic spectrum and instrument response matrix
	ConvolveSpectra(p_folded_spectrum, (*(*synth_obs_grb).p_source_spectrum), instrument_response, true);
	(*p_folded_spectrum).add_unc(0.1);
	//// We now have our synthetic observed spectrum 
	// Write convolved spectrum to text file
	(*p_folded_spectrum).WriteToTXT("data-file-dir/spec_obs.txt");

	// Make Data Analysis object
	DataAnalysis data_analysis = DataAnalysis();

	// Fit the synthetic observed spectrum and see if the input parameters are recovered by the fitting algorithm.
	// data_analysis.FitSpectrum_SynthGRB();
	// Or use an empirical function
	// data_analysis.add_fit_func("PL");
	data_analysis.add_fit_func("BPL");
	// Fit the synthetic observed spectrum with
	data_analysis.FitSpectrum();

	// Using the best fit parameters, print out the best fit spectrum
	Spectrum * p_mod_unf = new Spectrum(instrument_response.phot_energ_min, instrument_response.phot_energ_max,instrument_response.num_phot_bins);
	Spectrum * p_mod_fold = new Spectrum(instrument_response.chan_energ_min, instrument_response.chan_energ_max,instrument_response.num_chans);

	// Make spectrum using the best fit model parameters
	data_analysis.make_model_spectrum(p_mod_unf, data_analysis.best_fit_params);

	ConvolveSpectra(p_mod_fold, (*p_mod_unf), instrument_response, true);
	(*p_mod_fold).WriteToTXT("data-file-dir/spec_mod_emp.txt");
	*/

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* Testing ObsGRB */

	/*	
	// GBM NaI TTE Data file name
	std::string file_name = "./data-file-dir/glg_tte_n4_bn190114873_v00.fit";

	// Construct a GRB object and assign a pointer to it
	ObsGRB * p_obs_grb = new ObsGRB();

	// Load data from the TTE file, this will also make a spectrum and light curve (full time and full energy band)
	(*p_obs_grb).LoadTTEData(file_name);
	// If you want to make a Light Curve with a certain time interval and energy band.
	(*p_obs_grb).MakeTTELightCurve(-10,75,10,1000);
	// Write light curve to text file
	(*(*p_obs_grb).p_obs_light_curve).WriteToTXT("data-file-dir/190114C_n4_tte_light_curve.txt");


	// If you want to make a Spectrum with a certain time interval and energy band.
	// Background
	(*p_obs_grb).MakeTTESpectrum(-30,-20,10,1000);
	// Write spectrum at different times to text files
	(*(*p_obs_grb).p_obs_spectrum).WriteToTXT("data-file-dir/190114C_n4_tte_spec_bak.txt");
	// Rise
	(*p_obs_grb).MakeTTESpectrum(0,1,10,1000);  
	(*(*p_obs_grb).p_obs_spectrum).WriteToTXT("data-file-dir/190114C_n4_tte_spec_rise.txt");
	// Peak
	(*p_obs_grb).MakeTTESpectrum(2.5,3,10,1000);
	(*(*p_obs_grb).p_obs_spectrum).WriteToTXT("data-file-dir/190114C_n4_tte_spec_peak.txt");
	// Fall
	(*p_obs_grb).MakeTTESpectrum(6,7,10,1000);
	(*(*p_obs_grb).p_obs_spectrum).WriteToTXT("data-file-dir/190114C_n4_tte_spec_fall.txt");

	*/
	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* Testing multiple spectra+response combinations*/

	/*
	// Make instrument response object
	Response instrument_response_1 = Response();
	Response instrument_response_2 = Response();

	// Load instrument response 
	instrument_response_1.LoadRespFromFile("~/Research/instr-rsp-investigate/example-instr-files/glg_cspec_n6_bn160625945_v04.rsp");
	instrument_response_2.LoadRespFromFile("~/Research/instr-rsp-investigate/example-instr-files/glg_cspec_b0_bn160625945_v04.rsp");
	

	// Make "observed" GRB data from a synthetic GRB and the given instrument response matrix 
	
	// Use inputs from text file
	std::string JetParamFile = "./input-files/jet-params.txt";
	// Initialize synthetic GRB pointer
	SynthGRB * synth_obs_grb = new SynthGRB();
	(*synth_obs_grb).LoadJetParamsFromTXT(JetParamFile);

	// Simulate the Jet Dynamics
	(*synth_obs_grb).SimulateJetDynamics();

	// Make spectrum from jet simulation for each instrument response function
	(*synth_obs_grb).make_source_spectrum(instrument_response_1.phot_energ_min, instrument_response_1.phot_energ_max,instrument_response_1.num_phot_bins);
	// Make an a more realistic spectrum by adding fluctuations
	(*(*synth_obs_grb).p_source_spectrum).add_fluctuations();	
	// Make a pointer for each synthetic observed spectrum (one for response matrix)
	Spectrum * p_folded_spectrum_1 = NULL;
	// Convolve the spectrum with each response matrix	
	ConvolveSpectra(p_folded_spectrum_1, (*(*synth_obs_grb).p_source_spectrum), instrument_response_1, true);
	// Add uncertainty column
	(*p_folded_spectrum_1).add_unc(0.1);
	(*p_folded_spectrum_1).WriteToTXT("data-file-dir/spec_obs_1.txt");

	// Repeat for other instrument response 
	(*synth_obs_grb).make_source_spectrum(instrument_response_2.phot_energ_min, instrument_response_2.phot_energ_max,instrument_response_2.num_phot_bins);
	(*(*synth_obs_grb).p_source_spectrum).add_fluctuations();	
	Spectrum * p_folded_spectrum_2 = NULL;
	ConvolveSpectra(p_folded_spectrum_2, (*(*synth_obs_grb).p_source_spectrum), instrument_response_2, true);
	(*p_folded_spectrum_2).add_unc(0.1);	
	(*p_folded_spectrum_2).WriteToTXT("data-file-dir/spec_obs_2.txt");


	// Make Data Analysis object
	DataAnalysis data_analysis = DataAnalysis();

	// Set initial parameters for the fit
	data_analysis.set_init_params(10, 0.002, 0.4, 0.4, 0.1, 1e52, 0.15, 1e6, 0.03, 0.1, "step", "Default");
	// Define parameter space to explore
	std::vector<float> eps_e_vec = {0.33, 0.5};
	std::vector<float> eps_b_vec = {0.33, 0.5};
	std::vector<float> zeta_int_vec = {1e-3, 1.};
	std::vector<float> zeta_ext_vec = {1.};
	std::vector<double> E_dot_iso_vec = {1e52, 1e53};
	std::vector<float> theta_vec = {0.1, 0.2};
	std::vector<float> r_open_vec = {1e6};
	std::vector<float> eps_th_vec = {0.03};
	std::vector<float> sigma_vec = {0.1};
	std::vector<float> p_vec = {2.2, 2.5};
	// Define parameter space for the fitting algorithm to explore
	data_analysis.set_param_space(eps_e_vec, eps_b_vec, zeta_int_vec, zeta_ext_vec, E_dot_iso_vec, theta_vec, r_open_vec, eps_th_vec, sigma_vec, p_vec);

	// Load observations and response into the data_analysis object
	data_analysis.LoadSpecAndResp(p_folded_spectrum_1,&instrument_response_1);
	data_analysis.LoadSpecAndResp(p_folded_spectrum_2,&instrument_response_2);

	// Fit the synthetic observed spectrum and see if the input parameters are recovered by the fitting algorithm.
	data_analysis.FitSpectrum_SynthGRB();

	// // Print best fit parameters
	(*data_analysis.p_curr_model_params).PrintAllParams();


	// Make a new GRB using the best fit parameters
	std::string JetParamFile = "./input-files/jet-params.txt";
	// Initialize synthetic GRB pointer
	SynthGRB * model_grb = new SynthGRB();
	(*model_grb).LoadJetParamsFromTXT(JetParamFile);

	// Simulated the jet dynamics and create the spectrum
	(*model_grb).SimulateJetDynamics();
	
	// For instrument response 1
	(*model_grb).make_source_spectrum(instrument_response_1.chan_energ_min, instrument_response_1.chan_energ_max,instrument_response_1.num_chans);
	(*(*model_grb).p_source_spectrum).WriteToTXT("data-file-dir/spec_model_1.txt");
	// Spectrum * p_folded_model_spectrum_1 = NULL;
	// ConvolveSpectra(p_folded_model_spectrum_1, (*(*model_grb).p_source_spectrum), instrument_response_1, true);
	// (*p_folded_model_spectrum_1).WriteToTXT("data-file-dir/spec_model_conv_1.txt");

	// For instrument response 2
	(*model_grb).make_source_spectrum(instrument_response_2.chan_energ_min, instrument_response_2.chan_energ_max,instrument_response_2.num_chans);
	(*(*model_grb).p_source_spectrum).WriteToTXT("data-file-dir/spec_model_2.txt");
	// Spectrum * p_folded_model_spectrum_2 = NULL;
	// ConvolveSpectra(p_folded_model_spectrum_2, (*(*model_grb).p_source_spectrum), instrument_response_2, true);
	// (*p_folded_model_spectrum_2).WriteToTXT("data-file-dir/spec_model_conv_2.txt");

	*/


	//////////////////////////////////////////////////////////////////////////////////////////////////////////




	return 0;
}


