
/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2023-06-13

Function to create and use interactive Spectra and Light Curves

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

int interplot(int argc, char const *argv[])
{


	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* 
	Called during interactive plotting manipulation

	argv[1] = indicates if the update changes the plotted time interval or the energy band

	argv[2] = indicates if all the components should be calculated or just the summed total 
	argv[3] = indicates if a log scale should be used to calculate the light curve
	argv[4], argv[5] = input time interval
	argv[6], argv[7] = input energy interval

	*/

	// Directory path to data storage directory
	std::string dir="../files-data/synthetic-data/interplot-data/";

	if(argv[1]!=NULL)
	{

		bool logscale = false;
		if( (strcmp(argv[3],"True")==0) | (strcmp(argv[3],"true")==0) | (strcmp(argv[3],"1")==0) )
		{
			logscale = true;
		}
		
		if( strcmp(argv[1], "timechange")==0 )
		{

			float tmin = stof(argv[4]);
			float tmax = stof(argv[5]);

			float energ_min = stof(argv[6]);
			float energ_max = stof(argv[7]);
			float num_energ_bins = 10.*log10(energ_max/energ_min);

			SynthGRB synth_grb = SynthGRB();
			synth_grb.LoadJetParamsFromTXT("files-input/jet-params.txt");
			synth_grb.SimulateJetDynamics();

			// Total spectrum
			synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax);
			synth_grb.WriteSpectrumToTXT(dir+"interplot_spectrum.txt");

			// Make component spectra
			if(strcmp(argv[2], "true")==0)
			{
				// Thermal Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "TH");
				synth_grb.WriteSpectrumToTXT(dir+"interplot_spectrum_TH.txt");

				// Internal Shock Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "IS");
				synth_grb.WriteSpectrumToTXT(dir+"interplot_spectrum_IS.txt");

				// Forward Shock Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "FS");
				synth_grb.WriteSpectrumToTXT(dir+"interplot_spectrum_FS.txt");

				// Reverse Shock Component
				synth_grb.make_source_spectrum(energ_min, energ_max, num_energ_bins, tmin, tmax, "RS");
				synth_grb.WriteSpectrumToTXT(dir+"interplot_spectrum_RS.txt");
			}
			
			return 0; 
		}
		if( strcmp(argv[1], "energychange")==0 )
		{

			float tmin = stof(argv[4]);
			float tmax = stof(argv[5]);
			float dt = stof(argv[6]);

			float energ_min = stof(argv[7]);
			float energ_max = stof(argv[8]);

			SynthGRB synth_grb = SynthGRB();
			synth_grb.LoadJetParamsFromTXT("files-input/jet-params.txt");
			synth_grb.SimulateJetDynamics();

			synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "all", logscale);
			synth_grb.WriteLightCurveToTXT(dir+"interplot_light_curve.txt");
			

			// Make component spectra and light curves
			if(strcmp(argv[2], "true")==0)
			{
				// Plot the other light curves components 
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "TH", logscale);
				synth_grb.WriteLightCurveToTXT(dir+"interplot_light_curve_TH.txt");
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "IS", logscale);
				synth_grb.WriteLightCurveToTXT(dir+"interplot_light_curve_IS.txt");
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "FS", logscale);
				synth_grb.WriteLightCurveToTXT(dir+"interplot_light_curve_FS.txt");
				synth_grb.make_source_light_curve(energ_min, energ_max, tmin, tmax, dt, "RS", logscale);
				synth_grb.WriteLightCurveToTXT(dir+"interplot_light_curve_RS.txt");				
			}

			return 0; 
		}
	}

	return 0;
}