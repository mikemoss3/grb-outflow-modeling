/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-12-01

This object holds Time Tagged Event (TTE) data. 
Each event (or instrument count) has an associated time stamp and measured PHA. The PHA is the instrument energy channel the count was observed in.
For any instrument there is a relation between the energy channel an energy bound associated with the channel.   
*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// Import Custom Libraries
#include "Spectrum.hpp"
#include "LightCurve.hpp"
#include "TTEs.hpp"

using namespace std;

// TTEs constructor
TTEs::TTEs()
{
	/*
	TTE object constructor
	*/

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load GRB spectrum and light curve from TTE data designated by the given file name
int TTEs::LoadTTEData(std::string file_name)
{
	/*
	This method is used to load a spectrum and light curve from a FITS file defined by file_name.

	The given FITS file should be of the standard form of the TTE data reported by the Fermi/GBM instrument.
	*/

	// Declare all variable types 
	fitsfile *fptr; // FITS file object 
	int status = 0; // Error status
	double nulldouble = 0.0; // Flag used to ignore null elements 
	long num_chans,num_evts;
	int irow, colnum, anynull, hdutype; 

	// Open FITS file
	// First, we need to convert the file_name from a string to a character array 
	int n = file_name.length() + 1; // Find length of string
	char dumn_char_str[n]; // Make empty character array
	strcpy(dumn_char_str, file_name.c_str()); // Copy string to character array
	fits_open_file(&fptr, dumn_char_str, READONLY, &status); // Open FILE. The file will be pointed to by fptr
	
	/* Load EBOUNDS BinTable*/
	// Move to the HDU which stores the EBOUNDS data
	fits_movabs_hdu(fptr, 2, &hdutype, &status);

	// Get number of rows, this is the number of events
	fits_get_num_rows(fptr, &num_chans, &status);


	// Relevant column names:
	char colname_CHANNEL[] = "CHANNEL";
	char colname_E_MIN[] = "E_MIN";	
	char colname_E_MAX[] = "E_MAX";	

	// Temporary variables to store element values
	float tmp_CHANNEL = 0.;
	float tmp_E_MIN = 0.;
	float tmp_E_MAX = 0.;

	
	// Resize ebounds vectors (i.e., the TTE CHANS, ELO, and EHI vectors)
	resize_ebounds(num_chans);


	/* Read rows of channel data*/
	for (irow=1; irow<=num_chans; irow++) 
	{
		// Select CHANNEL column
		fits_get_colnum(fptr, CASEINSEN, colname_CHANNEL, &colnum, &status);
		// Read column data at row i
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_CHANNEL, &anynull, &status);
		// Store column data in appropriate vector
		TTE_CHANS.at(irow-1) = tmp_CHANNEL;

		// E_MIN column
		fits_get_colnum(fptr, CASEINSEN, colname_E_MIN, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_E_MIN, &anynull, &status);
		TTE_ELO.at(irow-1) = tmp_E_MIN;

		// E_MAX column
		fits_get_colnum(fptr, CASEINSEN, colname_E_MAX, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_E_MAX, &anynull, &status);
		TTE_EHI.at(irow-1) = tmp_E_MAX;
	}


	/* Read event times and PHA values*/
	// Move to the HDU which stores the event data
	fits_movabs_hdu(fptr, 3, &hdutype, &status);

	// Get number of rows, this is the number of events
	fits_get_num_rows(fptr, &num_evts, &status);

	// Extract the trigger time of the event. This will be subtracted from the time array to change from "detector time" to "time since trigger"
	double trig_time;
	char comment[FLEN_COMMENT];
	char char_trig_time[] = "TRIGTIME";
	fits_read_key(fptr, TDOUBLE, char_trig_time, &trig_time, comment, &status);

	// Relevant column names:
	char colname_TIME[] = "TIME";
	char colname_PHA[] = "PHA";	

	// Temporary variables to store element values
	double tmp_TIME = 0.;
	float tmp_PHA = 0.;

	// Resize the events vectors (i.e., the TTE TIME and PHA vectors)
	resize_events(num_evts);

	/* Read rows of channel data*/
	for (irow=1; irow<=num_evts; irow++) 
	{
		// Select TIME column
		fits_get_colnum(fptr, CASEINSEN, colname_TIME, &colnum, &status);
		// Read column data at row i
		fits_read_col(fptr, TDOUBLE, colnum, irow, 1, 1, &nulldouble, &tmp_TIME, &anynull, &status);
		// Store column data in appropriate vector
		TTE_TIME.at(irow-1) = tmp_TIME-trig_time;

		// PHA column
		fits_get_colnum(fptr, CASEINSEN, colname_PHA, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_PHA, &anynull, &status);
		TTE_PHA.at(irow-1) = tmp_PHA;
	}


	fits_close_file(fptr, &status);

	/* print any error messages */
	if (status)
	{
		fits_report_error(stderr, status);
		return status;
	}

	return status;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Resize the events vectors (i.e., the TTE TIME and PHA vectors)
void TTEs::resize_events(int num_evts)
{
	TTE_TIME.resize(num_evts);
	TTE_PHA.resize(num_evts);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Resize ebounds vectors (i.e., the TTE CHANS, ELO, and EHI vectors)
void TTEs::resize_ebounds(int num_chans)
{
	TTE_CHANS.resize(num_chans);
	TTE_ELO.resize(num_chans);
	TTE_EHI.resize(num_chans);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////


// Make spectrum from TTE data
void TTEs::make_spectrum(Spectrum * & p_spectrum)
{
	/*
	Method to make a spectrum from the TTE data and store it in the Spectrum object passed to this method.
	
	The process is as follows: each event in the TTE data is has a PHA reported, which indicates the instrument channel the event was observed in. 
	Using the relation between the instrument energy channels and the associated energy bounds, we can add each event to the correct energy bin and 
	therefore create a spectrum.

	Attributes:
	p_spectrum = pointer to a Spectrum object, this will be used to store the constructed spectrum
	*/

	/* Make the middle energy bin vector */
	// Make a vector to store the values
	std::vector<float> TTE_EMID(TTE_CHANS.size());
	// Calculate the middle energy bin values
	for(size_t i=0; i<TTE_CHANS.size(); ++i )
	{
		TTE_EMID.at(i) = (TTE_ELO.at(i) + TTE_EHI.at(i))/2.;
	}

	/* Calculate count in each energy bin */
	// Create spectrum rate array and rate uncertainty array
	std::vector<double> spectrum_rate(TTE_CHANS.size());
	std::vector<double> spectrum_unc(TTE_CHANS.size());
	// For each event, check which instrument channel it is in
	for(size_t i=0; i<TTE_PHA.size(); ++i)
	{
		// Find the corresponding energy channel to the event PHA 
		for(size_t j=0; j<TTE_CHANS.size(); ++j)
		{
			// If the channel matches the PHA of the event
			// And if the event is within the designated time interval. Add the event to the spectrum
			if( TTE_PHA.at(i) == TTE_CHANS.at(j) )
			{
				// Add a count to the corresponding spectrum bin
				spectrum_rate.at(j) += 1;
			}
		}
	}

	// Calculate spectrum uncertainty array
	for(size_t i = 0; i < TTE_CHANS.size(); ++i )
	{
		spectrum_unc.at(i) = sqrt(spectrum_rate.at(i));
	}

	/* Make Spectrum object from TTE data*/
	// Construct a Spectrum object 
	p_spectrum = new Spectrum(spectrum_rate, spectrum_unc, TTE_ELO, TTE_EMID, TTE_EHI);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make spectrum from TTE data
void TTEs::make_spectrum(Spectrum * & p_spectrum,float TSTART, float TEND, float EMIN, float EMAX)
{
	/*
	Method to make a spectrum from the TTE data and store it in the Spectrum object passed to this method.
	
	The process is as follows: each event in the TTE data is has a PHA reported, which indicates the instrument channel the event was observed in. 
	Using the relation between the instrument energy channels and the associated energy bounds, we can add each event to the correct energy bin and 
	therefore create a spectrum.

	Attributes:
	p_spectrum = pointer to a Spectrum object, this will be used to store the constructed spectrum
	TSTART, TEND = time interval bounds, all events that occur within this time interval will be used to make the spectrum
	EMIN, EMAX = energy band bounds, all events that occur within this energy band will be used to make the spectrum
	*/


	// Find the indices of the new energy band
	// EMIN
	int ind_emin = 0;
	while( TTE_ELO.at(ind_emin) < EMIN ){ ++ind_emin; }
	// EMAX
	int ind_emax = 0;
	while( TTE_EHI.at(ind_emax) < EMAX ){ ++ind_emax; }

	// New number of energy bins
	int new_num_chans = ind_emax - ind_emin;
	// Slice energy vectors to confine them to just the energy band
	std::vector<float> new_TTE_CHANS = slicing(TTE_CHANS, ind_emin, ind_emax);
	std::vector<float> new_TTE_ELO = slicing(TTE_ELO, ind_emin, ind_emax);
	std::vector<float> new_TTE_EHI = slicing(TTE_EHI, ind_emin, ind_emax);


	/* Make the middle energy bin vector */
	// Make a vector to store the values
	std::vector<float> TTE_EMID(new_num_chans);
	// Calculate the middle energy bin values
	for(int i=0; i<new_num_chans; ++i )
	{
		TTE_EMID.at(i) = (new_TTE_ELO.at(i) + new_TTE_EHI.at(i))/2.;
	}

	/* Calculate count in each energy bin */
	// Create spectrum rate array and rate uncertainty array
	std::vector<double> spectrum_rate(new_num_chans);
	std::vector<double> spectrum_unc(new_num_chans);
	// For each event, check which instrument channel it is in
	for(size_t i=0; i<TTE_PHA.size(); ++i)
	{
		// Find the corresponding energy channel to the event PHA 
		for(int j=0; j<new_num_chans; ++j)
		{
			// If the channel matches the PHA of the event
			// And if the event is within the designated time interval. Add the event to the spectrum
			if( (TTE_PHA.at(i) == new_TTE_CHANS.at(j) ) & (TTE_TIME.at(i) >= TSTART) & ( TTE_TIME.at(i) <= TEND) )
			{
				// Add a count to the corresponding spectrum bin
				spectrum_rate.at(j) += 1;
			}
		}
	}

	for(int i = 0; i < new_num_chans; ++i )
	{
		// Calculate spectrum uncertainty array
		spectrum_unc.at(i) = sqrt(spectrum_rate.at(i));
		// Normalize spectrum and the uncertainty by exposure time
		spectrum_rate.at(i) /= (TEND - TSTART);
		spectrum_unc.at(i) /= (TEND - TSTART);

	}

	/* Make Spectrum object from TTE data*/
	// Construct a Spectrum object 
	p_spectrum = new Spectrum(spectrum_rate, spectrum_unc, new_TTE_ELO, TTE_EMID, new_TTE_EHI);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Slices vectors between the new start and stop points
vector<float> TTEs::slicing(vector<float>& arr, int new_start, int new_end)
{ 
    // Starting and Ending iterators
    auto start = arr.begin() + new_start;
    auto end = arr.begin() + new_end + 1;
 
    // To store the sliced vector
    vector<float> result(new_end - new_start + 1);
 
    // Copy vector using copy function()
    copy(start, end, result.begin());
 
    // Return the final sliced vector
    return result;
}

// Slices vectors between the new start and stop points
vector<double> TTEs::slicing(vector<double>& arr, int new_start, int new_end)
{ 
    // Starting and Ending iterators
    auto start = arr.begin() + new_start;
    auto end = arr.begin() + new_end + 1;
 
    // To store the sliced vector
    vector<double> result(new_end - new_start + 1);
 
    // Copy vector using copy function()
    copy(start, end, result.begin());
 
    // Return the final sliced vector
    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make light curve from TTE data
void TTEs::make_light_curve(LightCurve * & p_light_curve, float dt)
{
	/*
	Method to make a light curve from the TTE data and store it in the LightCurve object passed to this method.

	The process is as follows: each event in the TTE data has a time stamp. Using a given time bin size, 
	all events with a time stamp within the time bin bounds can be summed to create a light curve.

	Attributes:
	p_light_curve = pointer to a LightCurve object. The calculated light curve will overwrite any information of the object.
	dt = time bin size
	*/

	/* Calculate count rate in each time bin */
	// Using the first event time and final event time, find the number of time bins of size dt
	int num_t_bins = (TTE_TIME.back() - TTE_TIME.at(0))/dt;

	// Define vectors to record time bins and count rate
	std::vector<double> lc_time(num_t_bins);
	std::vector<double> lc_rate(num_t_bins); 

	// Calculate the rate in each time bin and fill in the time values
	int j = 0;
	for(float i=0; i<num_t_bins; ++i)
	{
		lc_time.at(i) = TTE_TIME.at(0) + (dt*i);

		while( TTE_TIME.at(j) < lc_time.at(i) + dt)
		{
			lc_rate.at(i) += 1;
			++j;
		}
	}

	/* Make Light Curve object from TTE data*/
	p_light_curve = new LightCurve(lc_time, lc_rate);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make light curve from TTE data within a certain time interval and energy range
void TTEs::make_light_curve(LightCurve * & p_light_curve, float TSTART, float TEND, float EMIN, float EMAX,float dt)
{
	/*
	Method to make a light curve from the TTE data and store it in the LightCurve object passed to this method.

	The process is as follows: each event in the TTE data has a time stamp. Using a given time bin size, 
	all events with a time stamp within the time bin bounds can be summed to create a light curve.

	Attributes:
	p_light_curve = pointer to a LightCurve object. The calculated light curve will overwrite any information of the object.
	dt = time bin size
	TSTART, TEND = time interval bounds, all events that occur within this time interval will be used to make the spectrum
	EMIN, EMAX = energy band bounds, all events that occur within this energy band will be used to make the spectrum
	*/

	// Find the indices of the new time interval
	int ind_tstart = 0;
	while( TTE_TIME.at(ind_tstart) < TSTART ){ ++ind_tstart; }
	int ind_tend = 0;
	while( TTE_TIME.at(ind_tend) < TEND ){ ++ind_tend; }

	// Number of events within the time interval
	int new_num_evts = (ind_tend - ind_tstart);
	// Slice time vector and event vector to within the time interval
	std::vector<double> new_TTE_TIME = slicing(TTE_TIME, ind_tstart, ind_tend);
	std::vector<float> new_TTE_PHA = slicing(TTE_PHA, ind_tstart, ind_tend);
	// Number of time bins
	int num_t_bins = (TEND - TSTART)/dt;

	// Define vectors to record time bins and count rate
	std::vector<double> lc_time(num_t_bins);
	std::vector<double> lc_rate(num_t_bins); 

	// Calculate the rate in each time bin and fill in the time values
	int j = 0;
	for(float i=0; i<num_t_bins; ++i)
	{
		lc_time.at(i) = new_TTE_TIME.at(0) + (dt*i);

		while( (new_TTE_TIME.at(j) < lc_time.at(i) + dt) & (j < new_num_evts) )
		{
			// Only count events that occur within the designated energy band
			if ( ( TTE_EHI.at( TTE_PHA.at(j) ) > EMIN ) | ( TTE_ELO.at( TTE_PHA.at(j) ) < EMAX ) )
				lc_rate.at(i) += 1;
			++j;
		}
	}

	/* Make Light Curve object from TTE data*/
	p_light_curve = new LightCurve(lc_time, lc_rate);

}