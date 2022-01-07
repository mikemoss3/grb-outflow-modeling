/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-20

This script can be used to create the spectrum for given emission data.
This script calls methods from spectrum_funcs.cpp to calculate the spectrum of the source for each energy bin.

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include "fitsio.h"


// Import Custom Libraries
#include "Spectrum.hpp"
#include "cosmology.hpp"

using namespace std;


// Define Spectrum class and member functions 

// Spectrum constructors
Spectrum::Spectrum(float energ_min, float energ_max, int num_energ_bins)
{	
	/*
	Constructor with a defined number of energy bins between specified minimum and maximum energy bounds
	
	Attributes:
	energ_min = minimum energy bound of the spectrum
	energ_max = maximum energy bound of the spectrum
	num_energ_bins = number of energy bins
	*/
	this->energ_min = energ_min; // Minimum energy of the spectrum used to create the light curve 
	this->energ_max = energ_max; // Maximum energy of the spectrum used to create the light curve 
	this->num_energ_bins = num_energ_bins; // Number of energy bins within the energy range

	// Make energy axes the correct size
	set_energ_axes_length(num_energ_bins);
	// Make energy axis, in log scale
	make_energ_axes(true);
	// Zero spectrum
	ZeroSpectrum();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

Spectrum::Spectrum(std::vector<double> spectrum_rate, std::vector<double> spectrum_unc, std::vector<float> energ_lo, std::vector<float> energ_mid, std::vector<float> energ_hi)
{
	/*
	Constructor method to explicitly define all member variables of this Spectrum object.

	Attributes:
	spectrum_rate = Vector that records the spectrum rate in each energy bin
	spectrum_unc = Vector that records the uncertainty for the spectrum rate in each energy bin
	energ_lo = Vector of energy bin lower limits
	energ_mid = Vector of the centers of each energy bin
	energ_hi = Vector of energy bin upper limits
	*/

	// Set the member variables of this Spectrum object using information from the input vectors
	this->energ_min = energ_lo.at(0);
	this->energ_max = energ_hi.back();
	this->num_energ_bins = energ_mid.size();

	// Set the size of the member variable vectors of this object.
	set_energ_axes_length(num_energ_bins);

	// Set the member variable vectors for this Spectrum object equal to the input vectors
	this->spectrum_sum = std::accumulate(spectrum_rate.begin(), spectrum_rate.end(), 0.0);
	this->spectrum_rate = spectrum_rate;
	this->spectrum_unc = spectrum_unc;
	this->energ_lo = energ_lo;
	this->energ_mid = energ_mid;
	this->energ_hi = energ_hi;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
 
// Copy constructor 
Spectrum::Spectrum(const Spectrum& tmp_spec)
{
	/*
	Copy constructor.
	*/
	Spectrum(tmp_spec.energ_min, tmp_spec.energ_max, tmp_spec.num_energ_bins);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Spectrum member functions
// Fill in energy vectors 
void Spectrum::make_energ_axes(float energ_min, float energ_max, float num_energ_bins, bool logscale)
{
	/*
	Make the energy axis from a defined number of energy bins between specified minimum and maximum energy bounds. 

	Attributes:
	energ_min = minimum energy bound of the spectrum
	energ_max = maximum energy bound of the spectrum
	num_energ_bins = number of energy bins
	logscale == boolean, indicates if the  energy axes are in log scale or not
	*/

	this->energ_min=energ_min;
	this->energ_max=energ_max;
	this->num_energ_bins=num_energ_bins;

	set_energ_axes_length(num_energ_bins);

	make_energ_axes(logscale);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void Spectrum::make_energ_axes(bool logscale)
{
	/*
	Make the energy axis from the objects current minimum energy, maximum energy, and number of energy bins. 

	Attributes:
	logscale == boolean, indicates if the  energy axes are in log scale or not
	*/
	if(logscale == true)
	{
		// float logdE = log( energ_max/energ_min ) / num_energ_bins;
		float logdE =  ( log( energ_max) - log(energ_min) ) / num_energ_bins;

		for(float i=0; i< num_energ_bins; ++i)
		{
			energ_lo.at(i) = energ_min * exp(logdE*i);
			energ_hi.at(i) = energ_min * exp(logdE*(i+1));
			energ_mid.at(i) =  energ_min * ( (exp(logdE*i) + exp(logdE*(i+1)) ) /2.);
		}
	}
	else
	{
		float dE = (energ_max-energ_min) / num_energ_bins;
		for(float i=0; i<num_energ_bins; ++i)
		{
			energ_lo.at(i) = energ_min + (dE*i);
			energ_hi.at(i) = energ_min + (dE*(i+1));
			energ_mid.at(i) = energ_min + ( (2*i + 1) /2.);
		}		
	}
	// std::cout << "Spectrum has been set back to zeros." << "\n";
	ZeroSpectrum();
	spectrum_sum = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the length of the energy axes
void Spectrum::set_energ_axes_length(int num_energ_bins)
{
	spectrum_rate.resize(num_energ_bins); // Records the spectrum rate in each energy bin
	spectrum_unc.resize(num_energ_bins); // Records the rate uncertainty in each energy bin
	energ_lo.resize(num_energ_bins); // Energy bin lower limits
	energ_mid.resize(num_energ_bins); // Center of each energy bin
	energ_hi.resize(num_energ_bins); // Energy bin upper limits
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Zero out the spectrum (this also serves to set the length of the vector)
void Spectrum::ZeroSpectrum()
{
	if(spectrum_rate.size()>0)
	{
		for(int i=0; i < num_energ_bins; ++i)
		{
			spectrum_rate.at(i) = 0.;
			spectrum_unc.at(i) = 0.;
		}

	}
	else
	{
		for(int i=0; i < num_energ_bins; ++i)
		{
			spectrum_rate.push_back(0.);
			spectrum_unc.push_back(0.);
		}

	}
	spectrum_sum = 0.;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Add an uncertainty column
void Spectrum::add_unc(std::vector<double> unc_vector)
{
	/*
	Add uncertainty values for each element in the rate vector.
	Supply an uncertainty vector designated the uncertainty for each element. The supplied vector must be same length as spectrum rate vector
	*/

	// First, check if the supplied vector is the same length as the spectrum rate vector
	if(unc_vector.size() != spectrum_rate.size())
	{
		std::cout << "The supplied uncertainty vector must be the same length as the spectrum rate vector.\n";
		return;
	}

	this->spectrum_unc = unc_vector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Add an uncertainty column
void Spectrum::add_unc(float unc_factor)
{
	/*
	Add uncertainty values for each element in the rate vector.
	Supply an uncertainty factor, uncertainty = rate * factor.
	*/

	for(size_t i=0; i<spectrum_rate.size(); ++i)
	{
		this->spectrum_unc.at(i) = spectrum_rate.at(i)*unc_factor;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Add fluctuations to the spectrum
void Spectrum::add_fluctuations()
{
	// Make a random seed generator
	std::random_device rd{};
	std::mt19937 gen{rd()};
	// For each energy bin in the spectrum, fluctuate the count rate by a Gaussian 
	for(size_t i=0; i<spectrum_rate.size(); ++i)
	{
		// Make Gaussian distribution generator
		std::normal_distribution<double> distribution(spectrum_rate.at(i), spectrum_rate.at(i)/10.);

		// Fluctuate the spectrum according to a Gaussian
		spectrum_rate.at(i) = distribution(gen);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load spectrum from text file 
void Spectrum::LoadFromTXT(std::string in_file_name)
{
	// Temporary variables used to grab file data
	std::vector<float> tmp_energ_vec;
	std::vector<double> tmp_spec_vec;
	double tmp_spectrum_sum = 0.;

	// Load in the file
	ifstream file_spec_data(in_file_name);
	string line_spec_data;
	int iter=0;
	while ( getline( file_spec_data, line_spec_data) ) 
	{
		std::istringstream iss(line_spec_data);
		string col1_val;
		string col2_val;
		if ( iss >> col1_val >> col2_val)
		{
			// Grab the necessary data 
			tmp_energ_vec.push_back( stof(col1_val) );
			tmp_spec_vec.push_back( stof(col2_val) );
			tmp_spectrum_sum += stof(col2_val);
		}
		++iter;
	}

	num_energ_bins = iter; // Number of energy bins within the energy range
	energ_min = tmp_energ_vec.at(0); // Minimum energy of the spectrum used to create the light curve 
	energ_max = tmp_energ_vec.back(); // Maximum energy of the spectrum used to create the light curve 

	// Make energy vectors with correct sizes
	set_energ_axes_length(num_energ_bins);

	// Set the bin center energy vector to the data file's energy bin column 
	energ_mid = tmp_energ_vec;
	// Set the spectrum rate vector to the data file's spectrum column
	spectrum_rate = tmp_spec_vec;
	// Set the spectrum summation to the sum of the data file spectrum
	spectrum_sum = tmp_spectrum_sum;
	
	// Close files and free memory 
	file_spec_data.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load spectrum from FITS file
int Spectrum::LoadFromFITS(std::string in_file_name)
{
	/*
	This method is used to load a spectrum from a FITS file defined by in_file_name.

	The given FITS file should be structured in the following way:
	One HDU containing an EBOUNDS BinTable, which indicates the EMIN and EMAX for each instrument CHANNEL.
	And one HDU containing a Spectrum BinTable, which contains CHANNEL, RATE, and ERROR for the spectrum.
	*/

	// Declare all variable types 
	fitsfile *fptr; // FITS file object 
	int status = 0; // Error status
	char comment[FLEN_COMMENT], extname[FLEN_VALUE]; // FITS character[size]
	double nulldouble = 0.0; // Flag used to ignore null elements 
	long num_rows,number_chans;
	int irow, colnum, anynull, num_hdus, ihdu, hdutype; 
	int spectrum_ihdu, ebounds_ihdu;

	// Open FITS file
	int n = in_file_name.length() + 1;
	char dumn_char_str[n];
	strcpy(dumn_char_str, in_file_name.c_str());
	fits_open_file(&fptr, dumn_char_str, READONLY, &status);

	// Get number of HDUs 
	fits_get_num_hdus(fptr, &num_hdus, &status);
	// Skip hdu=1, this is the primary HDU
	for (ihdu=2; ihdu<=num_hdus; ihdu++)
	{
		// Move to current HDU 
		fits_movabs_hdu(fptr, ihdu, &hdutype, &status);
		// We only want Binary Table data units
		if (hdutype == BINARY_TBL) 
		{
			// Extract the extension name and the array name.
			fits_read_key(fptr, TSTRING, "EXTNAME", extname, comment, &status);
			
			// Check for errors 
			if (status){fits_report_error(stderr, status); return status;}

			// Keep track of which HDU corresponds to which extension 
			if (strcmp(extname, "SPECTRUM")==0){spectrum_ihdu = ihdu;}
			if (strcmp(extname, "EBOUNDS")==0){ebounds_ihdu = ihdu;}
		}
	}

	/* Read Ebounds  */
	// Move to Ebounds table 
	fits_movabs_hdu(fptr, ebounds_ihdu, &hdutype, &status);

	// Get number of rows
	fits_get_num_rows(fptr, &number_chans, &status);

	// Temporary variables to store element values
	// float tmp_CHAN=0.;
	float tmp_chan_energ_lo = 0.;
	float tmp_chan_energ_hi = 0.;
	std::vector<float> tmp_energ_lo(number_chans);
	std::vector<float> tmp_energ_mid(number_chans);
	std::vector<float> tmp_energ_hi(number_chans);

	// Column names :
	char colname_chan_energ_lo[] = "E_MIN";
	char colname_chan_energ_hi[] = "E_MAX";

	/* Read rows of channel data*/
	for (irow=1; irow<=number_chans; irow++) 
	{
		// Select chan_energ_lo column
		fits_get_colnum(fptr, CASEINSEN, colname_chan_energ_lo, &colnum, &status);
		// Read column data at row i
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_chan_energ_lo, &anynull, &status);
		// Store column data in appropriate vector
		tmp_energ_lo.at(irow-1) = tmp_chan_energ_lo;

		// chan_energ_hi
		fits_get_colnum(fptr, CASEINSEN, colname_chan_energ_hi, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_chan_energ_hi, &anynull, &status);
		tmp_energ_hi.at(irow-1) = tmp_chan_energ_hi;

		// Calculate the center of the bin
		tmp_energ_mid.at(irow-1) =  (tmp_chan_energ_lo+tmp_chan_energ_hi)/2. ;
	}

	/* Read spectrum data  */
	// Move to spectrum table 
	fits_movabs_hdu(fptr, spectrum_ihdu, &hdutype, &status);

	// Get number of rows
	fits_get_num_rows(fptr, &num_rows, &status);

	// Temporary variables to store element values
	float tmp_spectrum_rate = 0.;
	std::vector<double> tmp_spectrum_rate_vector(number_chans);
	double tmp_spectrum_sum = 0.;

	// Column names :
	char colname_count_rate[] = "RATE";

	/* Read rows of spectrum data*/
	for (irow=1; irow<=num_rows; irow++) 
	{
		// Select count rate Column
		fits_get_colnum(fptr, CASEINSEN, colname_count_rate, &colnum, &status);
		// Read column data at row i
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_spectrum_rate, &anynull, &status);
		// Store column data in appropriate vector
		tmp_spectrum_rate_vector.at(irow - 1) = tmp_spectrum_rate;
		tmp_spectrum_sum += tmp_spectrum_rate;

	}

	energ_min = tmp_energ_lo.at(0);
	energ_max = tmp_energ_hi.back();
	num_energ_bins = number_chans;

	set_energ_axes_length(num_energ_bins);

	energ_lo = tmp_energ_lo;
	energ_mid = tmp_energ_mid;
	energ_hi = tmp_energ_hi;

	spectrum_rate = tmp_spectrum_rate_vector;
	spectrum_sum = tmp_spectrum_sum;

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

void Spectrum::WriteToTXT(std::string out_file_name)
{
	// Write the light curve to a text file
	ofstream spec_file; // Construct file 
	spec_file.open(out_file_name); // Open text file with this name
	int i=0;
	// For each time bin, write the time and count rate to the file.
	while ( i < num_energ_bins)
	{
		spec_file << energ_mid.at(i);
		spec_file << " ";
		spec_file << spectrum_rate.at(i);
		spec_file << " ";
		spec_file << spectrum_unc.at(i);
		spec_file << "\n";		
		i++;
	}
	spec_file.close(); // Close file
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
int Spectrum::WriteToFITS(std::string out_file_name)
{
	// Declare all variable types 
	fitsfile *fptr; // FITS file object 
	int status = 0; // Error status
	char comment[FLEN_COMMENT], extname[FLEN_VALUE]; // FITS character[size]
	double nulldouble = 0.0; // Flag used to ignore null elements 
	long num_rows,number_chans;
	int irow, colnum, anynull, num_hdus, ihdu, hdutype; 
	int spectrum_ihdu, ebounds_ihdu;

	// Declaring character array
	std::string s = "!"+out_file_name;
    int n = s.length();
    char char_array[n + 1];
 
    // Copying the contents of the string to char array
    strcpy(char_array, s.c_str());

	// Open FITS file
	// fits_create_file(fitsfile **fptr, char *filename, > int *status)
	fits_create_file(&fptr, char_array, &status);

	// Get number of HDUs 
	fits_get_num_hdus(fptr, &num_hdus, &status);
	// Skip hdu=1, this is the primary HDU
	for (ihdu=2; ihdu<=num_hdus; ihdu++)
	{
		// Move to current HDU 
		fits_movabs_hdu(fptr, ihdu, &hdutype, &status);
		// We only want Binary Table data units
		if (hdutype == BINARY_TBL) 
		{
			// Extract the extension name and the array name.
			fits_read_key(fptr, TSTRING, "EXTNAME", extname, comment, &status);
			
			// Check for errors 
			if (status){fits_report_error(stderr, status); return status;}

			// Keep track of which HDU corresponds to which extension 
			if (strcmp(extname, "SPECTRUM")==0){spectrum_ihdu = ihdu;}
			if (strcmp(extname, "EBOUNDS")==0){ebounds_ihdu = ihdu;}
		}
	}

	fits_close_file(fptr, &status);
	return status;
}
*/

