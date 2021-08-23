/*
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-08-11

Definitions of the Response class 

*/

// Import Standard Libraries
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "fitsio.h"

using namespace std;

// Include custom libraries and headers
#include "Response.hpp"

// Define Response class and member functions 

// Spectrum constructor
Response::Response()
{

}
 
// Response member functions
// Set photon energy axis 
void Response::set_phot_energ(float phot_energ_min, float phot_energ_max, int num_phot_bins, bool logscale)
{

	this->phot_energ_min = phot_energ_min;
	this->phot_energ_max = phot_energ_max;
	this->num_phot_bins = num_phot_bins;

	if(logscale == true)
	{
		float logdE = log( phot_energ_max/phot_energ_min ) / num_phot_bins;
		for(int i=0; i<num_phot_bins;i++)
		{
			phot_energ_lo.push_back( phot_energ_min + exp(logdE*i) );
			phot_energ_hi.push_back( phot_energ_min + exp(logdE*(i+1)) );
			phot_energ_mid.push_back( (phot_energ_lo.at(i) + phot_energ_hi.at(i)) /2 );
		}
	}
	else
	{
		float dE = (phot_energ_max - phot_energ_min) / num_phot_bins;
		for(int i=0; i<num_phot_bins;i++)
		{
			phot_energ_lo.push_back( phot_energ_min + (dE*i) );
			phot_energ_hi.push_back( phot_energ_min + (dE*(i+1)) );
			phot_energ_mid.push_back( (phot_energ_lo.at(i) + phot_energ_hi.at(i)) /2 );
		}		
	}
	std::cout << "Response matrix has been reset to zeros." << "\n";
	make_empty_resp();
	set_n_grp();
	set_f_chan();

}
// Set instrument channel energy axis
void Response::set_chan_energ(float chan_energ_min, float chan_energ_max, int num_chans, bool logscale)
{
	this->chan_energ_min = chan_energ_min;
	this->chan_energ_max = chan_energ_max;
	this->num_chans = num_chans;

	if(logscale == true)
	{
		float logdE = log( chan_energ_min/chan_energ_min ) / num_chans;
		for(int i=0; i<num_chans;i++)
		{
			chan_energ_lo.push_back( chan_energ_min + exp(logdE*i) );
			chan_energ_hi.push_back( chan_energ_min + exp(logdE*(i+1)) );
			chan_energ_mid.push_back( (chan_energ_lo.at(i) + chan_energ_lo.at(i)) /2 );
		}
	}
	else
	{
		float dE = (chan_energ_min - chan_energ_min) / num_chans;
		for(int i=0; i<num_chans;i++)
		{
			chan_energ_lo.push_back( chan_energ_min + (dE*i) );
			chan_energ_hi.push_back( chan_energ_min + (dE*(i+1)) );
			chan_energ_mid.push_back( (chan_energ_lo.at(i) + chan_energ_lo.at(i)) /2 );
		}		
	}
	std::cout << "Response matrix has been reset to zeros." << "\n";
	make_empty_resp();
	set_n_chan();
}

void Response::set_n_grp()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		n_grp.push_back(1);
	}
}
void Response::set_f_chan()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		f_chan.push_back(0);
	}
}
void Response::set_n_chan()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		n_chan.push_back(num_chans);
	}
}

// If the shape of the response matrix is changed, the response matrix is reset to zeros.
void Response::make_empty_resp()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		std::vector<float> tmp_vec(num_chans);
		for(int j=0; j<num_chans; j++)
		{
			tmp_vec.at(j) = 0.;
		}
		prob_matrix.push_back(tmp_vec);
	}
}
// Make identity matrix
void Response::Identity()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		std::vector<float> tmp_vec(num_chans);
		for(int j=0; j<num_chans; j++)
		{
			if(i==j)
			{
				tmp_vec.at(i) = 1;
			}
			else
			{
				tmp_vec.at(j) = 0;
			}
		}
		prob_matrix.push_back(tmp_vec);
	}
}

// Decrease as 1/DeltaE^alpha from E_true
void Response::OverDeltaE(float alpha)
{
	float col_sum; // Used for normalizing the column
	float tmp_val;
	for(int i=0; i<num_phot_bins; i++)
	{
		col_sum = 0;
		std::vector<float> tmp_vec(num_chans);
		for(int j=0; j<num_chans; j++)
		{
			tmp_val = 1/(1+pow(abs( chan_energ_mid.at(j) - phot_energ_mid.at(i)),alpha));
			tmp_vec.at(j) = tmp_val;
			col_sum += tmp_val; 
		}
		// Normalize the 
		for(int j=0; j<num_chans; j++)
		{
			tmp_vec.at(j) /= col_sum;
		}
		prob_matrix.push_back(tmp_vec);
	}

}

// Load response matrix from file
int Response::LoadRespFromFile(std::string file_name)
{
	// Declare all variable types 
	fitsfile *fptr; // FITS file object 
	int status = 0; // Error status
	char comment[FLEN_COMMENT], extname[FLEN_VALUE]; // FITS character[size]
	double nulldouble = 0.0; // Flag used to ignore null elements 
	long num_rows,number_chans;
	int irow, colnum, anynull, num_hdus, ihdu, hdutype; 
	int response_ihdu, ebounds_ihdu;

	// Open FITS file
	int n = file_name.length() + 1;
	char dumn_char_str[n];
	strcpy(dumn_char_str, file_name.c_str());
	// dumn_char_str = file_name;
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
			if (strcmp(extname, "SPECRESP MATRIX")==0){response_ihdu = ihdu;}
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
	float tmp_chan_energ_lo=0.;
	float tmp_chan_energ_hi=0.;

	// Column names :
	// char colname_CHAN[] = "CHANNEL";
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
		chan_energ_lo.push_back(tmp_chan_energ_lo);

		// chan_energ_hi
		fits_get_colnum(fptr, CASEINSEN, colname_chan_energ_hi, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_chan_energ_hi, &anynull, &status);
		chan_energ_hi.push_back(tmp_chan_energ_hi);

		// Calculate the center of the bin
		chan_energ_mid.push_back( (tmp_chan_energ_lo+tmp_chan_energ_hi)/2 );
	}

	/* Read response matrix  */
	// Move to Response table 
	fits_movabs_hdu(fptr, response_ihdu, &hdutype, &status);

	// Get number of rows
	fits_get_num_rows(fptr, &num_rows, &status);

	// Temporary variables to store element values
	float tmp_phot_energ_lo=0.;
	float tmp_phot_energ_hi=0.;
	float tmp_n_grp=0.;
	float tmp_f_chan=0.;
	float tmp_n_chan=0.;
	float tmp_prob_matrix[number_chans];

	// Column names :
	char colname_phot_energ_lo[] = "ENERG_LO";
	char colname_phot_energ_hi[] = "ENERG_HI";
	char colname_n_grp[] = "N_GRP";
	char colname_f_chan[] = "F_CHAN";
	char colname_n_chan[] = "N_CHAN";
	char colname_prob_matrix[] = "MATRIX";

	/* Read rows of spectrum data*/
	for (irow=1; irow<=num_rows; irow++) 
	{
		// Select phot_energ_lo Column
		fits_get_colnum(fptr, CASEINSEN, colname_phot_energ_lo, &colnum, &status);
		// Read column data at row i
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_phot_energ_lo, &anynull, &status);
		// Store column data in appropriate vector
		phot_energ_lo.push_back(tmp_phot_energ_lo);

		// phot_energ_hi
		fits_get_colnum(fptr, CASEINSEN, colname_phot_energ_hi, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_phot_energ_hi, &anynull, &status);
		phot_energ_hi.push_back(tmp_phot_energ_hi);

		// Calculate ENERG_MID
		phot_energ_mid.push_back( (tmp_phot_energ_lo+tmp_phot_energ_hi)/2 );

		// N_GRP
		fits_get_colnum(fptr, CASEINSEN, colname_n_grp, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_n_grp, &anynull, &status);
		n_grp.push_back(tmp_n_grp);

		// F_CHAN
		fits_get_colnum(fptr, CASEINSEN, colname_f_chan, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_f_chan, &anynull, &status);
		f_chan.push_back(tmp_f_chan);

		// N_CHAN
		fits_get_colnum(fptr, CASEINSEN, colname_n_chan, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_n_chan, &anynull, &status);
		n_chan.push_back(tmp_n_chan);

		// MATRIX
		fits_get_colnum(fptr, CASEINSEN, colname_prob_matrix, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, number_chans, &nulldouble, &tmp_prob_matrix, &anynull, &status);
		std::vector<float> tmp_vec(tmp_prob_matrix, tmp_prob_matrix + number_chans );
		prob_matrix.push_back(tmp_vec);

	}

	phot_energ_min = phot_energ_lo.at(0);
	phot_energ_max = phot_energ_hi.back();
	num_phot_bins = num_rows;

	chan_energ_min = chan_energ_lo.at(0);
	chan_energ_max = chan_energ_hi.back();
	num_chans = number_chans;

	fits_close_file(fptr, &status);

	/* print any error messages */
	if (status)
	{
		fits_report_error(stderr, status);
		return status;
	}

	return status;
}

// Fold a spectrum with a response matrix
int Response::ConvolveSpectrum(Spectrum * folded_spectrum, const Spectrum & source_spectrum)
{
	// Fill in the energy vector
	(*folded_spectrum).energ_lo = chan_energ_lo;
	(*folded_spectrum).energ_mid = chan_energ_mid;
	(*folded_spectrum).energ_hi = chan_energ_hi;

	// Check if the number of source photon energy bins is the same as the number of
	// response photon energy bins 
	if(source_spectrum.num_energ_bins != num_phot_bins)
	{
		std::cout << "Source spectrum and Response Matrix do not have same number of photon bins.\n";
		return 1;
	}

	/* Convolve source spectrum with instrument response matrix */
	double tmp_col_sum = 0; // Keeps track of matrix multiplication sum 
	// For each instrument energy channel
	for( int i=0; i < num_chans; i++)
	{
		tmp_col_sum = 0; // Reset column sum to zero

		// For for each photon energy bin
		for( int j=0 ; j < source_spectrum.num_energ_bins; j++)
		{
			tmp_col_sum += prob_matrix.at(j).at(i) * source_spectrum.spectrum_rate.at(j);
		}
		
		// Write the summation to the corresponding energy channel of the folded spectrum 
		(*folded_spectrum).spectrum_rate.at(i) = tmp_col_sum;
		// Increase spectrum sum 
		(*folded_spectrum).spectrum_sum += tmp_col_sum;
	}

	return 0;
}