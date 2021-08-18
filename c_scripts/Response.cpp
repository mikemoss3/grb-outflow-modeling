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
void Response::set_ENERG(float E_phot_min, float E_phot_max, int num_phot_bins, bool logscale)
{
	if(logscale == true)
	{
		float logdE = log( E_phot_max/E_phot_min ) / num_phot_bins;
		for(int i=0; i<num_phot_bins;i++)
		{
			ENERG_LO.push_back( E_phot_min + exp(logdE*i) );
			ENERG_HI.push_back( E_phot_min + exp(logdE*(i+1)) );
			ENERG_MID.push_back( (ENERG_LO[i] + ENERG_HI[i]) /2 );
		}
	}
	else
	{
		float dE = (E_phot_max - E_phot_min) / num_phot_bins;
		for(int i=0; i<num_phot_bins;i++)
		{
			ENERG_LO.push_back( E_phot_min + (dE*i) );
			ENERG_HI.push_back( E_phot_min + (dE*(i+1)) );
			ENERG_MID.push_back( (ENERG_LO[i] + ENERG_HI[i]) /2 );
		}		
	}
	std::cout << "Response matrix has been reset to zeros." << "\n";
	make_empty_resp();
	set_N_GRP();
	set_F_CHAN();

}
// Set instrument channel energy axis
void Response::set_ECHAN(float E_chan_min, float E_chan_max, int num_chans, bool logscale)
{
	if(logscale == true)
	{
		float logdE = log( E_chan_max/E_chan_min ) / num_chans;
		for(int i=0; i<num_chans;i++)
		{
			ENERG_LO.push_back( E_chan_min + exp(logdE*i) );
			ENERG_HI.push_back( E_chan_min + exp(logdE*(i+1)) );
			ENERG_MID.push_back( (ENERG_LO[i] + ENERG_HI[i]) /2 );
		}
	}
	else
	{
		float dE = (E_chan_max - E_chan_min) / num_chans;
		for(int i=0; i<num_chans;i++)
		{
			ECHAN_LO.push_back( E_chan_min + (dE*i) );
			ECHAN_HI.push_back( E_chan_min + (dE*(i+1)) );
			ECHAN_MID.push_back( (ECHAN_LO[i] + ECHAN_HI[i]) /2 );
		}		
	}
	std::cout << "Response matrix has been reset to zeros." << "\n";
	make_empty_resp();
	set_N_CHAN();
}

void Response::set_N_GRP()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		N_GRP.push_back(1);
	}
}
void Response::set_F_CHAN()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		F_CHAN.push_back(0);
	}
}
void Response::set_N_CHAN()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		N_CHAN.push_back(num_chans);
	}
}

// If the shape of the response matrix is changed, the response matrix is reset to zeros.
void Response::make_empty_resp()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		std::vector<float> tmp_vec;
		for(int j=0; j<num_chans; j++)
		{
			tmp_vec.push_back(0.);
		}
		MATRIX.push_back(tmp_vec);
	}
}
// Make identity matrix
void Response::identity()
{
	for(int i=0; i<num_phot_bins; i++)
	{
		std::vector<float> tmp_vec;
		for(int j=0; j<num_chans; j++)
		{
			if(i==j)
			{
				tmp_vec.push_back(1);
			}
		}
		MATRIX.push_back(tmp_vec);
	}
}

// Decrease as 1/DeltaE^alpha from E_true
void Response::overDeltaE(float alpha)
{
	float col_sum; // Used for normalizing the column
	float tmp_val;
	for(int i=0; i<num_phot_bins; i++)
	{
		col_sum = 0;
		std::vector<float> tmp_vec;
		for(int j=0; j<num_chans; j++)
		{
			tmp_val = 1/(1+pow(abs( ECHAN_MID[j] - ENERG_MID[i]),alpha));
			tmp_vec.push_back(tmp_val);
			col_sum += tmp_val; 
		}
		// Normalize the 
		for(int j=0; j<num_chans; j++)
		{
			tmp_vec[j] /= col_sum;
		}
		MATRIX.push_back(tmp_vec);
	}

}

// Load response matrix from file
int Response::load_rsp_from_file(char file_name[])
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
	fits_open_file(&fptr, file_name, READONLY, &status);

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
	float tmp_ECHAN_LO=0.;
	float tmp_ECHAN_HI=0.;

	// Column names :
	// char colname_CHAN[] = "CHANNEL";
	char colname_ECHAN_LO[] = "E_MIN";
	char colname_ECHAN_HI[] = "E_MAX";

	/* Read rows of channel data*/
	for (irow=1; irow<=number_chans; irow++) 
	{
		// Select ECHAN_LO column
		fits_get_colnum(fptr, CASEINSEN, colname_ECHAN_LO, &colnum, &status);
		// Read column data at row i
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_ECHAN_LO, &anynull, &status);
		// Store column data in appropriate vector
		ECHAN_LO.push_back(tmp_ECHAN_LO);

		// ECHAN_HI
		fits_get_colnum(fptr, CASEINSEN, colname_ECHAN_HI, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_ECHAN_HI, &anynull, &status);
		ECHAN_HI.push_back(tmp_ECHAN_HI);

		// Calculate the center of the bin
		ECHAN_MID.push_back( (tmp_ECHAN_LO+tmp_ECHAN_HI)/2 );
	}

	/* Read response matrix  */
	// Move to Response table 
	fits_movabs_hdu(fptr, response_ihdu, &hdutype, &status);

	// Get number of rows
	fits_get_num_rows(fptr, &num_rows, &status);

	// Temporary variables to store element values
	float tmp_ENERG_LO=0.;
	float tmp_ENERG_HI=0.;
	float tmp_N_GRP=0.;
	float tmp_F_CHAN=0.;
	float tmp_N_CHAN=0.;
	float tmp_MATRIX[number_chans];

	// Column names :
	char colname_ENERG_LO[] = "ENERG_LO";
	char colname_ENERG_HI[] = "ENERG_HI";
	char colname_N_GRP[] = "N_GRP";
	char colname_F_CHAN[] = "F_CHAN";
	char colname_N_CHAN[] = "N_CHAN";
	char colname_MATRIX[] = "MATRIX";

	/* Read rows of spectrum data*/
	for (irow=1; irow<=num_rows; irow++) 
	{
		// Select ENERG_LO Column
		fits_get_colnum(fptr, CASEINSEN, colname_ENERG_LO, &colnum, &status);
		// Read column data at row i
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_ENERG_LO, &anynull, &status);
		// Store column data in appropriate vector
		ENERG_LO.push_back(tmp_ENERG_LO);

		// ENERG_HI
		fits_get_colnum(fptr, CASEINSEN, colname_ENERG_HI, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_ENERG_HI, &anynull, &status);
		ENERG_HI.push_back(tmp_ENERG_HI);

		// Calculate ENERG_MID
		ENERG_MID.push_back( (tmp_ENERG_LO+tmp_ENERG_HI)/2 );

		// N_GRP
		fits_get_colnum(fptr, CASEINSEN, colname_N_GRP, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_N_GRP, &anynull, &status);
		N_GRP.push_back(tmp_N_GRP);

		// F_CHAN
		fits_get_colnum(fptr, CASEINSEN, colname_F_CHAN, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_F_CHAN, &anynull, &status);
		F_CHAN.push_back(tmp_F_CHAN);

		// N_CHAN
		fits_get_colnum(fptr, CASEINSEN, colname_N_CHAN, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, 1, &nulldouble, &tmp_N_CHAN, &anynull, &status);
		N_CHAN.push_back(tmp_N_CHAN);

		// MATRIX
		fits_get_colnum(fptr, CASEINSEN, colname_MATRIX, &colnum, &status);
		fits_read_col(fptr, TFLOAT, colnum, irow, 1, number_chans, &nulldouble, &tmp_MATRIX, &anynull, &status);
		std::vector<float> tmp_vec(tmp_MATRIX, tmp_MATRIX + number_chans );
		MATRIX.push_back(tmp_vec);

	}

	E_phot_min = ENERG_LO.at(0);
	E_phot_max = ENERG_HI.back();
	num_phot_bins = num_rows;

	E_chan_min = ECHAN_LO.at(0);
	E_chan_max = ECHAN_HI.back();
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
