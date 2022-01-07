"""
Author: Michael Moss
Contact: mikejmoss3@gmail.com
Last Edited: 2021-07-21

This file contains methods to make simulated data into instrument ready data products. 
For instance, writing data into FITS file readable by XSPEC.

"""

import numpy as np
from astropy.io import fits
import sys

def write_fits_file(spec_file_name,outfile_name="outfile.fit"):
	"""
	Method to write input data into FITS file readable by XSPEC.
	"""

	# Load spectrum data
	spectrum_data = np.genfromtxt(spec_file_name,dtype=[("ENERGY",float),("RATE",float)])


	# Make Primary header
	hdu_primary = write_prim_header()

	# Make first header (contains data)
	hdu_data = write_data_header(col1_data=spectrum_data['ENERGY'],col2_data=spectrum_data['RATE'])

	hduL = fits.HDUList([hdu_primary,hdu_data])
	hduL.writeto(outfile_name)

def write_prim_header():

	hdu_prim = fits.PrimaryHDU()

	hdu_prim.header.set('SIMPLE', value='T',comment='file does conform to FITS standard')
	hdu_prim.header.set('BITPIX', value=8,comment='number of bits per data pixel')
	hdu_prim.header.set('NAXIS', value=0,comment='number of data axes')
	hdu_prim.header.set('EXTEND', value='T',comment='FITS dataset may contain extensions')
	
	return hdu_prim

def write_data_header(col1_data,col2_data,col1_label='ENERGY',col2_label='RATE'):

	# Write data into fits columns
	c1=fits.Column(name=col1_label,array=col1_data,format='J')
	c2=fits.Column(name=col2_label,array=col2_data,format='J')

	columns=fits.ColDefs([c1,c2])
	hdu=fits.BinTableHDU.from_columns(columns)

	# Write header key words
	hdu.header.set('XTENSION', value = 'BINTABLE', comment = 'binary table extension')
	hdu.header.set('BITPIX', value = 8, comment= '8-bit bytes')
	hdu.header.set('NAXIS', value = 2, comment= '2-dimensional binary table')
	hdu.header.set('NAXIS1', value = 12, comment= 'width of table in bytes')
	hdu.header.set('NAXIS2', value = len(col1_data), comment= 'number of rows in table')
	hdu.header.set('PCOUNT', value = 0, comment= 'size of special data area')
	hdu.header.set('GCOUNT', value = 1, comment= 'one data group (required keyword)')
	hdu.header.set('TFIELDS', value = 2, comment= 'number of fields in each row')
	hdu.header.set('TTYPE1', value = col1_label, comment = 'label for field   1')
	hdu.header.set('TFORM1', value = 'J', comment = 'data format of field: 4-byte INTEGER')
	hdu.header.set('TTYPE2', value = col2_label, comment = 'label for field   2')
	hdu.header.set('TFORM2', value = 'J', comment = 'data format of field: 4-byte INTEGER')
	hdu.header.set('EXTNAME', value = 'SPECTRUM', comment= 'name of this binary table extension')
	hdu.header.set('HDUCLASS', value = 'OGIP', comment= 'format conforms to OGIP standard')
	hdu.header.set('HDUCLAS1', value = 'SPECTRUM', comment= 'PHA dataset (OGIP memo OGIP-92-007)')
	hdu.header.set('HDUVERS1', value = '1.2.1', comment= 'Version of format (OGIP memo OGIP-92-007a)')
	hdu.header.set('TELESCOP', value = 'UNKNOWN ', comment= 'mission/satellite name')
	hdu.header.set('INSTRUME', value = 'UNKNOWN ', comment= 'instrument/detector name')
	hdu.header.set('CHANTYPE', value = 'PHA', comment= 'channel type (PHA, PI etc)')
	hdu.header.set('HISTORY', value = "Simulated data written in XSPEC form")
	hdu.header.set('RESPFILE', value = '        ', comment= 'associated redistrib matrix filename')
	hdu.header.set('ANCRFILE', value = '        ', comment= 'associated ancillary response filename')
	hdu.header.set('CORRFILE', value = '        ', comment= 'associated correction filename')
	hdu.header.set('CORRSCAL', value = -1., comment = 'correction file scaling factor')
	hdu.header.set('BACKFILE', value = '        ', comment= ' associated background filename')
	hdu.header.set('EXPOSURE', value = 1., comment = 'exposure (in seconds)')
	hdu.header.set('TLMIN1', value = 1, comment = 'Lowest legal channel number')
	hdu.header.set('TLMAX1', value = 5000, comment = 'Highest legal channel number')
	hdu.header.set('DETCHANS', value = 5000, comment = 'total number possible channels')
	hdu.header.set('POISSERR', value = 'T', comment = 'Pois. err assumed ?')
	hdu.header.set('AREASCAL', value = 1., comment = 'area scaling factor')
	hdu.header.set('BACKSCAL', value = 1., comment = 'background file scaling factor')

	return hdu


if __name__ == "__main__":
	write_fits_file(sys.argv[1])