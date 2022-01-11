# Platesolve reduced images using astrometry.net
# Seneca Heilesen
# 9/16/2021

import astropy
import astroquery
import photutils
import glob
import os
from astropy import wcs
from astropy.io import fits
from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet

# Set up astrometry.net
ast = AstrometryNet()

# Convert RA & dec to degrees
def radec_convert(sciencekey):
	hdu = fits.open(sciencekey)[0]
	ra_hms = hdu.header['RA']
	dec_dms = hdu.header['DEC']
	ra_hms = ra_hms.split(':',3)
	dec_dms = dec_dms.split(':',3)
	#print(ra_hms)
	#print(dec_dms)
	ra_deg = (float(ra_hms[0])+float(ra_hms[1])/60+float(ra_hms[2])/3600)*15
	dec_deg = float(dec_dms[0])+float(dec_dms[1])/60+float(dec_dms[2])/3600
	#print(ra_deg)
	#print(dec_deg)
	return ra_deg,dec_deg

# Platesolve reduced images
def platesolve(path, sciencekey, ra_dec_units=('degree','degree')):
    # Create directory for platesolved data
	try:
		os.mkdir(path+"platesolved_data")
	except FileExistsError:
		print("Warning: platesolved_data directory already exists")

    # Make list of files to platesolve
	files_to_solve = glob.glob(path+sciencekey)
    
    # Upload images to astrometry.net platesolving service
	for file in files_to_solve:
		hdu = fits.open(file)
		radeg,decdeg = radec_convert(file)
		
		wcs_header = ast.solve_from_image(file, center_ra=radeg, center_dec=decdeg, radius=5) # Solve from source list
        
        # Modify the header to include the new RA & dec
		print("original RA: "+hdu[0].header['RA'])
		print("original dec: "+hdu[0].header['DEC'])
		new_ra = wcs_header['CRVAL1']
		new_dec = wcs_header['CRVAL2']
		hdu[0].header['RA'] = new_ra
		hdu[0].header['DEC'] = new_dec
		print("new RA: "+str(hdu[0].header['RA']))
		print("new dec: "+str(hdu[0].header['DEC']))
        
        # Save image with new header
		filename = file.replace(path,"")
		#hdu[0].writeto(path+"platesolved_data/"+filename)
