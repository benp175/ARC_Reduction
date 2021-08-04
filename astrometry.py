# Platesolve reduced images using astrometry.net
# Seneca Heilesen
# 8/3/2021

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

# Platesolve reduced images
def platesolve(path, sciencekey):
    # Create directory for platesolved data
	try:
		os.mkdir(path+"platesolved_data")
	except FileExistsError:
		print("Warning: platesolved_data directory already exists")

    # Make list of files to platesolve
	files_to_solve = glob.glob(path+sciencekey)
    
    # Upload images to astrometry.net platesolving service
	for file in files_to_solve:
		hdul = fits.open(file)
		wcs_header = ast.solve_from_image(file) # Solve from source list
        
        # Modify the header to include the new RA & dec
		print("original RA: "+hdul[0].header['RA'])
		print("original dec: "+hdul[0].header['DEC'])
		new_ra = wcs_header['CRVAL1']
		new_dec = wcs_header['CRVAL2']
		hdul[0].header['RA'] = new_ra
		hdul[0].header['DEC'] = new_dec
		print("new RA: "+str(hdul[0].header['RA']))
		print("new dec: "+str(hdul[0].header['DEC']))
        
        # Save image with new header
		filename = file.replace(path,"")
		hdul[0].writeto(path+"platesolved_data/"+filename)