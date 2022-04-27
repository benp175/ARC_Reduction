#
#	arcticreduction.py
#
#	Pipeline for reducing arctic photometry. Specifically using quad amp with a fast readout.
#
#	Benjamin Proudfoot
#	11/17/2020
#

import astropy
from astropy.io import fits
from astropy import units as u
import numpy as np
from astropy.nddata import CCDData
import ccdproc
from astropy.modeling import models
#import matplotlib.pyplot as plt
import glob
from astropy.stats import mad_std
import os
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import Gaussian2DKernel
from datetime import datetime
import astrometry

# Input all files to reduce????
filename = "../Q4DD02/UT201115/lempo.0057.fits"

def inv_median(a):
	return 1 / np.median(a)

def overscan_remove(filename):
	raw = fits.open(filename)[0].data
	raw = CCDData(raw[1:2100,1:2048], unit=u.adu)
    
	# Do an overscan correction on the frame
	poly_model = models.Polynomial1D(3)
	oscan_cor = ccdproc.subtract_overscan(raw, fits_section = '[2047:2100,1:2048]', overscan_axis = 1, model = poly_model)
	oscan_cor = ccdproc.subtract_overscan(raw, fits_section = '[1:2,1:2048]', overscan_axis = 1, model = poly_model)
	trim = ccdproc.trim_image(oscan_cor, fits_section = '[3:2050,1:2048]')
	rawtrimmed = trim.data
	#rawtrimmed = CCDData(rawtrimmed, unit = u.adu)
    
	# Clean cosmic rays
	cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(rawtrimmed, sigclip=5, gain = 2.0)
	
	cr_cleaned = CCDData(cr_cleaned, unit = u.adu)

	return cr_cleaned, crmask


# Create master bias
def create_master_bias(path, key):
	# Read into a ccdproc image collection the files
	untrimmed_biases_files = glob.glob(path+key)
	trimmed_biases = []
	for file in untrimmed_biases_files:
		trimmed_biases.append(overscan_remove(file)[0])

	master_bias = ccdproc.combine(trimmed_biases, method = "median", sigma_clip=True, 
					sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
					sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
	return master_bias
    

# Create master dark
def create_master_dark(path, key, master_bias):
	# Read into a ccdproc image collection the files
	untrimmed_dark_files = glob.glob(path+key)
	trimmed_darks = []
	for file in untrimmed_dark_files:
		trimmed_darks.append(	ccdproc.subtract_bias(overscan_remove(file)[0], master_bias, add_keyword = False)	)

	master_dark = ccdproc.combine(trimmed_darks, method = "median", sigma_clip=True, 
					sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
					sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
	return master_dark
    


# Create master flat
def create_master_flat(path, key, master_bias, master_dark, darkexptime, sciencekey): #TODO: check median flat values 10000<x<50000s (all else fails==trim cornersin psf.py)
	untrimmed_flat_files = glob.glob(path+key)
	trimmed_flats = []
	flat_exp_time = []
	for file in untrimmed_flat_files:
		flat_hdu = fits.open(file)
		flat_filter = flat_hdu[0].header["FILTER"]
		img_filter = fits.open(sciencekey)[0].header["FILTER"]
		flat_data = flat_hdu[0].data
		mean_count = np.mean(flat_data)
		if flat_filter == img_filter:
			#print("Flat filter matches image filter")
			if mean_count < 10000:
				untrimmed_flat_files.remove(file)
				print("Undersaturated flat removed")
			elif mean_count > 60000:
				untrimmed_flat_files.remove(file)
				print("Oversaturated flat removed")
			else:
				trim = overscan_remove(file)[0]
				exptime = fits.open(file)[0].header["EXPTIME"]
				bias_corrected = ccdproc.subtract_bias(trim, master_bias, add_keyword = False)
				dark_corrected = ccdproc.subtract_dark(bias_corrected, master_dark, 
						       dark_exposure = darkexptime*u.second, data_exposure = exptime*u.second,
						       scale = True, add_keyword = False)
				trimmed_flats.append(dark_corrected)
				flat_exp_time.append(exptime)
		else:
			print("Flat removed: filter does not match image filter")
			untrimmed_flat_files.remove(file)

	master_flat = ccdproc.combine(trimmed_flats, method='median', scale=inv_median,
                                   sigma_clip=True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
                                   sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std)

	master_flat_mask = ccdproc.ccdmask(master_flat, findbadcolumns = True)
	return master_flat, master_flat_mask
    


# Apply cal frames
def reduce_data(path, biaskey, darkkey, flatkey, sciencekey, darkexptime):
	# Create directory for reduced data
	try:
		os.mkdir(path+"reduced_data")
	except FileExistsError:
		print("Warning: reduced_data directory already exists")
        
	unreduced_science_files = glob.glob(path+sciencekey)
    
	# Create master calibration frames
	master_bias = create_master_bias(path, biaskey)
	#master_bias.write(path+"reduced_data/master_bias.fits")
	print("Created master bias")

	master_dark = create_master_dark(path, darkkey, master_bias)
	#master_dark.write(path+"reduced_data/master_dark.fits")
	print("Created master dark")
    
	master_flat, mask = create_master_flat(path, flatkey, master_bias, master_dark, darkexptime, unreduced_science_files[0])
	#master_flat.write(path+"reduced_data/master_flat.fits")
	print("Created master flat")

	# Apply calibration frames!
	print("Applying calibration frames")

	for file in unreduced_science_files:
		trim, crmask = overscan_remove(file)
		imgheader = fits.open(file)[0].header
		nowstr = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
		imgheader.set("red_on", nowstr)
		exptime = imgheader["EXPTIME"] #TODO: copy over to other reduction script

		bias_corrected = ccdproc.subtract_bias(trim, master_bias, add_keyword = False)

        # ARC recommends only subtracting darks for exposures >~60s, so filter dark correction by exposure time.
		if exptime>=60:
			dark_corrected = ccdproc.subtract_dark(bias_corrected, master_dark, 
						       dark_exposure = darkexptime*u.second, data_exposure = exptime*u.second,
						       scale = True, add_keyword = False)
			del bias_corrected
			flat_corrected = ccdproc.flat_correct(dark_corrected, master_flat, add_keyword = False)
			del dark_corrected
		elif exptime<60:
			print("Skipped dark subtraction because of short exposure time.")
			flat_corrected = ccdproc.flat_correct(bias_corrected, master_flat, add_keyword = False)
			del bias_corrected
            
		print("Calibration frames applied")

		# Now fix any problems with bad pixels and cosmic rays using astropy's convolve function
		flat_corrected = flat_corrected.data
		flat_corrected[mask] = np.nan
		flat_corrected[crmask] = np.nan
		kernel = Gaussian2DKernel(x_stddev=9)
		final = interpolate_replace_nans(flat_corrected, kernel)
		reduced_image = CCDData(final, unit = u.adu)
		reduced_image.header = imgheader

		# Save the reduced image!
		filename = file.replace(path,"")
		reduced_image.write(path+"reduced_data/"+filename)
        
        #Platesolve the reduced image
		#astrometry.platesolve(path+"reduced_data/",filename)

		del reduced_image
