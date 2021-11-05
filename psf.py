#psf.py
#
#Code to solve ePSF of ARCTIC data
#
#Seneca Heilesen
#11/4/2021

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from photutils.detection import find_peaks
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.nddata import NDData
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder
from astropy.io import fits
from astropy.nddata import CCDData
from astropy import units as u
import glob

def epsf(path, sciencekey):
    
    # Create directory for ePSFs
    try:
        os.mkdir(path+"ePSF")
    except FileExistsError:
        print("Warning: ePSF directory already exists")
    
    #Make list of files to create PSFs for
    files_list = glob.glob(path+sciencekey)
    
    for file in files_list:
        #Load data
        hdu = fits.open(file)[0]
        data = hdu.data

        #Show the image
        #norm = simple_norm(data, 'sqrt', percent=99.)
        #plt.imshow(data, norm=norm, origin='lower', cmap='viridis')

        #Finds peaks of light->objects
        peaks_tbl = find_peaks(data, threshold=np.percentile(np.array(data),99.5), npeaks=25)
        peaks_tbl['peak_value'].info.format = '%.8g'  # for consistent table output  
        print(peaks_tbl)

        #Get stars in format for cutouts
        size = 25 #25x25 pixel image
        hsize = (size - 1) / 2
        x = peaks_tbl['x_peak']
        y = peaks_tbl['y_peak']
        
        #Exclude repeated peaks from saturated stars
        index = 0
        mask = []
        while index < len(x)-1:
            if abs(x[index] - x[index+1]) < 15:
                if abs(y[index] - y[index+1]) < 15:
                    mask.append(False)
            else:
                mask.append(True)
            index = index + 1
        index = 0
        mask.append(False)
        print(mask)
        print(len(mask)) # for testing
        print(len(x)) # for testing
    
        #Create table of good star positions
        stars_tbl = Table()
        stars_tbl['x'] = x[mask]  
        stars_tbl['y'] = y[mask]
        print(stars_tbl)
        mask = []
        
        #Subtract background (sigma_clipped median)
        mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)
        data = np.array(data)
        data = data - float(median_val)
        
        #Change object type of data to use in extract_stars
        nddata = NDData(data=data)
        
        #Create star cutouts
        stars = extract_stars(nddata, stars_tbl, size=25)
            
        #Show first 25 cutouts
        #nrows = 5
        #ncols = 5
        #fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
        #ax = ax.ravel()
        #for i in range(nrows*ncols):
        #    norm = simple_norm(stars[i], 'log', percent=99.)
        #    ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')
        #i = 1
        
        #Construct ePSF
        epsf_builder = EPSFBuilder(oversampling=1, maxiters=20, progress_bar=True)
        epsf, fitted_stars = epsf_builder(stars)
        epsffit = CCDData(epsf.data, unit=u.adu)
        
        #Save ePSF
        filename = file.replace(path,"")
        filename=filename.replace('.fits','.psf.fits')
        epsffit.write(path+"ePSF/"+filename)
        
        #Show constructed ePSF
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.figure()
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        
        #Save ePSF as png
        plt.savefig(path+'ePSF/'+filename.replace('.psf.fits','.psf.png'))