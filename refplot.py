#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 16:29:17 2018

@author: steffensenchyna
"""

import aplpy
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib 
import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube as sc
from reproject import reproject_interp
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table

#Get Headers from Fits file
hdu_list = fits.open('M33-ARM05_yclean_test.tc_final.image.pbcor.K_props.fits', memmap=True)
hdu1 = fits.open(get_pkg_data_filename('mom0.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('mom0a.fits'))[0]
mol_data = Table(hdu_list[1].data)

#Reading the Mom 0 files and reference frame
ref = fits.getdata('m33_noema.armdistance.fits', ext=0)
mol = fits.getdata('mom0.fits', ext=0)
atomic = fits.getdata('mom0a.fits', ext=0)

dec = mol_data['YCTR_DEG']
ra = mol_data['XCTR_DEG']
coor = np.column_stack((ra,dec))

#Reprojecting the atomic data to match the molecular data 
array, footprint = reproject_interp(hdu2, hdu1.header)

#Taking the world coordinates from the cloud catalouge and converting the pixel coordinates
w = WCS('mom0.fits')
pixradec = w.wcs_world2pix(coor, 1)
pixradec = np.around(pixradec)
pixra = pixradec[:,0]
pixra = pixra.astype(np.int64)
pixde = pixradec[:,1]
pixde = pixde.astype(np.int64)
refoffset = ref[pixde,pixra]

#Extracting data from cloud catalouge 
vm = mol_data['MVIR_MSUN']
lm = mol_data['MLUM_MSUN']
vol = vm/lm
vel = mol_data['VCTR_KMS']
vlsr = vel+125  #Velocity from local standard of rest
mass = lm

plt.figure(1)
plt.plot(ref,mol,'bo',markersize=0.5)
plt.xlabel('Offset from Spiral Arm')
plt.ylabel('Intensity (K)')
plt.title('Molecular Brightness Distribution in the Spiral Arm')
plt.ylim(bottom=0)
plt.savefig('moloffset.png',bbox_inches='tight')

plt.figure(2)
plt.plot(refoffset,mass,'o')
plt.xlabel('Offset from Spiral Arm')
plt.ylabel('Mass (M$\odot$)')
plt.title('Molecular Cloud Mass Distribution in the Spiral Arm')
plt.savefig('massoffset.png',bbox_inches='tight')

plt.figure(3)
plt.plot(refoffset,vol,'o')
plt.xlabel('Offset from Spiral Arm')
plt.ylabel('$M_{vir}$/$M_{lum}$')
plt.title('$M_{vir}$/$M_{lum}$ vs offset')
plt.savefig('vloffset.png',bbox_inches='tight')

plt.figure(4)
plt.plot(refoffset,vlsr,'o')
plt.xlabel('Offset from Spiral Arm')
plt.ylabel('Velocity (km/s)')
plt.title('Velocity Distribution in the Spiral Arm')
plt.savefig('veloffset.png',bbox_inches='tight')

plt.figure(5)
plt.plot(ref,array,'bo',markersize=0.5)
plt.xlabel('Offset from Spiral Arm')
plt.ylabel('Intensity (K)')
plt.title('Atomic Brightness Distribution in the Spiral Arm')
plt.ylim(bottom=0)
plt.savefig('atmoffset.png',bbox_inches='tight')




