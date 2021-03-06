#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:23:10 2018

@author: ssenchyna
"""

import aplpy
from astropy import units as u
import matplotlib.pyplot as plt
import pylab
import matplotlib 
import numpy as np

from spectral_cube import SpectralCube as sc
sc.allow_huge_operations=True
s = sc.read('M33-ARM05_yclean_test.tc_final.image.pbcor.K.fits')
mom0 = s.moment0()
mom0.write('mom0.fits', overwrite=True)
sbcube0 = s.spectral_slab(-169500 * u.m / u.s, -159500 * u.m / u.s)
sbcube0.write('sbcube0.fits', overwrite=True)
sbcube1 = s.spectral_slab(-159500 * u.m / u.s, -149500 * u.m / u.s)
sbcube1.write('sbcube1.fits', overwrite=True)
sbcube2 = s.spectral_slab(-149500 * u.m / u.s, -139500 * u.m / u.s)
sbcube2.write('sbcube2.fits', overwrite=True)
sbcube3 = s.spectral_slab(-139500 * u.m / u.s, -129500 * u.m / u.s)
sbcube3.write('sbcube3.fits', overwrite=True)
sbcube4 = s.spectral_slab(-129500 * u.m / u.s, -119500 * u.m / u.s)
sbcube4.write('sbcube4.fits', overwrite=True)
sbcube5 = s.spectral_slab(-119500 * u.m / u.s, -109500 * u.m / u.s)
sbcube5.write('sbcube5.fits', overwrite=True)
sbcube6 = s.spectral_slab(-99500 * u.m / u.s, -89500 * u.m / u.s)
sbcube6.write('sbcube6.fits', overwrite=True)
sbcube7 = s.spectral_slab(-89500 * u.m / u.s, -80500 * u.m / u.s)
sbcube7.write('sbcube7.fits', overwrite=True)


map = s.spatial_coordinate_map
ra = map[0][700//2][1024//2]
dec = map[1][700//2][1024//2]

cmap = matplotlib.cm.get_cmap('bwr')
norm = matplotlib.colors.Normalize(vmin=-150, vmax=-110)
thiscolor1 = cmap(norm(-145))
thiscolor2 = cmap(norm(-135))
thiscolor3 = cmap(norm(-125))
thiscolor4 = cmap(norm(-115))




fig = aplpy.FITSFigure('F475W_Arm_subim.fits')
fig.show_grayscale(vmin=0.1,vmax=1,stretch='log')
#fig.show_contour('sbcube0.fits',dimensions=[0,1],slices=[2],levels=[1, 1.5],linewidths=0.7,colors='#ff1493',layer='-169.5 km/s to -159.5 km/s') #pink              
#fig.show_contour('sbcube1.fits',dimensions=[0,1],slices=[2],levels=[1, 1.5],linewidths=0.7,colors='#ff0000',layer='-159.5 km/s to -149.5 km/s') #red 
fig.show_contour('sbcube2.fits',dimensions=[0,1],slices=[2],levels=2**np.linspace(0, 3,7),linewidths=0.7,colors=[thiscolor1],layer='-149.5 km/s to -139.5 km/s') #light red
fig.show_contour('sbcube3.fits',dimensions=[0,1],slices=[2],levels=2**np.linspace(0, 3,7),linewidths=0.7,colors=[thiscolor2],layer='-139.5 km/s to -129.5 km/s') #orange
fig.show_contour('sbcube4.fits',dimensions=[0,1],slices=[2],levels=2**np.linspace(0, 3,7),linewidths=0.7,colors=[thiscolor3],layer='-129.5 km/s to -119.5 km/s') #yellow
fig.show_contour('sbcube5.fits',dimensions=[0,1],slices=[2],levels=2**np.linspace(0, 3,7),linewidths=0.7,colors=[thiscolor4],layer='-119.5 km/s to -109.5 km/s') #green
#fig.show_contour('sbcube6.fits',dimensions=[0,1],slices=[2],levels=[1, 1.5],linewidths=0.7,colors='#87cefa',layer='-99.5 km/s to -89.5 km/s') #light blue
#fig.show_contour('sbcube7.fits',dimensions=[0,1],slices=[2],levels=[1, 1.5],linewidths=0.7,colors='#0000ff',layer='-89.5 km/s to -80.5 km/s') #blue
proxy = [plt.Line2D([],[],
					color=fig._layers[key].get_cmap().colors[0])
		 for key in fig._layers]
labels = fig._layers.keys()
fig._ax1.legend(proxy, labels, frameon=True, loc='upper left')                
       
                
fig.recenter(23.3882958333, 30.535200000000003, width=0.05, height=0.04)
fig.save('gasdatacmap1.pdf')





