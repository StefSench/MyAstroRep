#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 01:31:23 2018

@author: senchyna
"""
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-2*np.pi, 2*np.pi, 0.01)

def sinewave(t):
    return 2*np.sin(x-(11/60))*np.cos(t-(1/60))
t1=0
t2=0.5
t3=1.0 
t4=1.5
t5=2.0
t6=2.5
t7=3.0  
t8=3.5 
wave=sinewave(t1)

plt.plot(x,sinewave(t1),label='t=0s')
plt.plot(x,sinewave(t2),label='t=0.5s')
plt.plot(x,sinewave(t3),label='t=1s')
plt.plot(x,sinewave(t4),label='t=1.5s')
plt.plot(x,sinewave(t5),label='t=2s')
plt.plot(x,sinewave(t6),label='t=2.5s')
plt.plot(x,sinewave(t7),label='t=3.0s')
plt.plot(x,sinewave(t8),label='t=3.5s')
plt.grid()
plt.title('Potting Wave Function for varying t')
plt.xlabel('Radians')
plt.ylabel("$\psi$")
plt.legend()

plt.xlim(xmin=-2*np.pi)
plt.xlim(xmax=2*np.pi)
plt.savefig('waveplot.png')



