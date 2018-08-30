#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 01:31:23 2018

@author: senchyna
"""
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-2*np.pi, 2*np.pi, 0.01)
l= float( input("Please input a value for lambda, "))
def MySineWave(l):
    return np.sin(((2*np.pi)/l)*x)
    

plt.plot(x,MySineWave(l))
plt.title('Potting Sin Curve for some value $\lambda$')
plt.xlabel('Radians')
plt.ylabel("Sin($\lambda$)")
plt.show()






