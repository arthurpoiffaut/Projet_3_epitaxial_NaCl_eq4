# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:27:13 2020

@author: plume_2.0
"""
import numpy as np
import matplotlib.pyplot as plt

cristal=np.zeros((100,100))
plt.figure(figsize=(6,3))

for b in range(0,10):
    for a in range(np.size(cristal,1)):
        c=a+b
        if (c%2)==0:
            cristal[99-b,99-a]=1
        else:
            cristal[99-b,99-a]=2




plt.figure(figsize=(20,15))
plt.imshow(cristal)




 



