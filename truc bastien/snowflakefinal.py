# -*- coding: utf-8 -*-
"""
Created on Tue May 28 22:26:40 2019

@author: Paulina
"""

import numpy as np
import matplotlib.pyplot as plt
#import random 

maxIt = 250 #maximal iteration of growth steps
mx = 251 #width 
my = 251 #and height of the snowflake arrays

rho = 0.7 #water vapor
gamma = 0.0000005 #constant influence from the third dim. which can be understood 
                #similar as the vapor,
                #but from the normal dimension of the snowflake
alpha = 1.2 #diffusion constant

#implement arrays for the ice mass value, receptive and
#non-receptive sites and set all elements equal to the inital
#humidity value rho
ice = np.empty((mx,my)); ice[:] = rho #values of ice mass
rec = np.zeros((mx,my)) #receptive site
nonrec = np.zeros((mx,my)) ; nonrec[:] = rho #non-receptive site

#initial condition: set an ice seed at the origin
x0 = int((mx-1)/2) ; y0 = int((my-1)/2)
ice[x0][y0] = 1.0

#initialize the vectors for the neighbors
dx = [-1, 0, -1, 1, 0, 1]
dy = [-1, -1, 0, 0, 1, 1]

#create function which checks if the considered site is inside the 
#snowflake or not. use therefore a boolean function
def inImage(i,j,mx,my):
    inside = False
    if i >= 0 and i < mx and j >= 0 and j < my:
        inside = True
    return inside

#define boolean funuciton which checks which sites
#are receptive and which not
def receptive(i,j,dx,dy):
    recep = False
    if ice[i,j] >= 1.: #if the value of the site in ice is greater equal 1, 
                        #its part of the snowflake
        recep = True
    else: #if the value of one neighbour site is greater equal 1, 
            #it could be part of the snowflake
        for k in range(6): # check the 6 neighbours
            kx = i + dx[k] ; ky = j + dy[k]
            if inImage(kx,ky,mx,my):
                if ice[kx,ky] >= 1.: #check value of k-th neighbour site
                    recep = True
                    break #if one neighbour is receptive, break, bc then site
                            #is receptive as well
    return recep

#implement the average rule of the nonrec sites it simulates the
#diffusion if the water in the air according to the diffusion equ.
def average(i,j,dx,dy,alpha):
    sum_neigh_nonrec = 0. # sum of value of all neighbours of considered site
    for k in range(6):#check neighbors
        kx = dx[k] + i ; ky = dy[k] + j
        if inImage(kx,ky,mx,my): #if neighbour is inImage, add value to sum
            sum_neigh_nonrec += nonrec[kx,ky]
    summy = nonrec[i,j] + (alpha/12.) * \
            (sum_neigh_nonrec - 6.*nonrec[i,j]) #implement solution of the
                                                #diffusion equation
    return summy

#Now we apply the functions to create a snowflake
#depending on if the considered site is part of the snowflake 
#the code handles the sites as receptive or non-receptive with if conditions
#for loops over all the sites (i and j) and the number of iteration
  
for g in range(maxIt):
    #print("step " + str(g+1) + " out of " + str(maxIt)) 

    for i in range(mx):
        for j in range(my):
            #check if site is part of the snowflake or not
            RECEPT_ij = receptive(i,j,dx,dy)
            if RECEPT_ij: #if part of the snowflake
                rec[i,j] = ice[i,j] + gamma #add influence from 3rd dim.
                nonrec[i,j] = 0. #set the value in nonrec zero,
                                #bc it's part of the snowflake   
            else: #if not part of snowflake, it's non-receptive:
                nonrec[i,j] = ice[i,j] #set value of ice equal nonrec
                     
    #calculate the new value of the nonrec sites due to diffusion            
    for i in range(mx):
         for j in range(my):            
                new_hum_value = average(i,j,dx,dy,alpha) #value after diffusion
                ice[i,j] = rec[i,j] + new_hum_value #overwrite old ice value

#prepare plot
#therefore we have to rearrange the site distance in the plot due to the
#hexagonal lattice

#array which contains at the end hex. lattice points
hexagonal = np.empty((2,mx,my))
#hexagonal basis vectors
a1 = np.array([1., -1./np.sqrt(3.)])
a2 = np.array([1., 1./np.sqrt(3.)])
#create the position of the lattice points due to the basis vectors
for i in range(mx):
    for j in range(my):
        hexagonal[:,i,j] = int(i)*a1 + int(j)*a2
        
#create an array which contains only values which are part of the snowflake 
#e.g., if the value is greater equal 1
#the np.where function saves the indeces of the corresponding elements in 
#two tupels
index_snow = np.where(ice >= 1.)       
#snow is array with all thevalues of the sites which are in the snowflake
snow = hexagonal[:,index_snow[0],index_snow[1]]
fig = plt.figure()
plt.scatter(snow[0,:], snow[1,:], c =ice[index_snow[0],\
                 index_snow[1]]-1., s = 1., \
                    cmap = 'cubehelix', marker = 'H')
plt.axis("off")
plt.show()
fig.savefig('rho_' + str(rho)+ '_gamma_' + str(gamma) + \
            '_alp_' + str(alpha) + '_iter_' + str(maxIt)+\
            '_mxmy_' + str(mx) + '.png')