# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:38:56 2020

@author: Alexis
"""

import numpy as np
import matplotlib.pyplot as plt


dimx=100
dimy=100
dimz=13

cx=int((dimx+1)/2)
cy=int((dimy+1)/2)
cz=int((dimz+1)/2)

#Initialisation de la matrice test
mattest = np.zeros((dimx,dimy,dimz))

#Pour carré
#mattest[45:55,45:55]=1
#mattest[40:45,40:60]=1
#mattest[55:60,40:60]=1
#mattest[45:55,40:45]=1
#mattest[45:55,55:60]=1

#Pour hexagone
mattest[cx,cy,cz]=1
mattest[cx-1,cy,cz]=1
mattest[cx-1,cy-1,cz]=1
mattest[cx,cy-1,cz]=1
mattest[cx+1,cy,cz]=1
mattest[cx+1,cy+1,cz]=1
mattest[cx,cy+1,cz]=1



matrix=np.sum(mattest,axis=2) #image d'une matrice 2D seulement

plt.imshow(matrix,interpolation='spline16',cmap='viridis')


#

def hexplot(matflocon):#prend la matrice 2D contenant l'image du flocon et ''l'hexagonifie'' en scatter
    hexagonal = np.empty((2,dimx,dimy)) #array contenant les points hex
    a1 = np.array([1., -1./np.sqrt(3.)]) #vecteur de la base hexagonale
    a2 = np.array([1., 1./np.sqrt(3.)])
    for i in range(dimx):
        for j in range(dimy):
            hexagonal[:,i,j] = int(i)*a1 + int(j)*a2
    index_flocon = np.where(matflocon >= 1)  #0 si on veut tous les points de la matrice.
    flocon = hexagonal[:,index_flocon[0],index_flocon[1]]
    fig = plt.figure()
    plt.scatter(flocon[0,:], flocon[1,:], c =matflocon[index_flocon[0],\
                 index_flocon[1]]-1., s = 10., \
                    cmap = 'cubehelix', marker = 'H')
    plt.axis("off")
    fig.savefig('testfloconhex_'+'mx'+str(dimx)+'my'+str(dimy)+'.png') #Enregistrement de la figure au même endroit où le code .py est situé
    plt.show()


hexplot(matrix)