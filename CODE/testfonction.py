# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 19:15:52 2020

@author: plume_2.0
"""
import numpy as np 
#Import

import matplotlib.pyplot as plt 
#import scipy as sp
import scipy.constants as spc

def inbox(i,j,k,dimx,dimy,dimz):
    if i >= 0 and i < dimx and j >= 0 and j < dimy and k>=0 and k<dimz:
        in_box = True    
    else:
        in_box = False  
    return in_box

def frontbox(i,j,k,dimx,dimy,dimz):
    if i == 0 or i==(dimx-1) or j==0 or j==(dimy-1) or k == 0 or  k==(dimz-1):
        front_box = True
    else:
        front_box = False
    return front_box

# Fonction voisin qui te revoi les voisin d'un point
def voisin( i,j,k,dx,dy,dz):
    Vpos=np.zeros([8,3]).astype(int)
    for a in range(0,8):
        Vpos[a,0]=i+dx[a]
        Vpos[a,1]=j+dy[a]
        Vpos[a,2]=k+dz[a]
    return Vpos
#pourrai simplifier la lecture du code je vais pas reecrire tout les fonction avec mais 
# je le ferai eventuellement

#Fonction valide si le poin donner est dans le cristal     
def in_crystal(i,j,k,ice,dimx,dimy,dimz):
    i=int(i)#y a un bug il veul toujour des int ...
    j=int(j)
    k=int(k)
    if inbox(i,j,k,dimx,dimy,dimz):
        if (ice[i,j,k]==1)==True:
            in_ice=True 
           
        else:
            in_ice=False
            
    else:
        in_ice=False
    return in_ice

# fonction qui retourne la frontiere en transtion 
def frontiere(ice,dx,dy,dz,dimx,dimy,dimz):
    icepos=np.argwhere(ice==1)
    fronpos=[]
    for a in range(0,np.shape(icepos)[0]):
        Vpos=voisin(icepos[a][0],icepos[a][1],icepos[a][2],dx,dy,dz)
        for b in range(0,8):
            #print(infron(Vpos[b][0],Vpos[b][1],Vpos[b][2],fronpos))
            if  inbox(Vpos[b][0],Vpos[b][1],Vpos[b][2],dimx,dimy,dimz) and (in_crystal(Vpos[b][0],Vpos[b][1],Vpos[b][2],ice,dimx,dimy,dimz)==False) :
                fronpos.append(Vpos[b])
                
            if  len(fronpos)>0 and ((np.sum(np.sum((fronpos)==Vpos[b],axis=1)==3).any())>=2):
                fronpos.pop()

    
    
    return fronpos




# fonction in frontiere
def infron(i,j,k,fronpos):
    infr=[]
    for a in range(0,np.shape(fronpos)[0]):
        if fronpos[a][0]==i and  fronpos[a][1]==j and  fronpos[a][2]==k:
            infr.append(True)
        else:
            infr.append(False) 
    if sum(infr)>=1:
        return True
    else:
        return False



# Fonction qui calcule la constente de drain
def K(alpha,D,delta_dim,T,m):
    K=alpha*((2*delta_dim)/(np.sqrt(3)*D))*np.sqrt((spc.k*T)/(2*spc.pi*m))
    return K

# Fonction qui calcule nu_kin
def nukin(ratio_c,T,m):
    nu_kin=ratio_c*np.sqrt((spc.k*T)/(2*spc.pi*m))
    return nu_kin
    
    
#fonction qui te donne le bon coef de alpha en fonction des voisin
def valapha(i,j,k,dx,dy,dz,dimx,dimy,dimz,ice):
    Vpos=voisin(i,j,k,dx,dy,dz)
    #print()
    H=0 #nombre de voisin de glace horizotal initialisation
    V=0 #nombre de voisin de glace vertical initialisation
    for a in range(0,6):
        #print(in_crystal(Vpos[a][0],Vpos[a][1],Vpos[a][2],ice,dimx,dimy,dimz))
        if in_crystal(Vpos[a][0],Vpos[a][1],Vpos[a][2],ice,dimx,dimy,dimz)==True:
             H=H+1
    for b in range(6,8):
        #print(in_crystal(Vpos[b][0],Vpos[b][1],Vpos[b][2],ice,dimx,dimy,dimz))
        if in_crystal(Vpos[b][0],Vpos[b][1],Vpos[b][2],ice,dimx,dimy,dimz)==True:
             V=V+1
    #On donne la valeur de alpha (on set tu les val dansla fonction ou dans le code je sais pas)
    if  V==1 and H==0:
        alpha=0.01 # valeur pour vertical
    elif V<=1 and H==1:
        alpha=0.5 # valeur pour horizontal pour les autre  (pas sur revoir article)
    elif H==2 and V==0:
        alpha=0.1
    elif H==0 and V==0:
        alpha=0
    else:
        alpha=1
    #print("H")
    #print(H)
    #print("V")
    #print(V)
    return alpha 
#Note peu probablement etre amélioré genr just compter des vrai et des fau ?



#fonction qui calcule deltaV
def deltaV(alpha,Cvap,nu_kin,delta_dim,delta_t):
    acc=1
    dV=acc*alpha*nu_kin*Cvap*(delta_dim**2)*delta_t
    return dV


#fonction de la vitesse de surface normal

def nun(alpha,Cvap,nu_kin):
    #print("Cvap")
    #print(Cvap)
    #print("nu_kin")
    #print(nu_kin)
    #print("alpha")
    #print(alpha)
   
    nu_n=alpha*Cvap*nu_kin
    return nu_n

#fonction  de la hauteur croissance 
def delta_L(nu_n,temp):
     #pas sur i des probleme avec Cvap des foi je met just pour voir ce que sa 
    #donne passen le probleme des vap peur null pour plus tard( pou nu_n)
    dL=nu_n*temp
    if nu_n==0: #pas sur brian mais sa deverai marcher
        #print("temp")
        #print(temp)
        #print("v")
        #print(nu_n)
        dL=1
    return dL

# fonction qui calule le temp croissance pour une cellule se remplise
def tcroi(longeur,nu_n): #longeur qui reste a croite
    #pas sur i des probleme avec Cvap des foi je met just pour voir ce que sa 
    #donne passen le probleme des vap peur null pour plus tard( pou nu_n)
    if nu_n==0:
        tc=0
    else:
        tc=(longeur)/nu_n
    return tc
#nan


#Fonction pour croissance
###############################################################################

