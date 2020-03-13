# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:04:51 2020

@author: plume_2.0
"""

import numpy as np

#init des vecteur pour regarder les plus proche voisin
dx = [-1,0,-1,1,0,1,0,0]
dy = [-1,-1,0,0,1,1,0,0]
dz = [0,0,0,0,0,0,1,-1]
#FONCTION
###############################################################################

#fonction qui valide si le poin est dans la boit du croisence 

def inbox(i,j,k,dimx,dimy,dimz):
    if i >= 0 and i < dimx and j >= 0 and j < dimy and k>=0 and k<dimz:
        in_box = True    
    else:
        in_box = False  
    return in_box



#def diffusion(vap):
    
#fonction qui convertira la vapeur en glace
#def conversion(ice,fron,vap):
    

#fonction qui update ce qui fais partie du cristal
# def update_ice(fron,ice):   

    
    
#valide si le poin donner est dans le cristal     
def in_crystal(i,j,k,ice):
    if ice[i,j,k]==1:
        in_ice=True 
    else:
        in_ice=False
        
    return in_ice

# fonction qui retourne la frontiere en transtion 
def frontiere(ice,dx,dy,dz):
    icepos=np.argwhere(ice==1)
    fronpos=[]
    for a in range(0,np.shape(icepos)[0]):
        for b in range(0,8):
            pos=[icepos[a]+[dx[b],dy[b],dz[b]]] 
            #print(pos)
            if in_crystal(pos[0],pos[1],pos[2],ice) :
                pass
            else:
                fronpos.append(pos)
            
            #if in_crystal(pos[0],pos[1],pos[2],ice) :
                #rien
              #  pass
            #if np.shape(fronpos)[0]==1 or np.shape(fronpos)[0]==0 or np.shape(fronpos)[0]==2:
             #    fronpos.append(pos)  
            #elif np.shape(fronpos)[0]>1:
            #    if sum((np.sum(fronpos==pos,axis=0))==4)==1:
             #       print(sum((np.sum(fronpos==pos,axis=0))==4))
              #      pass
               # else:
                   
            #else:
              #  fronpos.append(pos)          #fauderai en enlever les doublon et les element qui son dans le cristall           
    #fronpos=np.delete(fronpos,0,axis=0)
    
    fronpos_clean= [] 
    for c in range(0,np.shape(fronpos)[0]):
        e=True
        for d in range(0,np.shape(fronpos)[0]):
            if np.array_equal(fronpos[c],fronpos[d]): 
                pass
            elif e:
                fronpos_clean.append(fronpos[c]) 
                e=False
    
    return fronpos_clean
#np.append(fronpos,icepos[a]+np.array([dx[b],dy[b],dz[b]]))

# fonction in frontiere
def infron(i,j,k,fronpos):
    
    for a in fronpos:
        if fronpos(a)[0]==i and  fronpos(a)[1]==j and  fronpos(a)[2]==k:
            infr=True
        else:
            infr=False  
        
        return infr


    
###############################################################################



#init variable 
dimx=16; #dimention x  ect
dimy=16;
dimz=16;
deltadim=1; # valeur de la dimention 1 micron
deltat=1; # valeur de la variation de temps
D=10**-5; #valeur du coe de diff   m^2/sece
Vcell=np.sqrt(3/2)*(deltadim**3);
nu_kin=133; # micro metre/sec
#voir poure les autre truc 


#autre para




#def mat 

# mat cristal 0 si pas dans le cristal
ice=np.zeros([dimx,dimy,dimz]);
# mat vapeur varie enfonction de la densiter deau en éta vapeur entre 0 et 1 ?
vap=np.zeros([dimx,dimy,dimz]);
#mat de la frontiere
fron=np.zeros([dimx,dimy,dimz]);











