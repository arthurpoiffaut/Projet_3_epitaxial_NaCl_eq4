# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:04:51 2020

@author: plume_2.0
"""

import numpy as np
import matplotlib.pyplot as plt 

#FONCTION
###############################################################################

#fonction qui valide si le poin est dans la boit du croisence 

def inbox(i,j,k,dimx,dimy,dimz):
    if i >= 0 and i < dimx and j >= 0 and j < dimy and k>=0 and k<dimz:
        in_box = True    
    else:
        in_box = False  
    return in_box



def diffusion(vap,ice,delta_tau,dx,dy,dz,dimx,dimy,dimz):
     vapout=vap
     #i=1
     for a in range(0,dimx):
         for b in range(0,dimy):
             for c in range(0,dimz):
                 sum1=0
                 sum2=0
                 #print(sum1)
                 #print(sum2)
                 if not in_crystal(a,b,c,ice,dimx,dimy,dimz):
                     for d in range(0,6):
                         #print(vap[[a]+dx[d],b+dy[d],c+dz[d]])
                         #print(in_crystal(a+dx[d],b+dy[d],c+dz[d],ice,dimx,dimy,dimz))
                         if in_crystal(a+dx[d],b+dy[d],c+dz[d],ice,dimx,dimy,dimz) == False:
                             sum1=(sum1+vap[a,b,c])
                         elif  (inbox(a+dx[d],b+dy[d],c+dz[d],dimx,dimy,dimz))==False:
                              sum1=(sum1+vap[a,b,c])
                         else:
                             sum1=(sum1+vap[a+dx[d],b+dy[d],c+dz[d]])
                             print(sum1)
                     for e in range(6,8):
                        #print(c+dz[e])#(in_crystal(a+dx[e],b+dy[e],c+dz[e],ice,dimx,dimy,dimz))
                        if  in_crystal(a+dx[e],b+dy[e],c+dz[e],ice,dimx,dimy,dimz):
                            #print(vap[a+dx[e],b+dy[e],c+dz[e]])
                            #pas bon doit etre la condition haaaaaaaaaaaa
                            sum2=sum2+vap[a,b,c]   
                        elif  (inbox(a+dx[e],b+dy[e],c+dz[e],dimx,dimy,dimz))==False:
                            sum2=sum2+vap[a,b,c]
                            #i=i+1
                            #print(i)
                        else: 
                            sum2=sum2+vap[a+dx[e],b+dy[e],c+dz[e]]
                            #i=i+1
                            #print(i)
                     print(sum1)
                     print(sum2)
                     vapout[a,b,c]=(2/3)*(delta_tau)*sum1+(delta_tau)*sum2+(1-6*delta_tau)*vap[a,b,c]
    
     return vapout


#fonction qui convertira la vapeur en glace
#def conversion(ice,fron,vap):
    

#fonction qui update ce qui fais partie du cristal
# def update_ice(fron,ice):   

    
    
#valide si le poin donner est dans le cristal     
def in_crystal(i,j,k,ice,dimx,dimy,dimz):
    if inbox(i,j,k,dimx,dimy,dimz):
        if ice[i,j,k]==1:
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
        for b in range(0,8):
            #print(icepos[a])
            #print(dx[b])
            #print(dy[b])
            #print(dz[b])
            pos=[icepos[a]+[dx[b],dy[b],dz[b]]] 
            #print(pos)
            #print(pos[0][0])
            #print(pos[0][1])
            #print(pos[0][2])
            #print(inbox(pos[0][0],pos[0][1],pos[0][2],dimx,dimy,dimz))
            if  ((inbox(pos[0][0],pos[0][1],pos[0][2],dimx,dimy,dimz)) == False):
                pass
            elif in_crystal(pos[0][0],pos[0][1],pos[0][2],ice,dimx,dimy,dimz) :
                pass
            else:
                #p=((inbox(pos[0][0],pos[0][1],pos[0][2],dimx,dimy,dimz)) == False)
                #print(p)
                fronpos.append(pos[0])
    #peut etre retourner une matrice a la place ?
    return fronpos





#np.append(fronpos,icepos[a]+np.array([dx[b],dy[b],dz[b]]))

#update frontiere

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
dimx=10; #dimention x  ect
dimy=10;
dimz=10;

#init des vecteur pour regarder les plus proche voisin
dz = [-1,0,-1,1,0,1,0,0]
dy = [-1,-1,0,0,1,1,0,0]
dx = [0,0,0,0,0,0,1,-1]


#val
delta_dim=10**(-6); # valeur de la dimention 1 micron
delta_t=10**(-9); # valeur de la variation de temps
D=10**-5; #valeur du coe de diff   m^2/sece
Vcell=np.sqrt(3/2)*(delta_dim**3);
nu_kin=133; # micro metre/sec
delta_tau=(D*delta_t)/(delta_dim)**2



#voir poure les autre truc 


#autre para




#def mat 


# je sais pas pourquoi mais sa chie avec l ordre dim x dim y dim z sa pas dans 
# le bonne ordre du moin pour les test en 2d je sais pas pourquoi mais 
# fauderai recrire tout les fonction dans cet ordre different je croi je suis pas sur 
# mais sa chage rien temp qu on symetrique pour le moment mais quand plot la sum devera etre chenge i guess


# mat cristal 0 si pas dans le cristal
ice=np.zeros([dimx,dimy,dimz])
# mat vapeur varie enfonction de la densiter deau en éta vapeur entre 0 et 1 ?
vap=np.zeros([dimx,dimy,dimz])
#mat de la frontiere
fron=np.zeros([dimx,dimy,dimz])

#test glace
ice[5,5,5]=1
ice[1,1,1]=1

vap[0,0,0]=1
#x et z invercer???
#comme on va fair les image pour le moment 
#plt.imshow(np.sum(ice,axis=0),interpolation='spline16', cmap='viridis')




