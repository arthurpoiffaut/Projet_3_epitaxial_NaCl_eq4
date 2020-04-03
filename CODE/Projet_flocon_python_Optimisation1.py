# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 13:30:44 2020

@author: plume_2.0
"""

#Import
import numpy as np
import matplotlib.pyplot as plt 
#import scipy as sp
import scipy.constants as spc
import time 
from joblib import Parallel, delayed # truc de multi proceesing


#FONCTION DE BASE
###############################################################################

#Fonction qui valide si le point est dans la boîte de croissance 

def inbox(pos,dim_box):
    if (pos >= 0).all() and (pos<dim_box).all():
        in_box = True    
    else:
        in_box = False  
    return in_box


#Fonction qui valide si on est dans la frontiere de la boite
def frontbox(pos,dim_box):
    if (pos==[0,0,0]).any() or (pos==(dim_box-1)).any() :
        front_box = True
    else:
        front_box = False
    return front_box


#Fonction voisin qui renvoit les voisins d'un point
def voisin( pos,mat_voi):
    Vpos=np.ones([8,3]).astype(int)
    Vpos=np.ones([8,3]).astype(int)*np.array(pos)
    Vpos=Vpos+mat_voi
    return Vpos

#Fonction valide si le point donné est dans le cristal 
def in_crystal(pos,ice,dim_box):
    if inbox(pos,dim_box):
        if (ice[pos[0],pos[1],pos[2]]==1)==True:
            in_ice=True 
           
        else:
            in_ice=False
            
    else:
        in_ice=False
    return in_ice


#Fonction qui retourne la frontiere en transition 
def frontiere(ice,mat_voi,dim_box): # se parallilse mallllllllllll fauderai la repencer mais cest du avec le pop pour en lever les doublon
        
    icepos=np.argwhere(ice==1)
    fronpos=np.zeros([1,3])
    
    for pos1 in icepos :
        Vpos=voisin(pos1,mat_voi)
        for pos2 in Vpos:
            if  inbox(pos2,dim_box) and (in_crystal(pos2,ice,dim_box)==False) :
                fronpos=np.concatenate((fronpos, [pos2]), axis=0)
                
            # elif  len(fronpos)>0 and ((np.sum(np.sum((fronpos)==pos2,axis=1)==3).any())>=2):
            #     fronpos.pop()
    
    fronpos=np.unique(fronpos,axis=0)
    fronpos=np.delete(fronpos,0,axis=0)

    
    
    return fronpos

#fonction qui valide si on est dans la frontiere
def infron(pos,fronpos):

    f=(fronpos==np.array([pos[0],pos[1],pos[2],np.nan]))
    # print(f)
    # print()
    p=((np.sum(f,axis=1))==3)
    # print(p)
    if p.any():
        pos2=np.argwhere(p)
        # print(pos2)
        return np.array([pos[0],pos[1],pos[2],pos2[0][0]])
    else:
        return np.nan


#Fonction qui calcule la constante de drain AVec les varible global ?
def K(alpha):
    K=alpha*((2*delta_dim)/(np.sqrt(3)*D))*np.sqrt((spc.k*T)/(2*spc.pi*m))
    return K

#fonction qui donne le bon coefficent de alpha en fonction des voisins
def valapha(pos,mat_voi,dim_box,ice):
    Vpos=voisin(pos,mat_voi)
    #print()
    H=0 #nombre de voisin de glace horizotal initialisation
    V=0 #nombre de voisin de glace vertical initialisation
    for i1 in range(0,6):#voisin horizontaux
        if in_crystal(Vpos[i1],ice,dim_box)==True:
             H=H+1
    for i2 in range(6,8): #voisin verticau
        if in_crystal(Vpos[i2],ice,dim_box)==True:
             V=V+1
    #On donne la valeur de alpha (on set tu les val dansla fonction ou dans le code je sais pas)
    if  V==1 and H==0:
        alpha=0 # valeur pour vertical
    elif V<=1 and H==1:
        alpha=1 # valeur pour horizatal et vertical petit 
    elif H==2 and V==0:
        alpha=1
    elif H==0 and V==0:
        alpha=0
    else:
        alpha=1

    return alpha 
#Note: peut probablement etre amélioré genre juste compter des vrais et des faux ?
#pas sur si besoin de vectoriser pas tend de gain a fair je croi 

#Fonction qui calcul la vitess de croissence normal a la surface
def nun(alpha,Cvap,nu_kin):
    #print("Cvap")
    #print(Cvap)
    #print("nu_kin")
    #print(nu_kin)
    #print("alpha")
    #print(alpha)
   
    nu_n=alpha*Cvap*nu_kin
    return nu_n


#Fonction de la hauteur de croissance 
def delta_L(nu_n,temp):

    dL=nu_n*temp
#    if nu_n==0: #pas sur brian mais ca devrait marcher
#        dL=0

    return dL

#fonction qui donne le temp de croissence
def tcroi(longeur,nu_n): #longueur restant à croître
    #pas sur i des problemes avec Cvap des fois je mets juste pour voir ce quecsa 
    #donne passen le probleme des vapeurs nulles pour plus tard( pou nu_n)
    if nu_n==0:
        tc=0
    else:
        tc=(longeur)/nu_n
    return tc



#FONCTION AVENCÉ
###############################################################################

#Fonction de relaxation 1, alpha independant de sigma 
#condition aux frotières de la boîte, sont une source constante
def relax1(vap,ice,fron_state,dim_box,mat_voi,sigma_limit,Pos):
    #pmultiprooooccees
    #def relax1_multi(pos):
    
    #peut metre les valeur local pour certin truc ici voir si sa aide?
    
    
    
    vap_out=np.zeros(np.shape(vap))
    
    for pos1 in Pos:
        
        sum1=0
        sum2=0
                
        if in_crystal(pos1,ice,dim_box)== False:
            if ((type(infron(pos1,fron_state)))== np.ndarray):   
#               print("True")
                Vpos=voisin(pos1,mat_voi)
                alpha=valapha(pos1,mat_voi,dim_box,ice)
                #print(alpha)
                Kval=K(alpha)
                
                for i1 in range(0,6):
                    if in_crystal(Vpos[i1],ice,dim_box) == True: #peut etre metre dans un for avec Vpos mais faus une conditon pour les 2 haut
                        sum1=sum1+vap[pos1[0],pos1[1],pos1[2]]
                        
                    elif  inbox(Vpos[i1],dim_box) == False:
                        sum1=sum1+vap[pos1[0],pos1[1],pos1[2]]
                    else:
                        sum1=sum1+vap[int(Vpos[i1][0]),int(Vpos[i1][1]),int(Vpos[i1][2])]
                for i2 in range(6,8):
                    if in_crystal(Vpos[i2],ice,dim_box) == True:
                        sum2=sum2+ +vap[pos1[0],pos1[1],pos1[2]]
                    elif inbox(Vpos[i2],dim_box) == False:
                                sum2=sum2+vap[pos1[0],pos1[1],pos1[2]]
                    else:
                        sum2=sum2+vap[int(Vpos[i2][0]),int(Vpos[i2][1]),int(Vpos[i2][2])]
                    
                    vap_out[pos1[0],pos1[1],pos1[2]]=((2/3)*sum1+sum2)/(Kval+6)
                            
            elif frontbox(pos1,dim_box) == True:
                vap_out[pos1[0],pos1[1],pos1[2]]=sigma_limit
            else:
                Vpos=voisin(pos1,mat_voi)
                for i3 in range(0,6):
                    if in_crystal(Vpos[i3],ice,dim_box) == True:
                            sum1=sum1+ vap[pos1[0],pos1[1],pos1[2]]
                    elif  inbox(Vpos[i3],dim_box)==False:
                                sum1=sum1+vap[pos1[0],pos1[1],pos1[2]]
                    else:
                            sum1=sum1+vap[int(Vpos[i3][0]),int(Vpos[i3][1]),int(Vpos[i3][2])] 
                            
                for i4 in range(6,8):
                    if in_crystal(Vpos[i4],ice,dim_box) == True:
                        sum2=sum2+ +vap[pos1[0],pos1[1],pos1[2]]
                    elif inbox(Vpos[i4],dim_box) == False:
                        sum2=sum2+vap[pos1[0],pos1[1],pos1[2]]
                    else:
                        sum2=sum2+vap[int(Vpos[i4][0]),int(Vpos[i4][1]),int(Vpos[i4][2])]
                        #vap_out[a,b,c]=(2/3)*(delta_tau)*sum1+(delta_tau)*sum2+(1-6*delta_tau)*vap[a,b,c]
                vap_out[pos1[0],pos1[1],pos1[2]]=((2/3)*sum1+sum2)/6
        else:
            vap_out[pos1[0],pos1[1],pos1[2]]=2#np.nan

    return vap_out
                    
#Fonction de relaxation 2, alpha dependant de sigma 
#Condition aux frontières de la boîte sont une source constante??









#Fonction qui trouve le temps minimal pour que l'une des cellules se remplisse
def tmin(fron_state,ice,vap,nu_kin,mat_voi,dim_box,Lcell):
    list_tc=[]
    for a in range(0,np.shape(fron_state)[0]):
        pos=[int(fron_state[a][0]),int(fron_state[a][1]),int(fron_state[a][2])]
        longeur0=fron_state[a][3]
        longeur=Lcell-longeur0
        Cvap=vap[pos[0],pos[1],pos[2]]
        alpha=valapha(pos,mat_voi,dim_box,ice)
        nu_n=nun(alpha,Cvap,nu_kin)
        tc=tcroi(longeur,nu_n)
        if tc!=0:    
            list_tc.append(tc)
        
        
    t_min=min(list_tc)
    return t_min
#dLcell

#fonction qui fais la croissence des frontiere
def croissance(ice,fron_state,vap,delta_dim,mat_voi,dim_box,nu_kin,t_min):
    fron_state_c=fron_state.copy()
    fron_state_ini=fron_state.copy()
    for i1 in range(0,np.shape(fron_state)[0]):
        pos=[int(fron_state[i1][0]),int(fron_state[i1][1]),int(fron_state[i1][2])]
        alpha=valapha(pos,mat_voi,dim_box,ice)

        Cvap=vap[pos[0],pos[1],pos[2]]

        nu_n=nun(alpha,Cvap,nu_kin)
        
       
        dL=delta_L(nu_n,t_min)

        fron_state_c[i1][3]=(fron_state_ini[i1][3]+dL)

    return fron_state_c
            
#fron_stat plus sur si les positions de dans c'est un peu con je suis pas sur
    


#Fonction qui update la matrice ice et la liste frontiere
def update_fron(fron_state,ice,mat_voi,dim_box,Lcell):
    fron_state_new=np.zeros([1,4])
    ice_new=ice

    for state1 in fron_state:
        
        if np.round(state1[3],8)>=Lcell:
            
            ice_new[int(state1[0]),int(state1[1]),int(state1[2])]=1
                        
    state_new=frontiere(ice_new,mat_voi,dim_box)

    for state2 in state_new:
        
        pos=infron(state2,fron_state)
        # print(np.array([np.concatenate((state2,np.array([0])))]))
        # print(fron_state_new)
        # print(pos)
        # print(np.concatenate((fron_state_new,state2), axis=0))
        if  (type(pos)== np.ndarray) :
            # print(np.array([[pos[0],pos[1],pos[2],fron_state[int(pos[3])][3]]]))
            fron_state_new=np.concatenate((fron_state_new,np.array([[pos[0],pos[1],pos[2],fron_state[int(pos[3])][3]]])), axis=0)
            #print([fronpos[b][0],fronpos[b][1],fronpos[b][2],fron_state[pos[3]][3]])
        else:
            state2=np.array([np.concatenate((state2,np.array([0])))])
            fron_state_new=np.concatenate((fron_state_new,state2), axis=0)
            #append([state2[0],state2[1],state2[2],0])
            
    fron_state_new=np.delete(fron_state_new,0,axis=0)

    return [fron_state_new,ice_new]

#np.concatenate((fronpos, [pos2]), axis=0)

#FONCTION IMAGE
###############################################################################





# INITIALISATION VARIABLE 
###############################################################################


#Initiation de la matrice pour regarder les plus proches voisins

mat_voi=np.array([[-1,-1,0],[0,-1,0],[-1,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[0,0,-1]])


# Constantes :


# valeur de la dimention (m)
delta_dim=50*10**(-6); 


# valeur de la variation de temps 
delta_t=10*10**(-12); 


#temperature en kelvin
T=258 

#masse d'une molecule d'eau
m=2.988*10**(-26)


#valeur du coe de diff   (m^2/s)
D=2*10**-5 

#ratio des c_sat et c_solide
ratio_c=1*10**(-6) #pas sur ici faudrait revoir c_sat/c_solide

#nu_kin=133; # micrometre/s en fonction aussi
#Fonction qui calcule nu_ki
nu_kin=ratio_c*np.sqrt((spc.k*T)/(2*spc.pi*m))

#concentration aux extremitées
sigma_limit=1 


#longueur d'une cellule 
Lcell=np.round((np.sqrt(3)/2)*(delta_dim),8)


#valeur du delta tau parametre sans dimention pour les step de temp
delta_tau=(D*delta_t)/(delta_dim)**2


#dimention dim_box[0]=dimx ect
dim_box=np.array([100,100,5]).astype(int) 

center=(np.floor((dim_box)/2)).astype(int)


#Definition de matrices 

# mat cristal 0 si pas dans le cristal
ice_ini=np.zeros(dim_box)
# mat vapeur varie en fonction de la densité d'eau en état vapeur entre 0 et 1 ?
vap=np.ones(dim_box)

ice_ini[center[0],center[1],center[2]]=1

#initalisation d un vecteur de parcour avec tout les position voulue
Pos=[]
for i1 in range(0,dim_box[0]):
    for i2 in range(0,dim_box[1]):
        for i3 in range(0,dim_box[2]):
            Pos.append(np.array([i1,i2,i3]))       






#initialisation de l eta initiale 

Vposini1=voisin(center,mat_voi)
for i1 in range(0,6):
    ice_ini[int(Vposini1[i1,0]),int(Vposini1[i1,1]),int(Vposini1[i1,2])]=1
    vap[int(Vposini1[i1,0]),int(Vposini1[i1,1]),int(Vposini1[i1,2])]=2    






#Initalisation de l etat de la frontiere 
fronpos_ini=frontiere(ice_ini,mat_voi,dim_box)

fron_state=np.zeros([np.shape(fronpos_ini)[0],1])

fron_state=np.concatenate((fronpos_ini,fron_state),axis=1).astype(float)
    
# for pos in fronpos_ini:
#     fron_state.append((np.concatenate((pos,np.array([0])),axis=1)).astype(float))





# SIMULATION
###############################################################################
plt.rcParams["figure.figsize"] = (10,10)
if True:
    framev=[]
    frameice=[]
    framet=[]
    #framef=[]
    vap_ini=vap
    frameice.append(ice_ini.copy())
    framev.append(vap.copy())  #juste pour pas avoir tout commenter à chaque fois que je change de quoi
    continuer1=True
    i3=0    
    while continuer1:
        print(i3)
        print(np.sum(ice_ini))
        print()
        ice_c=ice_ini
        
        continuer2=True
        i4=0
        t1=time.time()
        while  continuer2:
            vap_c=relax1(vap_ini,ice_c,fron_state,dim_box,mat_voi,sigma_limit,Pos)
            i4=i4+1
            #relax1(vap,ice,fron_state,dim_box,mat_voi,Pos):
            if ((abs((vap_c-vap_ini)/vap_ini))<0.01).all() : #and i4>=10:
                framev.append(vap_c.copy())
                vap_ini=vap_c
                continuer2=False
                print('nbr iteration pour convergence:')
                print(i4)
                print()
                t2=time.time()
                print('temp convergence:')
                print(t2-t1)
                print()
            else:
                vap_ini=vap_c

        t_min=tmin(fron_state,ice_c,vap_c,nu_kin,mat_voi,dim_box,Lcell)  
       # print(t_min)tmin(fron_state,ice,vap,nu_kin,mat_voi,dim_box,Lcell):
        fron_state_c=croissance(ice_c,fron_state,vap_c,delta_dim,mat_voi,dim_box,nu_kin,t_min)
        #(ice,fron_state,vap,delta_dim,mat_voi,dim_box,nu_kin,t_min)
        update=update_fron(fron_state_c,ice_c,mat_voi,dim_box,Lcell)
        # update_fron(fron_state,ice,mat_voi,dim_box,Lcell)
        
        fron_state=update[0]
        ice_ini=update[1]
        
        #crl !!!!!!!!!!!!!!!!!!!!
        #framef.append(fron_state.copy())
        frameice.append(ice_ini.copy())
        framet.append(t_min)
        i3=i3+1
        if (fron_state==np.array([dim_box[0]-1,dim_box[1]-1,dim_box[2],np.nan])).any() or (fron_state==np.array([1,1,0,np.nan])).any():
            print('fin')
            continuer1=False
        








