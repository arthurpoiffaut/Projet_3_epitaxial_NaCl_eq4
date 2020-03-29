# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:32:58 2020

@author: plume_2.0
"""
#Import
import numpy as np
import matplotlib.pyplot as plt 
#import scipy as sp
import scipy.constants as spc


#FONCTIONS
###############################################################################


#Fonctions simples
###############################################################################
#Fonction qui valide si le point est dans la boîte de croissance 

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

#Fonction voisin qui renvoit les voisins d'un point
def voisin( i,j,k,dx,dy,dz):
    Vpos=np.zeros([8,3]).astype(int)
    for a in range(0,8):
        Vpos[a,0]=i+dx[a]
        Vpos[a,1]=j+dy[a]
        Vpos[a,2]=k+dz[a]
    return Vpos
#pourrait simplifier la lecture du code je vais pas reecrire toutes les fonctions avec mais 
# je le ferai eventuellement

#Fonction valide si le point donné est dans le cristal     
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

#Fonction qui retourne la frontiere en transition 
#def frontiere(ice,dx,dy,dz,dimx,dimy,dimz):
#    icepos=np.argwhere(ice==1)
#    fronpos=[]
#    for a in range(0,np.shape(icepos)[0]):
#        Vpos=voisin(icepos[a][0],icepos[a][1],icepos[a][2],dx,dy,dz)
#        for b in range(0,8):
#            #print(infron(Vpos[b][0],Vpos[b][1],Vpos[b][2],fronpos))
#            if  ((inbox(Vpos[b][0],Vpos[b][1],Vpos[b][2],dimx,dimy,dimz)) == False):
#                pass
#            elif a==0 and b==0:
#                fronpos.append(Vpos[b])
##            elif infron(Vpos[b][0],Vpos[b][1],Vpos[b][2],fronpos):
##                fronpos.pop()
##                print("pop")
#            elif  ((np.sum((fronpos)==Vpos[b],axis=1)==3).any()==True):
#                pass
#                #print("hello")
#            elif in_crystal(Vpos[b][0],Vpos[0][1],Vpos[b][2],ice,dimx,dimy,dimz) :
#                pass
#          #plus pop^^^^
#
#            else:
#                #print(Vpos[b])
#                fronpos.append(Vpos[b])
#                #print(np.sum((np.array(fronpos)==pos).all(axis=1))>1)
#                #print(a)
#            
#           
#            
#            #if  ((np.sum((fronpos)==Vpos[b],axis=1)==3).any()==True) :
#            #    print(fronpos.pop())
#                #fronpos.pop()
#    
#    
#    return fronpos


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

#Fonction in frontiere
def infron(i,j,k,fronpos):
    infr=[]
    pos=[]
    for a in range(0,np.shape(fronpos)[0]):
        if fronpos[a][0]==i and  fronpos[a][1]==j and  fronpos[a][2]==k:
            infr.append(True)
            pos=[fronpos[a][0],fronpos[a][1],fronpos[a][2],a]
        else:
            infr.append(False) 
    if sum(infr)>=1:
        return [pos[0],pos[1],pos[2],pos[3]]
    else:
        return 0



#Fonction qui calcule la constante de drain
def K(alpha,D,delta_dim,T,m):
    K=alpha*((2*delta_dim)/(np.sqrt(3)*D))*np.sqrt((spc.k*T)/(2*spc.pi*m))
    return K

#Fonction qui calcule nu_kin
def nukin(ratio_c,T,m):
    nu_kin=ratio_c*np.sqrt((spc.k*T)/(2*spc.pi*m))
    return nu_kin
    
    
#fonction qui donne le bon coefficent de alpha en fonction des voisins
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
        alpha=0.1 # valeur pour vertical
    elif V<=1 and H==1:
        alpha=0.5 # valeur pour horizontal pour les autre  (pas sur revoir article)
    elif H==2 and V==0:
        alpha=0.2
    elif H==0 and V==0:
        alpha=0
    else:
        alpha=1
    #print("H")
    #print(H)
    #print("V")
    #print(V)
    return alpha 
#Note: peut probablement etre amélioré genre juste compter des vrais et des faux ?



#Fonction qui calcule deltaV
def deltaV(alpha,Cvap,nu_kin,delta_dim,delta_t):
    acc=1
    dV=acc*alpha*nu_kin*Cvap*(delta_dim**2)*delta_t
    return dV


#Fonction de la vitesse de surface normale

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
     #pas sur i des problemes avec Cvap des fois je mets juste pour voir ce que ca 
    #donne passen le probleme des vapeurs nulles pour plus tard( pou nu_n)
    dL=nu_n*temp
    if nu_n==0: #pas sur brian mais ca devrait marcher
        #print("temp")
        #print(temp)
        #print("v")
        #print(nu_n)
        dL=1
    return dL

#Fonction qui calule le temps de croissance pour qu'une cellule se remplisse
def tcroi(longeur,nu_n): #longueur restant à croître
    #pas sur i des problemes avec Cvap des fois je mets juste pour voir ce quecsa 
    #donne passen le probleme des vapeurs nulles pour plus tard( pou nu_n)
    if nu_n==0:
        print("ha")
        tc=0
    else:
        tc=(longeur)/nu_n
    return tc
#nan


#Fonction pour croissance
###############################################################################


         

#Fonction de relaxation 1, alpha independant de sigma 
#condition aux frotières de la boîte, sont une source constante
def relax1(vap,ice,delta_tau,delta_dim,D,m,sigma_limit,dx,dy,dz,dimx,dimy,dimz):
    vap_out=np.zeros(np.shape(vap))
    fronpos=frontiere(ice,dx,dy,dz,dimx,dimy,dimz)
    for a in range(0,dimx):
        for b in range(0,dimy):
            for c in range(0,dimz):
                sum1=0
                sum2=0
                #if in_crystal(a,b,c,ice,dimx,dimy,dimz)== True:
                 #   print(a)
                  #  print(b)
                   # print(c)
                if in_crystal(a,b,c,ice,dimx,dimy,dimz)== False:
                    if infron(a,b,c,fronpos)!=0:
                       # print("True")
                        Vpos=voisin(a,b,c,dx,dy,dz)
                        alpha=valapha(a,b,c,dx,dy,dz,dimx,dimy,dimz,ice)
                        #print(alpha)
                        Kval=K(alpha,D,delta_dim,T,m)
                        for d in range(0,6):
                            if in_crystal(Vpos[d][0],Vpos[d][1],Vpos[d][2],ice,dimx,dimy,dimz) == True:
                                sum1=sum1+ +vap[a,b,c]
                            elif  inbox(Vpos[d][0],Vpos[d][1],Vpos[d][2],dimx,dimy,dimz) == False:
                                sum1=sum1+vap[a,b,c]
                            else:
                                sum1=sum1+vap[int(Vpos[d][0]),int(Vpos[d][1]),int(Vpos[d][2])]
                        for e in range(6,8):
                            if in_crystal(Vpos[e][0],Vpos[e][1],Vpos[e][2],ice,dimx,dimy,dimz) == True:
                                sum2=sum2+ +vap[a,b,c]
                            elif inbox(Vpos[e][0],Vpos[e][1],Vpos[e][2],dimx,dimy,dimz) == False:
                                sum2=sum2+vap[a,b,c]
                            else:
                                sum2=sum2+vap[int(Vpos[e][0]),int(Vpos[e][1]),int(Vpos[e][2])]
                        vap_out[a,b,c]=((2/3)*sum1+sum2)/(Kval+6)
                            
                    elif frontbox(a,b,c,dimx,dimy,dimz) == True:
                        vap_out[a,b,c]=sigma_limit
                    else:
                        #copy()
                        #print("oui")
                        Vpos=voisin(a,b,c,dx,dy,dz)
                        for f in range(0,6):
                            if in_crystal(Vpos[f][0],Vpos[f][1],Vpos[f][2],ice,dimx,dimy,dimz) == True:
                                sum1=sum1+ +vap[a,b,c]
                            elif  inbox(Vpos[f][0],Vpos[f][1],Vpos[f][2],dimx,dimy,dimz)==False:
                                sum1=sum1+vap[a,b,c]
                            else:
                                #print("oui")
                                sum1=sum1+vap[int(Vpos[f][0]),int(Vpos[f][1]),int(Vpos[f][2])] 
                            
                        for g in range(6,8):
                             if in_crystal(Vpos[g][0],Vpos[g][1],Vpos[g][2],ice,dimx,dimy,dimz) == True:
                                sum2=sum2+ +vap[a,b,c]
                             elif inbox(int(Vpos[g][0]),int(Vpos[g][1]),int(Vpos[g][2]),dimx,dimy,dimz) == False:
                                sum2=sum2+vap[a,b,c]
                             else:
                                sum2=sum2+vap[int(Vpos[g][0]),int(Vpos[g][1]),int(Vpos[g][2])]
                        #vap_out[a,b,c]=(2/3)*(delta_tau)*sum1+(delta_tau)*sum2+(1-6*delta_tau)*vap[a,b,c]
                        vap_out[a,b,c]=((2/3)*sum1+sum2)/6
                else:
                    vap_out[a,b,c]=np.nan

    return vap_out
                    
#Fonction de relaxation 2, alpha dependant de sigma 
#Condition aux frontières de la boîte sont une source constante??




#Fonction qui trouve le temps minimal pour que l'une des cellules se remplisse
def tmin(fron_state,delta_dim,ratio_c,T,m,dx,dy,dz,dimx,dimy,dimz,ice,Lcell):
    list_tc=[]
    for a in range(0,np.shape(fron_state)[0]):
        longeur0=fron_state[a][3]
        longeur=Lcell-longeur0
        Cvap=vap[int(fron_state[a][0]),int(fron_state[a][1]),int(fron_state[a][2])]
        alpha=valapha(int(fron_state[a][0]),int(fron_state[a][1]),int(fron_state[a][2]),dx,dy,dz,dimx,dimy,dimz,ice)
        nu_kin=nukin(ratio_c,T,m)
        nu_n=nun(alpha,Cvap,nu_kin)
        tc=tcroi(longeur,nu_n)
        #print(tc)
        list_tc.append(tc)
    #print(list_tc)
    t_min=min(list_tc)
    return t_min
#dLcell


def croissance(dx,dy,dz,dimx,dimy,dimz,ice,fron_state,vap,ratio_c,T,m,delta_dim,t_min):
    fron_state_c=fron_state.copy()
    fron_state_ini=fron_state.copy()
    for a in range(0,np.shape(fron_state)[0]):
        alpha=valapha(int(fron_state[a][0]),int(fron_state[a][1]),int(fron_state[a][2]),dx,dy,dz,dimx,dimy,dimz,ice)
        #print("postion:")
        #print(fron_state[a])
        #print("alpha")
        #print(alpha)
        Cvap=vap[int(fron_state[a][0]),int(fron_state[a][1]),int(fron_state[a][2])]
        #print("Cvap")
        #print(Cvap)
        nu_kin=nukin(ratio_c,T,m)
        nu_n=nun(alpha,Cvap,nu_kin)
        #t_min=tmin(fron_state,delta_dim,ratio_c,T,m,dx,dy,dz,dimx,dimy,dimz,ice)
        dL=delta_L(nu_n,t_min)
        #print(dL)
        #print(dL)
        #print(fron_state[a][3]+a*1000000*dV)
        #new_value=fron_state[a][3]+dV
        #print(dV)
        fron_state_c[a][3]=(fron_state_ini[a][3]+dL).copy()
        #print(fron_state[a][3]+dL)
    #fron_state_c=fron_state
    return fron_state_c
            
#fron_stat plus sur si les positions de dans c'est un peu con je suis pas sur
    



#Fonction qui update la matrice ice et la liste frontiere
def update_fron(fron_state,Lcell,ice,dx,dy,dz,dimx,dimy,dimz):
    fron_state_new=[]
    ice_new=ice
    #p=[]
    for a in range(0,np.shape(fron_state)[0]):
        #print(a)
        #print(np.shape(fron_state)[0])
        #print(fron_state[a][3])
        #if fron_state[a][2]==4 or fron_state[a][2]==6:
                #print(fron_state[a])
        if np.round(fron_state[a][3],8)>=Lcell:
            #if fron_state[a][2]==4 or fron_state[a][2]==6:
                #print(fron_state[a])
                
            #print(fron_state[a][3])
            #print("ice")
            #print()
            #print(fron_state[a])
            #p.append(a)
            ice_new[int(fron_state[a][0]),int(fron_state[a][1]),int(fron_state[a][2])]=1
            
#            Vpos=voisin(int(fron_state[a][0]),int(fron_state[a][1]),int(fron_state[a][2]),dx,dy,dz)
#            for b in range(0,6):
#                if (infron(Vpos[b][0],Vpos[b][1],Vpos[b][2],fronpos)== False) or (in_crystal(Vpos[b][0],Vpos[b][1],Vpos[b][2],ice_new,dimx,dimy,dimz)==False):
#                    fron_state_new.append(np.array([Vpos[b][0],Vpos[b][1],Vpos[b][2],0]))
#                    
#            for c in range(6,8):
#                if (infron(Vpos[c][0],Vpos[c][1],Vpos[c][2],fronpos)== False) or (in_crystal(Vpos[c][0],Vpos[c][1],Vpos[c][2],ice_new,dimx,dimy,dimz)==False):
#                    fron_state_new.append(np.array([Vpos[c][0],Vpos[c][1],Vpos[c][2],0]))
#            
    fronpos=frontiere(ice_new,dx,dy,dz,dimx,dimy,dimz)
    #print(np.shape(fronpos)[0])
    #print(np.shape(fron_state)[0])
    for b in range(0,np.shape(fronpos)[0]):
        pos=infron(fronpos[b][0],fronpos[b][1],fronpos[b][2],fron_state)
        if  pos !=0 :
            #print(pos)
            fron_state_new.append([fronpos[b][0],fronpos[b][1],fronpos[b][2],fron_state[pos[3]][3]])
            #print([fronpos[b][0],fronpos[b][1],fronpos[b][2],fron_state[pos[3]][3]])
        else:
            fron_state_new.append([fronpos[b][0],fronpos[b][1],fronpos[b][2],0])
    #fron_state_new=frontiere(ice_new,dx,dy,dz,dimx,dimy,dimz)
    #for d in p:
        #ice[int(fron_state[d][0]),int(fron_state[d][1]),int(fron_state[d][2])]=1
        #print("ice")
        #print(fron_state_new[d])
        #fron_state_new[d]=1
    return [fron_state_new,ice_new]



#Prends la matrice 2D (''écrasée'') contenant l'image du flocon et l'hexagonifie
#Sauvegarde aussi l'image obtennue au même endoit où le code .py est enregistré
def hexplot(matflocon):
    hexagonal = np.empty((2,dimx,dimy)) #array contenant les points hex
    n = np.array([1., -1./np.sqrt(3.)]) #vecteurs de la base hexagonale
    m = np.array([1., 1./np.sqrt(3.)])
    for i in range(dimx):
        for j in range(dimy):
            hexagonal[:,i,j] = int(i)*n + int(j)*m
    index_flocon = np.where(matflocon >= 1)  #0 sion veut tous les points de la matrice.
    flocon = hexagonal[:,index_flocon[0],index_flocon[1]]
    fig = plt.figure()
    plt.scatter(flocon[0,:], flocon[1,:], c =matflocon[index_flocon[0],\
                 index_flocon[1]]-1., s = 500., \
                    cmap = 'cubehelix', marker = 'H') # s est la taille des points scatter
    plt.axis("off")
    fig.savefig('testfloconhex_'+'mx'+str(dimx)+'my'+str(dimy)+'.png')
    plt.show()

# faudraiT voir  ce qu'on peut faire avec cela
#Fonction pour un affichage en hexagone #2 
def hexplot2(frameflocon,dimx,dimy,name):
    #xi=np.linspace(0,dimx-1,dimx)
    #yi=np.linspace(0,dimy-1,dimy)
    hexagonal = np.empty((2,dimx,dimy)) #array contenant les points hex
    n = np.array([1., -1./np.sqrt(3.)]) #vecteurs de la base hexagonale
    m = np.array([1., 1./np.sqrt(3.)])
    for a in range(0,np.shape(frameflocon)[0]):
        x=[]
        y=[]
        z=[]
        f='frame'+str(a)
        #matflocon=frameflocon[a][:,:,7]  # pour la vapeur dernier depen du centre
        matflocon=np.sum(frameflocon[a],axis=2) #pour la glace
        for b in range(0,dimx):
            for c in range(0,dimy):
                hexagonal[:,b,c] = int(b)*n + int(c)*m
                x.append( hexagonal[0,b,c])
                y.append( hexagonal[1,b,c])
                z.append(matflocon[b,c])
        
        index_flocon = np.where(matflocon >= 1) #depende de ce qu on veu prendre
         
        maxx=max(hexagonal[0,index_flocon[0],index_flocon[1]])
        maxy=max(hexagonal[1,index_flocon[0],index_flocon[1]])
        
        minx=min(hexagonal[0,index_flocon[0],index_flocon[1]])
        miny=min(hexagonal[1,index_flocon[0],index_flocon[1]])
        
        fig = plt.figure()
        
        plt.hexbin(x,y,z,gridsize=(dimx-1,dimy-1))
        plt.xlim(minx-5,maxx+5)
        plt.ylim(miny-5,maxy+5)
        fig.savefig('testfloconhex_'+name+'mx'+str(dimx)+'my'+str(dimy)+f+'.png')



###############################################################################
#Initialisation
###############################################################################





#Initiation des vecteurs pour regarder les plus proches voisins
dx = [-1,0,-1,1,0,1,0,0]
dy = [-1,-1,0,0,1,1,0,0]
dz = [0,0,0,0,0,0,1,-1]


# Constantes :

delta_dim=30*10**(-6); # valeur de la dimention (m)

delta_t=10*10**(-9); # valeur de la variation de temps 

T=258 #temperature en kelvin

m=2.988*10**(-26)#masse d'une molecule d'eau

D=2*10**-5; #valeur du coe de diff   (m^2/s)

Vcell=(np.sqrt(3)/2)*(delta_dim**3);

#nu_kin=133; # micrometre/s en fonction aussi

ratio_c=1*10**(-6) #pas sur ici faudrait revoir c_sat/c_solide

#c_solide=1*10**(-6)

sigma_limit=1 #concentration aux extremitées

#longueur d'une cellule 
Lcell=np.round((np.sqrt(3)/2)*(delta_dim),8)

delta_tau=(D*delta_t)/(delta_dim)**2
# att on doit le faire tendre vers 0 cest con 
# ca va etre fucking long a rouler....

#Dimension en x, y et z du cristal
dimx=100
dimy=100
dimz=13
#=======
#dimx=100; #dimention x  ect
#dimy=100;
#dimz=4;




#Definition de matrices 

# mat cristal 0 si pas dans le cristal
ice_ini=np.zeros([dimx,dimy,dimz])
# mat vapeur varie en fonction de la densité d'eau en état vapeur entre 0 et 1 ?
vap=np.zeros([dimx,dimy,dimz])
#mat de la frontiere
#fron=np.zeros([dimx,dimy,dimz])


#État initial
###############################################################################

vap[:,:,:]=1
#Point au centre du cristal (seed)
cx=int((dimx+1)/2)
cy=int((dimy+1)/2)
cz=int((dimz+1)/2)
#un hexagone au centre 

#Essai pour mettre un hexagone au centre. Il s'agirait d'un hexagone sur grille carrée(cubique)
#et donc décalé. Pas un vrai hexagone car la grille n'est pas hexagonifiée à l'état initial
#Forme hexagonale que en x et y. Épaisseur en z constante 
#Voir fichier FonctionHex.py

#ice_ini[cx,cy,cz]=1
#ice_ini[cx-1,cy,cz]=1
#ice_ini[cx-1,cy-1,cz]=1
#ice_ini[cx,cy-1,cz]=1
#ice_ini[cx+1,cy,cz]=1
#ice_ini[cx+1,cy+1,cz]=1
#ice_ini[cx,cy+1,cz]=1


ice_ini[cx,cy,cz]=1
vap[cx,cy,cz]=0
Vposini1=voisin(cx,cy,cz,dx,dy,dz)
# peut les changer pour des listes directement

for i1 in range(0,6):
    ice_ini[int(Vposini1[i1,0]),int(Vposini1[i1,1]),int(Vposini1[i1,2])]=1
    Vposini2=voisin(int(Vposini1[i1,0]),int(Vposini1[i1,1]),int(Vposini1[i1,2]),dx,dy,dz)    
    for i2 in range(0,6):
        ice_ini[int(Vposini2[i2,0]),int(Vposini2[i2,1]),int(Vposini2[i2,2])]=1
        Vposini3=voisin(int(Vposini2[i2,0]),int(Vposini2[i2,1]),int(Vposini2[i2,2]),dx,dy,dz)
        for i3 in range(0,6):
            ice_ini[int(Vposini3[i3,0]),int(Vposini3[i3,1]),int(Vposini3[i3,2])]=1
            Vposini4=voisin(int(Vposini3[i3,0]),int(Vposini3[i3,1]),int(Vposini3[i3,2]),dx,dy,dz)
            for i4 in range(0,6):
                ice_ini[int(Vposini4[i4,0]),int(Vposini4[i4,1]),int(Vposini4[i4,2])]=1
                Vposini5=voisin(int(Vposini4[i4,0]),int(Vposini4[i4,1]),int(Vposini4[i4,2]),dx,dy,dz)
                for i5 in range(0,6):
                    ice_ini[int(Vposini5[i5,0]),int(Vposini5[i5,1]),int(Vposini5[i5,2])]=1
                   
        
        
        
#ice_ini=ice
#pas trop sur pour la condition initiale de la glace mais on verra bien 

#Initalisation de la frontiere 
fronpos_ini=frontiere(ice_ini,dx,dy,dz,dimx,dimy,dimz)
fron_state=[]

for i2 in range(0,np.shape(fronpos_ini)[0]):
    fron_state.append((np.concatenate((fronpos_ini[i2],np.array([0])),axis=0)).astype(float))








#test

#ice[1,1,1]=1
#ice[1,2,1]=1
#print(K(1,D,delta_dim,T,m)*delta_tau)

###############################################################################
#Simulation
###############################################################################

plt.rcParams["figure.figsize"] = (10,10)
if False:
    framev=[]
    frameice=[]
    framet=[]
    framef=[]
    frameice.append(ice_ini.copy())
if False: #juste pour pas avoir tout commenter à chaque fois que je change de quoi
    for i3 in range(0,300):
        print(i3)
        print(np.sum(ice_ini))
        print()
        ice_c=ice_ini
        vap_ini=vap
        continuer=True
        i4=0
        while  continuer:
            vap_c=relax1(vap_ini,ice_c,delta_tau,delta_dim,D,m,sigma_limit,dx,dy,dz,dimx,dimy,dimz)
            
            i4=i4+1
            if (abs(vap_c-vap_ini).all()<0.0005) and i4>=5:
                framev.append(vap_c.copy())
                vap_ini=vap_c
                continuer=False
                print()
                print(i4)
                print()
            else:
                vap_ini=vap_c

        t_min=tmin(fron_state,delta_dim,ratio_c,T,m,dx,dy,dz,dimx,dimy,dimz,ice_c,Lcell)  
        
        fron_state_c=croissance(dx,dy,dz,dimx,dimy,dimz,ice_c,fron_state,vap,ratio_c,T,m,delta_dim,t_min)
        
        update=update_fron(fron_state_c,Lcell,ice_c,dx,dy,dz,dimx,dimy,dimz)
        
        fron_state=update[0]
        ice_ini=update[1]
        
        #crl !!!!!!!!!!!!!!!!!!!!
        framef.append(fron_state.copy())
        frameice.append(ice_ini.copy())
        framet.append(t_min.copy())
        
        if (fron_state==np.array([dimx-2,dimy-2,dimz-2,np.nan])).any() or (fron_state==np.array([1,1,1,np.nan])).any():
            break





#Image
###############################################################################



#for i in range(0,200):
#    vap=relax1(vap,ice,delta_tau,delta_dim,D,m,sigma_limit,dx,dy,dz,dimx,dimy,dimz)
plt.rcParams["figure.figsize"] = (10,10)
#plt.imshow(np.sum(ice,axis=2),interpolation='spline16', cmap='viridis')

#floc=np.sum(ice,axis=2)
#figg=plt.figure()
#plt.imshow(floc,interpolation='spline16', cmap='Blues')



