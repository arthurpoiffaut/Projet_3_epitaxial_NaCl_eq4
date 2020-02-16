import numpy as np
import matplotlib.pyplot as plt
import time as time

TempsI=time.clock()#Temps initial
tabM=[2500,5000,10000,20000]#Tableau pour différents nombre de marcheur
plt.figure(4,figsize=(10,5))


for i in tabM:#Boucle pour tracer figure avec différent nombres de marcheurs

    plt.figure(1,figsize=(10,5))
    
    N=256 # taille du réseau
    M=i # nombre de marcheurs
    nIterMax=100000 # nombre maximal d'itérations temporelles
    
    x=np.random.random_integers(1,N,M) # positions initiales des marcheurs
    y=np.random.random_integers(1,N,M) # positions initiales des marcheurs
    
    statusMobile=np.ones(M,dtype='bool') # True pour marcheur mobile
    grilleFixe=np.zeros([N+2,N+2],dtype='bool') # True pour marcheur fixe
    grilleFixe[:,0]=True #sites collants a y=0, au bas du reseau
    
    TabFixe=[]#Liste pour ajouter nombre marcheur fixe
    nFixe=0 # nombre de marcheurs fixes
    iTer=0 # nombre d'itération
    while (nFixe<M) and (iTer<nIterMax): # boucle temporelle, i.e. sur les pas
        m,=statusMobile.nonzero() # indices des marcheurs mobiles
        
        condition=np.random.random(m.size).round()#Condition pour avoir 0 ou 1
        dirX,=np.where(condition==0)
        dirY,=np.where(condition==1)# choix d'un pas vers un des 4 sites voisins
        pas=np.random.random(m.size).round()*2-1
        x[m[dirX]]=np.clip(x[m[dirX]]+pas[dirX],1,N)#On ajoute le pas avec clip
        y[m[dirY]]=np.clip(y[m[dirY]]+pas[dirY],1,N)#pour éviter sortir de grille

 #pour les marcheurs venant de se déplacer, on vérifie si un voisin est collant
 #(arithm. booleenne: True+True=True, True+False=True, et False+False=False)
        voisinFixe=grilleFixe[x[m]-1,y[m]-1]+\
        grilleFixe[x[m] ,y[m]-1]+\
        grilleFixe[x[m]+1,y[m]-1]+\
        grilleFixe[x[m]+1,y[m] ]+\
        grilleFixe[x[m]+1,y[m]+1]+\
        grilleFixe[x[m] ,y[m]+1]+\
        grilleFixe[x[m]-1,y[m]+1]+\
        grilleFixe[x[m]-1,y[m] ]
 #à ce stade voisinFixe est un tableau 1D contenant autant d'éléments qu'il y a de
 #marcheurs mobiles, chaque élément prenant la valeur True ou False selon que le
 #marcheur ait ou non un voisin fixe et collant
        
        k=m[voisinFixe.nonzero()[0]]#Indice des marcheur qui doivent être fixés
        nFixe+=k.size#On ajuste le nombre de marcheur fixes
        if k.size>0: #Boucle pour changer les marcheurs mobiles à fixe
                grilleFixe[x[k],y[k]]=1    
                statusMobile[k]=0

        TabFixe.append(nFixe)
        iTer+=1
        

    plt.subplot(1,2,1)#On trace les fractales
    plt.imshow(grilleFixe.T,origin='lower')
    plt.title('Fractale avec marcheurs')
    
    plt.subplot(1,2,2)#On trace les courbes de croissance
    plt.plot(range(iTer),TabFixe)
    plt.title('Courbe de croissance pour marcheurs')
    plt.xlabel("Nombre d'itérations effectuées")
    plt.ylabel('Nombre de marcheurs fixe')
    plt.show()
    
TF=time.clock()#Temps final
DelT=TF-TempsI#Temps écoulé
print('Temps écoulé: ',DelT)