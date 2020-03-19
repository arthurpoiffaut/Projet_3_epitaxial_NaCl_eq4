function [ frontpos ] = frontiere( ice,dx,dy,dz,dimx,dimy,dimz )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


end

% def frontiere(ice,dx,dy,dz,dimx,dimy,dimz):
%     icepos=np.argwhere(ice==1)
%     fronpos=[]
%     for a in range(0,np.shape(icepos)[0]):
%         for b in range(0,8):
%             #print(icepos[a])
%             #print(dx[b])
%             #print(dy[b])
%             #print(dz[b])
%             pos=[icepos[a]+[dx[b],dy[b],dz[b]]] 
%             #print(pos)
%             #print(pos[0][0])
%             #print(pos[0][1])
%             #print(pos[0][2])
%             #print(inbox(pos[0][0],pos[0][1],pos[0][2],dimx,dimy,dimz))
%             if  ((inbox(pos[0][0],pos[0][1],pos[0][2],dimx,dimy,dimz)) == False):
%                 pass
%             elif in_crystal(pos[0][0],pos[0][1],pos[0][2],ice,dimx,dimy,dimz) :
%                 pass
%             else:
%                 #p=((inbox(pos[0][0],pos[0][1],pos[0][2],dimx,dimy,dimz)) == False)
%                 #print(p)
%                 fronpos.append(pos[0])
%     #peut etre retourner une matrice a la place ?
%     return fronpos