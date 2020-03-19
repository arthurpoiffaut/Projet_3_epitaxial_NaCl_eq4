function [ vapout ] = diffusion_vap( vap,ice,delta_tau,dx,dy,dz,dimx,dimy,dimz)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
vapout=vap;

for i1=1:1:dimx;
    for i2=1:1:dimy;
        for i3=1:1:dimz;
            sum1=0;
            sum2=0;
            if in_crystal(i1,i2,i3,ice,dimx,dimy,dimz)==0;
                V=voisin(i1,i2,i3,dx,dy,dz);
                for i4=1:1:6;
                    if in_crystal(V(i4,1),V(i4,2),V(i4,3),ice,dimx,dimy,dimz) == 1;
                        sum1=sum1+vap(i1,i2,i3);
                    elseif in_box(V(i4,1),V(i4,2),V(i4,3),dimx,dimy,dimz)==0;
                        sum1=sum1+vap(i1,i2,i3);
                    else
                        sum1=sum1+vap(V(i4,1),V(i4,2),V(i4,3));
                    end
                end
                for i5=7:1:8;
                    if in_crystal(V(i5,1),V(i5,2),V(i5,3),ice,dimx,dimy,dimz) == 1;
                        sum2=sum2+vap(i1,i2,i3);
                    elseif in_box(V(i5,1),V(i5,2),V(i5,3),dimx,dimy,dimz)==0;
                        sum2=sum2+vap(i1,i2,i3);
                    else
                        sum2=sum2+vap(V(i5,1),V(i5,2),V(i5,3));
                    end
                end
%                 if sum1~=0;
%                     sum1
%                     sum2
%                end
%             else;
%                 disp(i1)
%                 disp(i2)
%                 disp(i3)
            vapout(i1,i2,i3)=(2/3)*(delta_tau)*sum1+(delta_tau)*sum2+(1-6*delta_tau)*vap(i1,i2,i3);
            end
               
                
            end
            
        end
        
    end
    
end



%
% def diffusion(vap,ice,delta_tau,dx,dy,dz,dimx,dimy,dimz):
%      vapout=vap
%      #i=1
%      for a in range(0,dimx):
%          for b in range(0,dimy):
%              for c in range(0,dimz):
%                  sum1=0
%                  sum2=0
%                  #print(sum1)
%                  #print(sum2)
%                  if not in_crystal(a,b,c,ice,dimx,dimy,dimz):
%                      for d in range(0,6):
%                          #print(vap[[a]+dx[d],b+dy[d],c+dz[d]])
%                          #print(in_crystal(a+dx[d],b+dy[d],c+dz[d],ice,dimx,dimy,dimz))
%                          if in_crystal(a+dx[d],b+dy[d],c+dz[d],ice,dimx,dimy,dimz) == False:
%                              sum1=(sum1+vap[a,b,c])
%                          elif  (inbox(a+dx[d],b+dy[d],c+dz[d],dimx,dimy,dimz))==False:
%                               sum1=(sum1+vap[a,b,c])
%                          else:
%                              sum1=(sum1+vap[a+dx[d],b+dy[d],c+dz[d]])
%                              print(sum1)
%                      for e in range(6,8):
%                         #print(c+dz[e])#(in_crystal(a+dx[e],b+dy[e],c+dz[e],ice,dimx,dimy,dimz))
%                         if  in_crystal(a+dx[e],b+dy[e],c+dz[e],ice,dimx,dimy,dimz):
%                             #print(vap[a+dx[e],b+dy[e],c+dz[e]])
%                             #pas bon doit etre la condition haaaaaaaaaaaa
%                             sum2=sum2+vap[a,b,c]   
%                         elif  (inbox(a+dx[e],b+dy[e],c+dz[e],dimx,dimy,dimz))==False:
%                             sum2=sum2+vap[a,b,c]
%                             #i=i+1
%                             #print(i)
%                         else: 
%                             sum2=sum2+vap[a+dx[e],b+dy[e],c+dz[e]]
%                             #i=i+1
%                             #print(i)
%                      print(sum1)
%                      print(sum2)
%                      vapout[a,b,c]=(2/3)*(delta_tau)*sum1+(delta_tau)*sum2+(1-6*delta_tau)*vap[a,b,c]
%     
%      return vapout