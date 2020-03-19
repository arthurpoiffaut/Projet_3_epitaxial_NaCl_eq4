function  [Vpos] = voisin( i,j,k,dx,dy,dz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%init Vpos
Vpos=zeros(8,3);
for i1= 1:1:length(dx);
    Vpos(i1,1)=i+dx(i1);
    Vpos(i1,2)=j+dy(i1);
    Vpos(i1,3)=k+dz(i1);
end



end

