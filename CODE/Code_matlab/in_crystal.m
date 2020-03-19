function [in_ice] = in_crystal(i,j,k,ice,dimx,dimy,dimz)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if in_box(i,j,k,dimx,dimy,dimz)==1;
    if ice(i,j,k)==1;
        in_ice=1;
    else;
        in_ice=0;
    end
else
    in_ice=0;
end


end

