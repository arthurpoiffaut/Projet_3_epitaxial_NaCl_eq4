function [in] = in_box( i,j,k,dimx,dimy,dimz)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    if i > 0 && i <= dimx & j > 0 && j <= dimy && k>0 && k<=dimz;
        in = 1 ;   
    else
        in = 0  ;
        
    end
end

