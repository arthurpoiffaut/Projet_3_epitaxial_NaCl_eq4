function [ y ] = calculetY( dt,yi_1,yi_2,k,m,alpha,g)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

y=(-(k-2*m/dt^2)*yi_1-(m/dt^2 - alpha/(2*dt))*yi_2 - m*g)/(m/dt^2 + alpha/(2*dt));


end

