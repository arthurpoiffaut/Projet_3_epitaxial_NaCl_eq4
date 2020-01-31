function [ y ] = calculetY( dt,yi,yi2,k,m,alpha,g)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

y=-(k-2m/dt)yi2-(m/dt^2-alpha/(2*dt))yi-m*g


end

