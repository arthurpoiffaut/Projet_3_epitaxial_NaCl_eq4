function [ y ] = calculetY( dt,yi_1,yi_2,k,m,alpha,g)
%This function calculates one itération of the y_i values


y=((2*m/dt^2 - k)*yi_1 + ((alpha/2 - m/dt)/dt)*yi_2 - m*g) *  (dt/( alpha/2 + m/dt));


end

