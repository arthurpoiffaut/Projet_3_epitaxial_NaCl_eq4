function [y,time] = calculetAllY( y,time,y0,v0,dt,T,k,m,alpha,g,timeMultiplier)
% This function thought all dt's and defines a y and tie vector for all of
% them.


for i1 = 1:length(dt)
        t=(0:dt(i1):timeMultiplier*T);
        time{i1}=t;
        y{i1}=zeros(length(t),1);
        y{i1}(1)=y0;
        y{i1}(2)=(1-(k*dt(i1)^2)/(2*m))*y0 + dt(i1)*(1-(alpha*dt(i1)/(2*m)))*v0 - (g*dt(i1)^2)/2 ;
        for i2= 3:length(t)
            y{i1}(i2) = calculetY( dt(i1),y{i1}(i2-1),y{i1}(i2-2),k,m,alpha,g);
        end
    end