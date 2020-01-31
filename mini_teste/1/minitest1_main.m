%main

%initialisation
k = 1
g = 1
m = 1
dt= 0.01
T = 100

y0 =1
y1 = (1-(k*dt^2)/2m)*y0 + dt*(1-(alpha*dt/2*m))*v0 - (g*dt^2)/2 ;


y = ones(1,T)
y = calculetY( dt,yi,yi2,k,m,alpha,g)