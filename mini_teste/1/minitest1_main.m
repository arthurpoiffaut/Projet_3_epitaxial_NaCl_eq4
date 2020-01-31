%main

%initialisation
k = 1
g = 1
m = 1
dt= 0.01
numberOfSteps = 100
y0 = 1
v0k = 1




y = ones(1,numberOfSteps)
y[1] =1
y[2] = (1-(k*dt^2)/2m)*y0 + dt*(1-(alpha*dt/2*m))*v0 - (g*dt^2)/2 ;

for i = (1,numberOfSteps)
y[i] = calculetY( dt,yi,yi2,k,m,alpha,g)