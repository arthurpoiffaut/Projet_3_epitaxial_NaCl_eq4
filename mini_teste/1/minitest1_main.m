%main

%initialisation
alpha =1;
k = 1;
g = 1;
m = 1;
dt= 0.01;
numberOfSteps = 10;
y0 = 1;
v0 = 1;




y = ones(1,numberOfSteps)
y(1) =y0;
y(2) = (1-(k*dt^2)/2*m)*y0 + dt*(1-(alpha*dt/2*m))*v0 - (g*dt^2)/2 ;

% for i = (1,numberOfSteps)
%     y[i] = calculetY( dt,yi,yi2,k,m,alpha,g)