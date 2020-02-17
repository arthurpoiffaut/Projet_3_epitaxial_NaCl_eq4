function U = getUp1(Up,gp1,gp,A,M,ksi,dx,dt)

U = ( (M - ksi*dt/(alpha*dx^2))\ ((M+ ((1-ksi)*dt/(alpha*dx^2))*A)*Up - (dt/alpha)*(ksi*gp1 + (1-ksi)*gp) ));