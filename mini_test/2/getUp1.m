function U = getUp1(Up,b,A,M,ksi,dx,dt,alpha)

U1 =  (M - ksi.*dt./(alpha*dx^2).*A);
U2 =   (M+ ((1-ksi).*dt./(alpha*dx^2)).*A)*Up ;
U3 = (dt/alpha).*(b) ;

U = U1\(U2-U3);

end