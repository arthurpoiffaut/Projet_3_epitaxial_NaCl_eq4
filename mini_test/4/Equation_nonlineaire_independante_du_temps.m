clear
% �quation diff�rentielle: d^2 u/dx^2=g(x) sur x=(a,b)
% Conditions aux limites g�n�rales:
% x=a: c1*du/dx+c2*u+c3=-si*T^4
% x=b: d1*du/dx+d2*u+d3=0

% �quation de transfert de chaleur d^2 T/dx^2=-S(x)/k sur x=(0,L)
% dans un mur d'isolation thermique
L=0.3; %m ; �paisseur du mur
k=0.85;%W/(m*K); La conductivit� thermique de la brique
si=5.67e-8; %constante de Stefan-Boltzmann

To=293; %oK
% Condition radiative � x=0 (face externe du mur): -k*dT/dx=-si*(T^4-To^4)
% !!! Condition radiative est implement�e SEULEMENT sur la face externe du mur !!!
c1=-k; c2=0; c3=-si*To^4;
% Condition de Neumann � x=L (face interne du mur): dT/dx=0
d1=1; d2=0; d3=0;

%(N+1) n�uds dans la maille
% Nmax=10000 pour 1G de m�moire

% Nar1=[100]; %dx=3mm
% %Nar1=[2:10 20:10:100 200:100:2000];
%
% Nar(1:2:2*length(Nar1)-1)=Nar1;
% Nar(2:2:2*length(Nar1))=2*Nar1;

N=100;
dx=L/N; % Pas de discr�tisation
x=(0:dx:L)';
T=To*ones(size(x)); % Approximation initiale

% Sourse volumique de chaleur q[W/m^3] d'�paisseur dL
% La source est int�gr�e dans la partie int�rieure du mur
dL=0.1;
q=2000;% W/m^3;
S=q*exp(-((x-L)/dL).^2);

figure(1);plot(x,T,'r');
xlabel('x [m]');ylabel('T [K]');hold
% boucle de convergence
ci=0;Err=[];tol=1e-12;flag=1;
while (flag==1)
    ci=ci+1;
    
    M=diag(-2*ones(1,N+1),0)+diag(ones(1,N),-1)+diag(ones(1,N),1);
    M(1,1)=2*c2*dx-3*c1;M(1,2)=4*c1;M(1,3)=-c1;
    M(N+1,N+1)=3*d1+2*d2*dx;M(N+1,N)=-4*d1;M(N+1,N-1)=d1;
    
    % !!! Condition radiative est implement�e SEULEMENT sur la face externe du mur !!!
    b=S/k*dx^2; b(1)=T(1)^4*2*dx*si+2*c3*dx; b(N+1)=2*d3*dx;
    
    F=M*T+b;
    Err=[Err sum(abs(F))/(N+1)];
    display(['�tape=', num2str(ci), '   ;   Err=' num2str(Err(end))])
    if (Err(end)<tol) break; end
    pause
    
    % !!! Condition radiative est implement�e SEULEMENT sur la face externe du mur !!!
    J=M;J(1,1)=M(1,1)+T(1)^3*8*dx*si;
    
    dT=-J\F;
    T=T+dT;
    
    plot(x,T,'b')
end
hold

figure(2);semilogy((1:ci),Err,':o');
xlabel('�tape');ylabel('Err');