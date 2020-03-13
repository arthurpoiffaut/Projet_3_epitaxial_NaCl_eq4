function [Tmax,Err,ci,T] = findTmax(q)
L=0.3; %m ; Épaisseur du mur
k=0.85;%W/(m*K); La conductivité thermique de la brique
si=5.67e-8; %constante de Stefan-Boltzmann
h=20;%

%TL=841.90; NONONONONON
To=293; %oK
% Condition radiative à x=0 (face externe du mur): -k*dT/dx=-si*(T^4-To^4)
% !!! Condition radiative est implementée SEULEMENT sur la face externe du mur !!!
c1=-k;
c2=0;
c3=-si*To^4;
% Condition de Neumann à x=L (face interne du mur): dT/dx=0
d1=-k;
d2=0;
d3=-si*To^4;

N=100;
dx=L/N; % Pas de discrétisation
x=(0:dx:L)';
T=To*ones(size(x)); % Approximation initiale

% Sourse volumique de chaleur q[W/m^3] d'épaisseur dL
% La source est intégrée dans la partie intérieure du mur
dL=0.05;
S=q*exp(-((x-L)/dL).^2);

ci=0;
Err=[];
tol=1e-12;
flag=1;

S=q*exp(-((x-L)/dL).^2);
while (flag==1)
    ci=ci+1;
    M=diag(-2*ones(1,N+1),0)+diag(ones(1,N),-1)+diag(ones(1,N),1);
    
    % condition x=0
    M(1,1)=-3*c1;
    M(1,2)=4*c1;
    M(1,3)=-c1;
    
    % condition x=L
    M(N+1,N+1)=-3*d1;
    M(N+1,N)=4*d1;
    M(N+1,N-1)=-d1;
    
    %Condition radiative est implementée SEULEMENT sur la face externe du mur !!!
    b=(S/k)*dx^2;
    %condition x=0
    b(1)=-2*h*dx*(To-T(1))+2*dx*(T(1)^4*si+c3);%%-2*dx^2*(To-T(1))+2*dx*(T(1)^4*si+c3);
    %conditon x=L
    b(N+1)=-2*h*dx*(To-T(N+1))+2*dx*(T(N+1)^4*si+d3);%-2*dx^2*(To-T(N+1))+2*dx*(T(N+1)^4*si+d3);%=2*dx*(T(N+1)^4*si+d3);%-2*dx^2*(To-T(N+1))+%2*dx*(T(N+1)^4*si+d3);%-2*dx*(T(N+1)^4*si+d3);%2*dx^2*(To-T(N+1))-2*dx*(T(N+1)^4*si+d3);%2*d3*dx;
    
    F=M*T+b;
    Err=[Err sum(abs(F))/(N+1)];
    %display(['Étape=', num2str(ci), '   ;   Err=' num2str(Err(end))])
    if (Err(end)<tol)
        flag=0 ; % metre flag ==0 ? jai la fleme
    end
    
    %Condition radiative est implementée SEULEMENT sur la face externe du mur !!!
    J=M;
    J(1,1)=M(1,1)+T(1)^3*8*dx*si+2*h*dx; %probablement la mon erreure mais je compren pas
    J(N+1,N+1)=M(N+1,N+1)+T(N+1)^3*8*dx*si+2*h*dx;
    
    dT=-J\F;
    T=T+dT;
end

Tmax=max(T);
% fprintf('%d \n',Tmax)

