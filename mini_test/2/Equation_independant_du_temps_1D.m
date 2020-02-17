clear 
% Équation différentielle: d^2 u/dx^2=g(x) sur x=(a,b)
% Conditions aux limites générales:
% x=a: c1*du/dx+c2*u+c3=0
% x=b: d1*du/dx+d2*u+d3=0

% Équation de transfert de chaleur d^2 T/dx^2=-S(x)/k sur x=(0,L)
% dans un mur d'isolation thermique
L=0.3; %m ; Épaisseur du mur

k=1;h=10; %Mini quiz 4, 2016
% k=1;%W/(m*K); La conductivité thermique de la brique
% h=1; %W/(m^2*K); Coefficient de transfert thermique pour l'interface plane entre l'air et solide.


ksi = 0.5;
Cv = 1000;
rho = 2000;

dL=0.10; 
q=2000;% W/m^3;



% Condition convective (de Robin) à x=0 (face externe du mur): -k*dT/dx=h(Ta-T)
Ta=-10; %oC
c1=-k; c2=h; c3=-h*Ta;
% Condition de Robin
Ti=20;
d1=-k; d2=-h; d3=h*Ti;

%(N+1) nœuds dans la maille
% Nmax=10000 pour 1G de mémoire


Nar1=[100]; %dx=3mm
%Nar1=[2:10 20:10:100];% 200:100:2000 3000:1000:5000]; % Matrice pleine
%Nar1=[2:10 20:10:100 200:100:2000 3000:1000:4000 5000:5000:100000]; %Matrice creuse

Nar(1:2:2*length(Nar1)-1)=Nar1;
Nar(2:2:2*length(Nar1))=2*Nar1;


%%
ci=0;Tmax=[];figure(1);time=[];
for N=Nar
    display(N)
    ci=ci+1;
    dx=L/N; %Pas de discrétisationà
    if N ==100
        dxReal = dx;
    end
    x=(0:dx:L)';
    
    % Sourse volumique de chaleur q[W/m^3] d'épaisseur dL
    % La source est intégrée dans la partie intérieure du mur
    S=q.*heaviside(x-(L-dL));
    
    if (0)
        % matrice pleine
        A=diag(-2*ones(1,N+1),0)+diag(ones(1,N),-1)+diag(ones(1,N),1);
    else
        
        % matrice creuse
        i=[(1:N+1) (1:N) (2:N+1)];
        j=[(1:N+1) (2:N+1) (1:N)];
        s=[-2*ones(1,N+1) ones(1,N) ones(1,N)];
        A=sparse(i,j,s,N+1,N+1);
    end
    A(1,1)=2*c2*dx-3*c1;A(1,2)=4*c1;A(1,3)=-c1;
    A(N+1,N+1)=3*d1+2*d2*dx;A(N+1,N)=-4*d1;A(N+1,N-1)=d1;
    b=-S/k*dx^2; b(1)=-2*c3*dx; b(N+1)=-2*d3*dx;
    
    if N ==100
        AReal = A;
    end
    
    tic
    u=A\b;
    time=[time toc];
    
    plot(x,u); if (ci==1) hold; end
    Tmax=[Tmax max(u)];
end

axis([x(1) x(end) Ta 30])
xlabel('x [m]')
ylabel('T_{eq}(x) [^oC]')
hold

Err=abs(Tmax(1:2:end)-Tmax(2:2:end));

figure(2)
loglog(L./Nar1,Err,'o')
xlabel('dx')
ylabel('Err(dx)=|Tmax(dx)-Tmax(dx/2)|')

figure(3)
loglog(Nar,time,':o')
xlabel('N')
ylabel('temps [s]')

Tmax_eq=Tmax(end-1);

%%
M = diag(ones(1,N+1));
M(1,1) = 0;
M(end,end) = 0;
M = sparse(M);

A = AReal
dx = dxReal

alpha = Cv*rho/k
dt = 1* (alpha*dx^2)




