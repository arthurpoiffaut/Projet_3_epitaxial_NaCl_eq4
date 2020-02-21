clear 
% ï¿½quation diffï¿½rentielle: d^2 u/dx^2=g(x) sur x=(a,b)
% Conditions aux limites gï¿½nï¿½rales:
% x=a: c1*du/dx+c2*u+c3=0
% x=b: d1*du/dx+d2*u+d3=0

% ï¿½quation de transfert de chaleur d^2 T/dx^2=-S(x)/k sur x=(0,L)
% dans un mur d'isolation thermique
L=0.3; %m ; ï¿½paisseur du mur

<<<<<<< HEAD
k=1;h=10; %Mini quiz 4, 2016
% k=1;%W/(m*K); La conductivitï¿½ thermique de la brique
=======
k=1;h=30; %Mini quiz 4, 2016
% k=1;%W/(m*K); La conductivité thermique de la brique
>>>>>>> master
% h=1; %W/(m^2*K); Coefficient de transfert thermique pour l'interface plane entre l'air et solide.
ksi = 0.65;
rho = 3000;
dL=0.20;

<<<<<<< HEAD
<<<<<<< HEAD
% Condition convective (de Robin) ï¿½ x=0 (face externe du mur): -k*dT/dx=h(Ta-T)
=======

ksi = 0.5;
Cv = 1000;
rho = 2000;

dL=0.10; 
=======
 
>>>>>>> master
q=2000;% W/m^3;
Cv = 1000;


% Condition convective (de Robin) à x=0 (face externe du mur): -k*dT/dx=h(Ta-T)
>>>>>>> master
Ta=-10; %oC
c1=-k; c2=h; c3=-h*Ta;
% Condition de Robin
Ti=20;
d1=-k; d2=-h; d3=h*Ti;

%(N+1) nï¿½uds dans la maille
% Nmax=10000 pour 1G de mï¿½moire


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
<<<<<<< HEAD
    dx=L/N; %Pas de discrï¿½tisation
=======
    dx=L/N; %Pas de discrétisationà
    if N ==100
        dxReal = dx;
    end
>>>>>>> master
    x=(0:dx:L)';
    
<<<<<<< HEAD
    % Sourse volumique de chaleur q[W/m^3] d'ï¿½paisseur dL
    % La source est intï¿½grï¿½e dans la partie intï¿½rieure du mur
    dL=0.10; 
    q=2000;% W/m^3;
=======
    % Sourse volumique de chaleur q[W/m^3] d'épaisseur dL
    % La source est intégrée dans la partie intérieure du mur
>>>>>>> master
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
title('Distribution de Température d''équilibre')
hold

Err=abs(Tmax(1:2:end)-Tmax(2:2:end));

figure(2)
loglog(L./Nar1,Err,'o')
xlabel('dx')
ylabel('Err(dx)=|Tmax(dx)-Tmax(dx/2)|')
title('Erreur sur la température max de l''équilibre')

figure(3)
loglog(Nar,time,':o')
xlabel('N')
ylabel('temps [s]')
title('Temps de calcul en fonction du pas')

Tmax_eq=Tmax(end-1);

%% Méthode dépendante du temps
clc
% Initialisation des paramètres
A = AReal;
dx = dxReal;
N =100;
alpha = Cv*rho/k;
dt = 1* (alpha*dx^2);
x=(0:dx:L)';
S=q.*heaviside(x-(L-dL));
b= -S/k; b(1)=-2*c3/dx; b(N+1)=-2*d3/dx;
t_final = 10e6;

type getUp1.m

P = t_final/dt;

M = diag(ones(1,N+1));
M(1,1) = 0;
M(end,end) = 0;
M = sparse(M);

%Calcul de la matrice U
U = ones(N+1);
U(:,1) = Ta*ones(N+1,1);



for i = 1:P
    U(:,i+1) = getUp1(U(:,i),b,A,M,ksi,dx,dt,alpha);
end


TmaxTau = max(U(:,1)) + 0.99*(Tmax_eq - max(U(:,1)));


MaxU=max(U);


i=1;
while MaxU(i) <= TmaxTau
    i=i+1;    
end
Result=i-1; %(i-1) car on a fait i+1 lorsqu'on avait MaxU(i)=TmaxTau

Tau = Result*dt; %(s)

fprintf('La valeur du pas d''espace est %.3e m \n',dx)
fprintf('La valeur du pas de temps est %.3e s \n',dt)
fprintf('La température d''équilibre est %.3e C \n',Tmax_eq)
fprintf('Erreur sur la température d''équilibre est %.3e C \n',Err)
fprintf('La valeur de temps d''équilibrage est %.3e s \n',Tau)

figure(4)
timeVec = 0:dt:t_final;
plot(timeVec,MaxU)
hold
plot(timeVec,TmaxTau*ones(1,length(timeVec)))
hold
xlabel('t [s]')
ylabel('Tmax(t) [C]')
title('Température maximale en fonction du temps')
legend('Temp. max (t)','Temp. de tau','location','best')



            
        

