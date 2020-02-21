% Distribution de la temp�rature dans un appartement d'un immeuble aux plusieurs �tages


% �quation de transfert de chaleur:
% k*(d^2 T(x,y)/dx^2 + d^2 T(x,y)/dy^2)+S=0


% Conditions aux limites:

% (1) Condition convective (de Robin) � x=0 et � x=Lx (faces externes du mur):
% -k*dT(x=0,y)/dx=-h*(T-Ta)
% -k*dT(x=L,y)/dx=h*(T-Ta)
Ta=-20; %oC

% (2) Condition de Dirichlet sur le plafond et sur le plancher
% T(x, y=0 ou y=Ly)=Tp
Tp=20; %oC

% Dimensions d'appartement
Lx=4; %[m]
Ly=2.4;  %[m]

% Parametres d'un mur d'isolation thermique
Lm=0.4; %m ; �paisseur du mur en brique
km=0.85;%W/(m*K); La conductivit� thermique de la brique
h=1; %W/(m^2*K); Coefficient de transfert thermique sur les surfaces ext�rieures du mur

% Param�tres de l'air qui remplit l'appartement
ka=0.024;

d_ar=[];tini_ar=[];tinv_ar=[];mem_ar=[];Tm_ar=[];
for fact=[2, 1, 1/2, 1/4]
    d=0.1*fact; %Pas de discr�tisation en [m]
    d_ar=[d_ar d]
    Nx=round(Lx/d)+1;
    Ny=round(Ly/d)+1;
    
    tic
    % Initialisation de la source de chaleur et de la conductivit� thermique
    S=zeros(Ny,Nx); k=zeros(Ny,Nx);
    for i=1:Ny
        y=(i-1)*d;
        for j=1:Nx
            x=(j-1)*d;
            
            % Sourse volumique de chaleur q[W/m^3] d'�paisseur dL.
            % La source est int�gr�e dans les parties int�rieures du mur � x=Lm et � x=Lx-Lm et
            % il occupe les tiers du mur dans la direction verticale
            dL=0.1;
            q=1e3;% W/m^3;
            if (x<=Lm)&&(y<=Ly/3+Lm)&&(y>Lm)
                % � l'int�rieur de l'�l�ment chauffant
                S(i,j)=q*exp(-((x-Lm)/dL).^2);
            elseif (x>=(Lx-Lm))&&(y<=Ly/3+Lm)&&(y>Lm)
                % � l'int�rieur de l'�l�ment chauffant
                S(i,j)=q*exp(-((Lx-Lm-x)/dL).^2);
            else
                % � l'ext�rieur de l'�l�ment chauffant
                S(i,j)=0;
            end
            
            % L'espace de vie de l'appartement est d�limit� par
            % les parois d'�paisseur Lm � tous les quatre c�t�s
            if (x<=Lm)||(x>=(Lx-Lm))||(y<=Lm)||(y>=(Ly-Lm))
                % � l'int�rieur du mur
                k(i,j)=km;
            else
                % � l'int�rieurde de l'appartement
                k(i,j)=ka;
            end
        end
    end
    
    M=zeros(Nx*Ny,Nx*Ny);
    b=zeros(Nx*Ny,1);
    
    for i=1:Ny
        for j=1:Nx
            % remplir la ligne pl de la matrice M
            pl=(i-1)*Nx+j;
            
            if ((i>1)&(i<Ny))&((j>1)&(j<Nx))
                % noeud qui est strictement � l'int�rieur de la cellule de simulation
                pc=pl;M(pl,pc)=-4; % contribution de noeud (i,j)
                pc=(i-1)*Nx+j-1;M(pl,pc)=1; % contribution de noeud (i,j-1)
                pc=(i-1)*Nx+j+1;M(pl,pc)=1; % contribution de noeud (i,j+1)
                pc=(i-2)*Nx+j;M(pl,pc)=1; % contribution de noeud (i-1,j)
                pc=(i)*Nx+j;M(pl,pc)=1; % contribution de noeud (i+1,j)
                b(pl)=-d^2*S(i,j)/k(i,j);
            elseif (i==1)
                % noeud sur le plafond y=0
                pc=pl;M(pl,pc)=1; % contribution de noeud (1,j)
                b(pl)=Tp;
            elseif (i==Ny)
                % noeud sur le plancher y=Ly
                pc=pl;M(pl,pc)=1; % contribution de noeud (Ny,j)
                b(pl)=Tp;
            elseif (j==1)
                % noeud � la surface externe du mur x=0
                pc=pl;M(pl,pc)=3+2*d*h/k(i,j); % contribution de noeud (i,1)
                pc=(i-1)*Nx+j+1;M(pl,pc)=-4; % contribution de noeud (i,2)
                pc=(i-1)*Nx+j+2;M(pl,pc)=1; % contribution de noeud (i,3)
                b(pl)=2*d*h*Ta/k(i,j);
            elseif (j==Nx)
                % noeud � la surface externe du mur x=Nx
                pc=pl;M(pl,pc)=3+2*d*h/k(i,j); % contribution de noeud (i,Nx)
                pc=(i-1)*Nx+j-1;M(pl,pc)=-4; % contribution de noeud (i,Nx-1)
                pc=(i-1)*Nx+j-2;M(pl,pc)=1; % contribution de noeud (i,Nx-2)
                b(pl)=2*d*h*Ta/k(i,j);
            else
                display('Erreur dans la d�finition de la matrice de coefficients');
            end
        end
    end
    tini_ar=[tini_ar toc]
    
    tic
    %T=M\b;
    [L,U]=lu(M);T=U\(L\b);
    tinv_ar=[tinv_ar toc]
    
    mem_ar=[mem_ar 8*(Nx*Ny)^2]
    
    Tr=reshape(T,Nx,Ny)';
    
    Tm_ar=[Tm_ar Tr(round(Ly/d/2+1),round(Lx/d/2+1))]
    
end

figure(1)
h=pcolor((0:d:Lx),(0:d:Ly),S);set(h,'LineStyle','none')
colorbar
xlabel('x [m]'); ylabel('y [m]'); title('S(x,y) [W/m^3]')
axis equal
axis tight

figure(2)
h=pcolor((0:d:Lx),(0:d:Ly),k);set(h,'LineStyle','none')
colorbar
xlabel('x [m]'); ylabel('y [m]'); title('k(x,y) [W/(m^2\cdotK)]')
axis equal
axis tight

figure(3)
h=pcolor((0:d:Lx),(0:d:Ly),Tr);set(h,'LineStyle','none')
colorbar
xlabel('x [m]'); ylabel('y [m]'); title('T(x,y) [^oC]')
axis equal
axis tight

figure(4)
loglog(d_ar,mem_ar/1024^3,':o')
axis([1e-3 1 1e-4 128])
axis square
xlabel('\delta [m]')
ylabel('mem [Gb]')

figure(5)
loglog(d_ar(2:1:end),abs(Tm_ar(1:end-1)-Tm_ar(2:1:end)),':o')
axis tight
xlabel('\delta [m]')
ylabel('Err [^oC]')

figure(6)
loglog(d_ar,tini_ar,':o')
axis([1e-3 1 1e-3 10000])
axis square
xlabel('\delta [m]')
ylabel('t_{ini} [s]')

figure(7)
loglog(d_ar,tinv_ar,':o')
axis([1e-3 1 1e-3 10000])
axis square
xlabel('\delta [m]')
ylabel('t_{inv} [s]')
