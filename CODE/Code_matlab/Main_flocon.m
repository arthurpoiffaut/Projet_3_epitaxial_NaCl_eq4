% MAIN

%% Initialisation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%init variable 
dimx=50; %dimention x  ect
dimy=50;
dimz=5;

%init des vecteur pour regarder les plus proche voisin
dx = [-1,0,-1,1,0,1,0,0];
dy = [-1,-1,0,0,1,1,0,0];
dz = [0,0,0,0,0,0,1,-1];


% val
delta_dim=10^(-6); % valeur de la dimention 1 micron
delta_t=10^(-9); % valeur de la variation de temps
D=10^(-5); %valeur du coe de diff   m^2/sece
Vcell=sqrt(3/2)*(delta_dim^3);
nu_kin=133; % micro metre/sec
delta_tau=(D*delta_t)/(delta_dim)^2;  %(plus peti que 1/6)

%voir poure les autre truc 


%autre para




%def mat 
% mat cristal 0 si pas dans le cristal
ice=zeros(dimx,dimy,dimz);
% mat vapeur varie enfonction de la densiter deau en éta vapeur entre 0 et 1 ?
vap=zeros(dimx,dimy,dimz);
%mat de la frontiere
fron=zeros(dimx,dimy,dimz);

%test glace



%teste pour difusion

ice(1,10,:)=1;
ice(2,10,:)=1;
ice(3,10,:)=1;
ice(4,10,:)=1;
ice(5,10,:)=1;
ice(6,10,:)=1;
ice(7,10,:)=1;
ice(8,10,:)=1;
ice(9,10,:)=1;
ice(10,10,:)=1;
vap(:,:,:)=0.03;
vap(1,1,:)=1;
vap(end,end,:)=1;

for i1=1:1:5000;
    vap= diffusion_vap(vap,ice,delta_tau,dx,dy,dz,dimx,dimy,dimz);
end
vapout2d=sum(vap,3)
imagesc(vapout2d)
