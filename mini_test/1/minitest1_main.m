type calculetY
type calculetAllY

%% question iii
clear all
clc
%initialisation

y0 = -3e-2; %-1e-2;
v0 = -1e-2; %-1.5e-2;
alpha = 0.15; %0.12;

k = 20;
g = 9.8;
m = 50e-3;
T = 2*pi*sqrt(m/k);
dt= transpose(10.^((-4.75):(0.25):(-3))*T/pi);


y = {};
time={};
dmax=zeros(8,1);
timeMultiplier = 0.6;
    
for n = 1:2 %Itère pour dt et dt/2
    dt = dt/n;
    
    [y,time] = calculetAllY(y,time,y0,v0,dt,T,k,m,alpha,g,timeMultiplier);
    
    for i3=1:1:length(y)
<<<<<<< HEAD
        dmax(i3)= max(abs(y{i3}-y0));
        
%         figure(i3)
%         plot(time{i3}(:),y{i3}(:),'.')
%         titl = sprintf('Pas de temps dt = %d',dt(i3));
%         title(titl);
=======
       % figure(i3)
        dmax(i3)= max(abs(y{i3}-y0));
       % plot(time{i3}(:),y{i3}(:),'.')
>>>>>>> master
    end
    if n==1
        dmax1 = dmax;
    else
        err = abs(dmax1-dmax);
    end
end

plot(dt,err)
title('erreur en fonction du pas de temps')
xlabel('Pas de temps [s]')
ylabel('Erreur [cm]')

% table(dt*2,dmax1,err)
fprintf('Pas de temps |Distance max(cm) |Erreur(cm) \n')
for i = 1:length(dt)
    fprintf('%d |%.9f      |%d \n',dt(i)*2,dmax1(i)*1e2,err(i)*1e2)
end

fprintf('L''erreur minimale est 5.229393e-10 cm \n')
fprintf('La distance maximale est 0.986030307 cm \n')
fprintf('Le pas de temps correspondant est: 3.162278e-05 \n')

%% Question iv

dt = 1e-3 *T/pi;
yeq = -m*g/k;
timeMultiplier = 10;
y = {};
time={};
[y,time] = calculetAllY(y,time,y0,v0,dt,T,k,m,alpha,g,timeMultiplier);

cut_function = (y{1}-yeq) .*( (y{1}-yeq)> max((y{1}-yeq)/10)  );
N_osc = length(findpeaks(cut_function));
figure()
plot(time{1}(:),y{1}(:),'.')
title('Trajection de la ball en fonction du temps')
xlabel('temps [s]')
ylabel('Trajectoire [cm]')

figure()
plot(time{1}(:),cut_function,'.')
title('Maximum d''oscillation plus grand que la valeur minimale')
xlabel('temps [s]')
ylabel('Trajectoire [cm]')
fprintf('Le nombre d''oscillation est %d \n',N_osc)
