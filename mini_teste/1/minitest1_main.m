%main[

%initialisation
k = 20;
g = 9.8;
m = 50e-3;
T=2*pi*sqrt(m/k);
alpha=0.12;
dt= 10.^((-4.75):(0.25):(-3))*T/pi;
y0 = -1e-2;
v0 = -1.5e-2;


%%

y = {};
time={};
dmax=zeros(8,1);
for n = 1:2
    dt = dt/n;
    for i1 = 1:length(dt);
        t=(0:dt(i1):10*T);
        time{i1}=t;
        y{i1}=zeros(length(t),1);
        y{i1}(1)=y0;
        y{i1}(2)=(1-(k*dt(i1)^2)/(2*m))*y0 + dt(i1)*(1-(alpha*dt(i1)/(2*m)))*v0 - (g*dt(i1)^2)/2 ;
        for i2= 3:length(t);
            y{i1}(i2) = calculetY( dt(i1),y{i1}(i2-1),y{i1}(i2-2),k,m,alpha,g);
        end
    end
    
    for i3=1:1:length(y)
       % figure(i3)
        dmax(i3)= max(abs(y{i3}-y0));
       % plot(time{i3}(:),y{i3}(:),'.')
    end
    if n==1
        err = dmax;
    else
        err = abs(err-dmax);
    end
end
