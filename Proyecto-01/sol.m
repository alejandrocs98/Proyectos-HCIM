% Resolver sistemas de ecuaciones diferenciales ordinarias
global K1 K2 K3 K4 K5 K6 K7 K8   
K1=0.15;
K2=1;
K3=1;
K4=1;
K5=1;
K6=1;
K7=1;
K8=2.5;
[t,X]=ode45(@FUN,[0 30],[0 0 0 0 0.5 0.5]);
plot(t,X)