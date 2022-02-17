% Análisis de estabilidad
global K1 K2 K3 K4 K5 K6 K7 K8 
K1=0.15;
K2=1;
K3=1;
K4=1;
K5=1;
K6=1;
K7=1;
K8=2.5;

%matriz analitica
syms x1 x2 x3 x4 x5 x6 
f=[K1-K2*x1*x5;K2*x1*x5-K3*x2+K4*x3;K3*x2-K4*x3-K5*x3*x5+K6*x4;K5*x3*x5-K6*x4-K7*x4;-K2*x1-K5*x3+K8*x6;K2*x1+K5*x3-K8*x6];	
v=[x1,x2,x3,x4,x5,x6];
tr=jacobian(f,v) 

%evaluacion en valores de estado estacionario
x1=0.5202;
x2=1.1955;
x3=1.0454;
x4=0.1506;
x5=0.2880;
x6=0.6262;
m=eval(tr)
%valores propios
jc=eig(m)