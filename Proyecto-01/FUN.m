function px=FUN(t,x);
global K1 K2 K3 K4 K5 K6 K7 K8  
px=zeros(6,1);
px(1)=K1-K2*x(1)*x(5);
px(2)=K2*x(1)*x(5)-K3*x(2)+K4*x(3);
px(3)=K3*x(2)-K4*x(3)-K5*x(3)*x(5)+K6*x(4);
px(4)=K5*x(3)*x(5)-K6*x(4)-K7*x(4);
px(5)=-K2*x(1)-K5*x(3)+K8*x(6);
px(6)=K2*x(1)+K5*x(3)-K8*x(6);



end

