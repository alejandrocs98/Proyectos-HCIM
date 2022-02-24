function dc=Threonine(t,x);
global ATP ADP NADP NADPH K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18 K19 K20 K21 K22 K23 K24 K25 K26 K27 K28 K29 K30 K31 K32 K33 K34 K35 K36 VAKI VAKIII VASD VHDH VHK VTS 
dc=zeros(6,1);

v_AKI = (VAKI*(x(1)*ATP - (x(2)*ADP)/K1))/((K2*((1+(x(6)/K6)^.K7)/)));
v_AKIII = ;
v_ASD = ;
v_HDH = ;
v_HK = ;
v_TS ;

dc(1) = -v_AKI - v_AKIII;
dc(2) = v_AKI + v_AKIII - v_ASD;
dc(3) = v_ASD - v_HDH;
dc(4) = v_HDH - v_HK;
dc(5) = v_HK - v_TS;
dc(6) = v_TS;

end

