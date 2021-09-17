clc; clear all; close all;
%% Parcial T-2C
s=tf('s');%S es ahora de una TF
Gps=10/(5*s+6);
Gp=Gps*exp(-0.5*s);
%% Condiciones de diseño
ts=3;
Mp=10;
sigma=4/ts;
wd=(-pi*sigma)/log(Mp/100);
%% rlocus
rlocus(Gps/s)
hold on
plot(-sigma,wd,'hm')
hold off
%% Condicion de angulo 
syms x;
polos=pole(Gps/s);
ceros=zero(Gps/s);
Spolos=(180-atand(wd/(abs(polos(1))-sigma)))+(180-atand(wd/(abs(polos(2))-sigma)));
sA=-Spolos;
cA=atand(wd/(x-sigma))+sA;
% 180*(2*landa+1)/n-m%Para sabe a que igualar
a=double(solve(cA==-450));
%% Condición de Magnitud
%10/6=1.667,6=>6(5/6s+1)
cM=x*1.6667*((sqrt(wd^2+(a-sigma)^2))/(sqrt(wd^2+(abs(polos(1))-sigma)^2)*sqrt(wd^2+(abs(polos(2))-sigma)^2)));
kd=double(solve(cM==1));
%% a
Gc=kd*(s+a)/s;
%% Smith
tm=0.5;
kp=10/6;
tau=5/6;
[Num,Den]=pade(tm,1);
ret=tf(Num,Den);
Gr=(kp*ret)/(tau*s+1);
CL=feedback(Gc,Gps);
%% b
ceq=feedback(CL,-Gr);
%% Discretización
infoy=stepinfo(feedback(Gp,1));
% infoy=stepinfo(Gp);
Tr=infoy.RiseTime;
T=Tr/15;
Gpz = c2d(Gr,T,'tustin');
Gcz=c2d(ceq,T,'tustin');
%% Diferencias
[NGz,DGz]=tfdata(Gpz,'v');
[NGcz,DGcz]=tfdata(Gcz,'v');
Tamano=1000;
uk=ones(1,Tamano);
yk(1)=0;
yk(2)=0;
yk(3)=0;
yk(4)=0;
yk(5)=0;
ref = 20;
for k = 5: Tamano
e(k) = ref-yk(k-1);
uk(k) = NGcz(1)*e(k)+NGcz(2)*e(k-1)+NGcz(3)*e(k-2)+NGcz(4)*e(k-3)+NGcz(5)*e(k-4)-...
        DGcz(2)*uk(k-1)-DGcz(3)*uk(k-2)-DGcz(4)*uk(k-3)-DGcz(5)*uk(k-4);
yk(k) = NGz(1)*uk(k)+NGz(2)*uk(k-1)+NGz(3)*uk(k-2)-...
        DGz(2)*yk(k-1)-DGz(3)*yk(k-2);  

% yk(k) = NGz(1)*uk(k)+NGz(2)*uk(k-1)+NGz(3)*uk(k-2)-...
%         DGz(2)*yk(k-1)-DGz(3)*yk(k-2);  
% e(k) = ref-yk(k);
% uk(k) = NGcz(1)*e(k)+NGcz(2)*e(k-1)+NGcz(3)*e(k-2)+NGcz(4)*e(k-3)+NGcz(5)*e(k-4)-...
%         DGcz(2)*uk(k-1)-DGcz(3)*uk(k-2)-DGcz(4)*uk(k-3)-DGcz(5)*uk(k-4);

% yk(k) = NGz(1)*uk(k)+NGz(2)*uk(k-1)+NGz(3)*uk(k-2)-...
%         DGz(2)*yk(k-1)-DGz(3)*yk(k-2);  
% e(k) = ref-yk(k);
% uk(k) = NGcz(1)*e(k)+NGcz(2)*e(k-1)-...
%         DGcz(2)*uk(k-1);
end
%% Inciso C
clear s
syms s
E=kp*exp(-tm*s)/(tau*s+1);
Esc=double(limit((s/(1+E))*(1/s)))
P=(4.304*s^4+39.95*s^3+126.9*s^2+161.8*s+71.47)/(4.137*s^4+43.88*s^3+116.3*s^2+83.56*s);
Ecc=double(limit((s/(1+(E*P)))*(1/s^2)))
%% Gráficas
step(ref*Gps,ref*feedback(Gc*Gps,1),ref*feedback(ceq*Gr,1),ref*feedback(Gcz*Gpz,1))
title('Parcial Teorico')
legend('Planta','PI','Smith','Discreto')
figure
plot(yk,'k')
title('Parcial Teorico')
legend('Diferencias')