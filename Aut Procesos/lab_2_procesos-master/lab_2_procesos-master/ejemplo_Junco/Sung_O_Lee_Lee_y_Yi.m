%% Con datos experimentales - Flujo
clc; clear all; close all;
load('Datos.mat')
s=tf('s');%S es ahora de una TF
tf = Data(1:158,1);
y=Data(4:161,3);
Sy=size(y);
u=Data(4:161,4);
dy=max(y)-min(y);
du=max(u)-min(u);
dy25=dy*0.25;
dy50=dy*0.5;
dy75=dy*0.75;
vdy=[dy25; dy50; dy75];
vdt= [Funcion(y,tf,dy25) Funcion(y,tf,dy50) Funcion(y,tf,dy75)]';
% disp([vdy_p vdt_p]);
%Modelo de primer orden más tiempo muerto (POMTM)
kp=dy/du;
tau=0.9102*(vdt(3)-vdt(1));
tm=1.2620*vdt(1)-0.2620*vdt(3);
GI=(kp*exp(-tm*s))/(tau*s+1);
[yi,ti]=step(20000*GI,tf);
error = sum(((y-yi).^2));
error_m=error/Sy(1);
% Modelo de polo doble más tiempo muerto (PDMTM)
tau2=0.5776*(vdt(3)-vdt(1));
tm2=1.552*vdt(1)-0.5552*vdt(3);
% Modelo de segundo orden más tiempo muerto (SOMTM)
%Método simplificado (SOMTMs)
tm3=tm2;
alpha=(vdt(2)-tm2-1.4362*tau2)/(1.9844*tau2-vdt(2)+tm2);
tau3=(2*tau2)/(1+alpha);
tau31=tau3;
tau32=alpha*tau3;
GI3=(kp*exp(-tm3*s))/((tau3*s+1)*(alpha*tau3*s+1));
[yi3,ti]=step(20000*GI3,tf);
error3 = sum(((y-yi3).^2));
error_m3=error3/Sy(1);
infoy=stepinfo(feedback(GI3,1));
Tr3=infoy.RiseTime;
T3=Tr3/15;
GIp3=pade(GI3,1);
GIz3=c2d(GIp3,T3,'tustin');
clc; clear all; close all;
load('Datos.mat')
s=tf('s');%S es ahora de una TF
tf = Data(1:158,1);
y=Data(4:161,3);
Sy=size(y);
u=Data(4:161,4);
dy=max(y)-min(y);
du=max(u)-min(u);
dy25=dy*0.25;
dy50=dy*0.5;
dy75=dy*0.75;
vdy=[dy25; dy50; dy75];
vdt= [Funcion(y,tf,dy25) Funcion(y,tf,dy50) Funcion(y,tf,dy75)]';
% disp([vdy_p vdt_p]);
kp=dy/du;
tau=0.9102*(vdt(3)-vdt(1));
tm=1.2620*vdt(1)-0.2620*vdt(3);
GI=(kp*exp(-tm*s))/(tau*s+1);
[yi,ti]=step(20000*GI,tf);
error = sum(((y-yi).^2));
error_m=error/Sy(1);
% Modelo de polo doble más tiempo muerto (PDMTM)
% kp=dy/du;
% tau=0.5776*(vdt(3)-vdt(1));
% tm=1.552*vdt(1)-0.5552*vdt(3);
% GI=(kp*exp(-tm*s))/(tau*s+1)^2;
% [yi,ti]=step(20000*GI,tf);
% error = sum(((y-yi).^2));
% error_m=error/Sy(1);
%Ecuaciones de sintonización
%% Tablas
a=[1.435 1.357 1.495];
b=[-0.921 -0.947 -0.945];
c=[0.878 0.842 1.101];
d=[-0.749 -0.738 -0.771];
e=[0.482 0.381 0.560];
f=[1.137 0.995 1.006];
%% 1
Kc=a(1)*((tm/tau)^b(1))/kp;
Ti=(1/c(1))*((tm/tau)^(-d(1)))*tau;
Td=e(1)*((tm/tau)^(f(1)))*tau;
% PID=Kc+(Ti/s)+Td*s;
PID=Kc*(1+(1/(Ti*s))+((Td*s)/(1+tau*s)));
infoy=stepinfo(feedback(GI,1));
Tr=infoy.RiseTime;
T=Tr/15;
GIp=pade(GI,1);
GIz=c2d(GIp,T,'tustin');
PIDz=c2d(PID,T,'tustin');
step(20000*GIz)
Tamano=300;
uk = ones(1,Tamano);
yk(1)=0;
yk(2)=0;
yk(3)=0;
ref = 1;
for k = 3: Tamano
    yk(k) = ((-1.316e-05)*uk(k))+((6.543e-07)*uk(k-1))+((1.381e-05)*uk(k-2))-((-1.817)*yk(k-1))-(0.8235*yk(k-2));
    e1(k) = ref -yk(k);
    uk(k) = (6033*e1(k))-((1.184e04)*e1(k-1))+(5808*e1(k-2))+(1.864*uk(k-1))-(0.8645*uk(k-2));   
end
TiempoMin=Tamano*T/60;
figure
plot(yk)
%% 2
Kc2=a(2)*((tm/tau)^b(2))/kp;
Ti2=(1/c(2))*((tm/tau)^(-d(2)))*tau;
Td2=e(2)*((tm/tau)^(f(2)))*tau;
% PID=Kc+(Ti/s)+Td*s;
PID2=Kc2*(1+(1/(Ti2*s))+((Td2*s)/(1+tau*s)));
% infoy=stepinfo(feedback(GI,1));
% Tr=infoy.RiseTime;
% T=Tr/15;
% GIp=pade(GI,1);
% GIz=c2d(GIp,T,'tustin');
PIDz2=c2d(PID2,T,'tustin');
% step(20000*GIz)
% Tamano=1000;
uk2 = ones(1,Tamano);
yk2(1)=0;
yk2(2)=0;
yk2(3)=0;
% ref = 1;
for k = 3: Tamano
    yk2(k) = ((-1.316e-05)*uk2(k))+((6.543e-07)*uk2(k-1))+((1.381e-05)*uk2(k-2))-((-1.817)*yk2(k-1))-(0.8235*yk2(k-2));
    e2(k) = ref -yk2(k);
    uk2(k) = (3818*e2(k))-((7429)*e2(k-1))+(3616*e2(k-2))+(1.864*uk2(k-1))-(0.8645*uk2(k-2));   
end
% TiempoMin2=Tamano*T/60;
hold on
plot(yk2)
%% 3
Kc3=a(3)*((tm/tau)^b(3))/kp;
Ti3=(1/c(3))*((tm/tau)^(-d(3)))*tau;
Td3=e(3)*((tm/tau)^(f(3)))*tau;
% PID=Kc+(Ti/s)+Td*s;
PID3=Kc3*(1+(1/(Ti3*s))+((Td3*s)/(1+tau*s)));
% infoy=stepinfo(feedback(GI,1));
% Tr=infoy.RiseTime;
% T=Tr/15;
% GIp=pade(GI,1);
% GIz=c2d(GIp,T,'tustin');
PIDz3=c2d(PID3,T,'tustin');
% step(20000*GIz)
% Tamano=1000;
uk3 = ones(1,Tamano);
yk3(1)=0;
yk3(2)=0;
yk3(3)=0;
% ref = 1;
for k = 3: Tamano
    yk3(k) = ((-1.316e-05)*uk3(k))+((6.543e-07)*uk3(k-1))+((1.381e-05)*uk3(k-2))-((-1.817)*yk3(k-1))-(0.8235*yk3(k-2));
    e3(k) = ref -yk3(k);
    uk3(k) = (5649*e3(k))-((1.106e04)*e3(k-1))+(5415*e3(k-2))+(1.864*uk3(k-1))-(0.8645*uk3(k-2));   
end
% TiempoMin2=Tamano*T/60;
% hold on
plot(yk3)
hold off
title (['Comparaci�n PID Ideal Lopez, Miller, Smith y Murril'])
legend('IAE','ITAE','ISE')