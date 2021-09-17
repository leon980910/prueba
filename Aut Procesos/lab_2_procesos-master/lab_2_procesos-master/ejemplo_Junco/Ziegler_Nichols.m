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
Kc=1.2*tau/(kp*tm);
Ti=2*tm;
Td=0.5*tm;
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
% ref = 1;
ref = 150;
for k = 3: Tamano
    yk(k) = ((-1.316e-05)*uk(k))+((6.543e-07)*uk(k-1))+((1.381e-05)*uk(k-2))-((-1.817)*yk(k-1))-(0.8235*yk(k-2));
    e(k) = ref -yk(k);
    uk(k) = (3733*e(k))-(7319*e(k-1))+(3588*e(k-2))+(1.864*uk(k-1))-(0.8645*uk(k-2));   
end
TiempoMin=Tamano*T/60;
figure
plot(yk)