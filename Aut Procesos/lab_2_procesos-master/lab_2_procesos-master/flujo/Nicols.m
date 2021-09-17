clc
clear all
close all 
%%
% data
load('Datos_Modulo_2.mat')
ti = TIME;
fi = FUJO;
pi = PRESION;
ui = BOMBA;
% plot(t,f)
% figure
% plot(t,p)
% figure
% plot(t,u)
%recorte de data
cut_data_i = 170;
cut_data_f = 160;
t = ti(1:end-cut_data_i-cut_data_f);
y = fi(cut_data_i:end-cut_data_f-1);
u = ui(cut_data_i:end-cut_data_f-1);
%%
% identificación
s = tf('s');
vector = [t,y];
% calculo de variables
dy = max(y) - min(y);
du = max(u);
% du = max(u) - min(u);

dy25 = dy*0.25;
dy50 = dy*0.50;
dy75 = dy*0.75;
% interpolar
dt25 = CalcularPunto(vector,dy25)
dt50 = CalcularPunto(vector,dy50)
dt75 = CalcularPunto(vector,dy75)
kp = dy / du;
tao =  0.9102*(dt75 - dt25)
tm = 1.2620*dt25 - 0.2620*dt75
%%
Gp1 = kp*exp(-tm*s)/(tao*s+1)
yi = lsim(Gp1,u,t);
plot(t,yi)
hold on
plot(t,y,'r')
legend('Orden 1','Planta','Location','SouthEast')
%%
% Calculo del indice de aproximacion
% definir un tiempo para ambas señales
data_error = y-yi;
error = sum(((y-yi).^2))
tam_dato = size(t);
error_promedio = error/tam_dato(1)
%%
% Ecuaciones de Sintonizacion
Kc=1.2*tao/(kp*tm);
Ti=2*tm;
Td=0.5*tm;
% PID=Kc+(Ti/s)+Td*s;
PID=Kc*(1+(1/(Ti*s))+((Td*s)/(1+tao*s)));
infoy=stepinfo(feedback(Gp1,1));
Tr=infoy.RiseTime;
T=Tr/15;
GIp=pade(Gp1,1);
GIz=c2d(GIp,T,'tustin');
PIDz=c2d(PID,T,'tustin');
step(du*GIz)
title(['Identificación orden 1 con tiempo muerto Error: ',num2str(error_promedio),' '])
legend('Orden 1 ','Planta','dis Tustin','Location','SouthEast')
xlabel('Tiempo[s]') 
ylabel('Caudal [L/h]')
figure
%%
step(feedback(PID*GIp,1))
legend('PID Nicols','Location','SouthEast')

Tamano=300;
uk = ones(1,Tamano);
yk(1)=0;
yk(2)=0;
yk(3)=0;
ref = 1;
for k = 3: Tamano
%     planta
    yk(k) = ((-1.316e-05)*uk(k))+((6.543e-07)*uk(k-1))+((1.381e-05)*uk(k-2))-((-1.817)*yk(k-1))-(0.8235*yk(k-2));
%     error
    e(k) = ref -yk(k);
%     Controlador
    uk(k) = (3733*e(k))-(7319*e(k-1))+(3588*e(k-2))+(1.864*uk(k-1))-(0.8645*uk(k-2));  
end
figure
plot(yk)



