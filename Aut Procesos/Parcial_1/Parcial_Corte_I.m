%%Parcial de Procesos
% Datos
clc
clear all
close all
t = [1 2 3 4 5 6 7 12 22 32 42]';
t = t*60;
y = [0 0 2.2 3.9 5.3 6.3 7.1 9.2 9.9 9.9 10]';
s = tf('s');
vector = [t,y];
dy = max(y)-min(y);
du = 60;
u  = [du du du du du du du du du du du]';
dy25 = dy*0.25+min(y);
dy50 = dy*0.50+min(y);
dy75 = dy*0.75+min(y);
dt25 = CalcularPunto(vector,dy25)
dt50 = CalcularPunto(vector,dy50)
dt75 = CalcularPunto(vector,dy75)
% -----------------------------
% Modelo de primer orden más tiempo muerto (POMTM)
% interpolación
kp = dy / du;
tao =  0.9102*(dt75 - dt25);
tm =  1.2620*dt25 - 0.2620*dt75
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%Planta Identificada
Gp1 = (kp*exp(-tm*s))/(tao*s+1)

%%
% Ecuaciones de Sintonizacion
Kc=1.2*tao/(kp*tm);
Ti=2*tm;
Td=0.5*tm;

PID=Kc*(1+(1/(Ti*s))+((Td*s)/(1+tao*s)));

%%
% Calculo Error
clear s
syms s
%Error sin control
C=pade(0.1667*exp(-116*s)/(260.8*s+1));
Eop=limit((s/(1+C))*(1/s))

%Error con control
P=(1.197e06*s^2+7996*s+16.26)/(6.026e04*s^2+231*s);
Ecl=limit((s/(1+(C*P)))*(1/s^2))


%%Punto 2
% Calculo para el transmisor a 10 mA
x1 = 4;
y1 = 20;
x2 = 20;
y2 = 180;
m = (y2-y1)/(x2-x1);
b = y1 - m*x1;

x=10;
y_10mA = m*x + b;

x_p = [x1 x x2];
y_p = [y1 y_10mA y2];
plot(x_p,y_p)
hold on 
plot(x1,y1,'ro')
hold on 
plot(x2,y2,'ro')
hold on 
plot(x,y_10mA,'ro')
grid on 
title('Curva Caracteristica y =' +string(m) + 'x ' + string(b))
xlabel('Corriente [mA]')
ylabel('Temperatura [°C]')
