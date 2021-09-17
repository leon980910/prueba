clc
clear all
close all 
% Ejemplo del paper
% Funcion original
s = tf('s');
t = 0:0.18:180;
Gpo = (1.25*exp(-0.25*s))/((16*s + 1)*(4*s + 1)*(2*s + 1)*(s + 1))
[y,t] = step(Gpo,t);
vector = [t,y];
step(Gpo);
% calculo de variables
dy = max(y) - min(y);
du = 1;
dy25 = dy*0.25;
dy50 = dy*0.50;
dy75 = dy*0.75;
% -----------------------------
% interpolar
dt25 = CalcularPunto(vector,dy25)
dt50 = CalcularPunto(vector,dy50)
dt75 = CalcularPunto(vector,dy75)
kp = dy / du;
tao =  0.9102*(dt75 - dt25)
tm = 1.2620*dt25 - 0.2620*dt75
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
Gp1 = kp*exp(-tm*s)/(tao*s+1)
% Gp1 = (1.25*exp(-6.69*s))/(16.94*s + 1)

step(Gp1)
hold on
step(Gpo)
% ------------------------------------------------------------------------
% Calculo del indice de aproximacion
% definir un tiempo para ambas señales
tf = 0:0.18:180;
% Real
[yt,tf] = step(Gpo,tf);
% identificada
[yi,tf] = step(Gp1,tf);
error = sum(((yt-yi).^2))


