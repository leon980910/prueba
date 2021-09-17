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
tao1 = 0.5776*(dt75-dt25);%t'
tm1 = 1.5552*dt25-0.5552*dt75;
tm2 = tm1;
a = (dt50-tm1-1.4362*tao1)/(1.9844*tao1-dt50+tm1);
tao2 = 2*tao1/(1+a);%t´´
t1 = tao2;
t2 = tao2*a;
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
Gp1 = (kp*exp(-tm2*s))/((tao2*s + 1)*(a*tao2*s + 1)) ;



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