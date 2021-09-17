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
a = (-0.6240*dt25 + 0.9866*dt50-0.3626*dt75)/(0.3533*dt25-0.7036*dt50 + 0.3503*dt75);
tao2 =   (dt75-dt25)/(0.9866 + 0.7036*a);%???
t1 = tao2;
t2 = tao2*a;
tm2 = dt75-(1.3421 + 1.3455*a)*tao2;


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