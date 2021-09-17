%%Diseno de controlador  PD
clc
clear all
close all
format long
s = tf('s');
%%Planta de ejemplo
G = 1080/(s*(s+6)*(s+18));
[numt,dent]=tfdata(G, 'v');
% Criterios de Diseño
ts = 0.2; 
Mp = 10 ; %%overshoot

sigma = 4/ts;
wd = -pi*sigma/(log(Mp));

% Calculo del punto a compensado
% Notas importantes p1 y p2 son los polos de la planta 
% si los polos son menores a 
p1 = 18;
p2 = 6;
syms a ;
eqn = atan(abs(wd)/(a-sigma)) - atan(abs(wd)/(p1-sigma))-pi -atan(abs(wd)/(p2-sigma))-pi +atan(abs(wd)/(sigma)) == -pi;
a = double(solve(eqn,a))


% Calculo de la ganacia del controlador
syms kd ;
eqn_1 = kd*(1080)*(sqrt(wd^2+(a-sigma)^2))/((sqrt(wd^2+sigma^2))*(sqrt(wd^2+(p1-sigma)^2))*(sqrt(wd^2+(p2-sigma)^2))) ==1 ;
kd = double(solve(eqn_1,kd))

%%estrucutra del PD
s = tf('s');
PD = kd*(s+a);

step(feedback(PD*G,1))