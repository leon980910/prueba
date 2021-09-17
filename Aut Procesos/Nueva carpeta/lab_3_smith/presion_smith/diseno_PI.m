%Diseño PI de presion
%%Diseno de controlador  PD
clc
clear all
close all


%%Planta de ejemplo

% Criterios de Diseño
ts = 19 ; 
Mp = 5 ; %%overshoot

sigma = 4/ts;
wd = -pi*sigma/(log(Mp));

% Calculo del punto a compensado
% Notas importantes p1 y p2 son los polos de la planta 
% si los polos son menores a 
p1 = 0.4;
syms a ;
eqn = atan(abs(wd)/(a-sigma)) - atan(abs(wd)/(p1-sigma))  - (pi-atan(abs(wd)/(sigma))) == -pi;
a = double(solve(eqn,a))


% Calculo de la ganacia del controlador
syms kd ;
eqn_1 = kd*(0.0001062)*(sqrt(wd^2+(a-sigma)^2))/((sqrt(wd^2+sigma^2))*(sqrt(wd^2+(p1-sigma)^2))) ==1 ;
kd = double(solve(eqn_1,kd))

%%estrucutra del PI
s = tf('s');
PI = kd*(s+a)/s;

load('planta.mat')% G,Gp1,tm
load('controller.mat')%Controlador sintonizado por control desing
step(2*feedback(PI*G,1))