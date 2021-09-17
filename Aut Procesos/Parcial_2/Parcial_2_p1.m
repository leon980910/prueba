% plantas de diseño
clc
clear all
close all
s = tf('s');
G1 = 0.2/((5*s+1)*(50*s+1));
G2 = 0.002/5*s;
% Os = 0% y ts = 6 segundos

% Diseño del controlador 1
% Criterios de Diseño
ts = 7; 
Mp = 0.01 ; %%overshoot

sigma = 4/ts;
wd = -pi*sigma/(log(Mp));

% Calculo del punto a compensado
% Notas importantes p1 y p2 son los polos de la planta 
% si los polos son menores a 
p1 = 1;
syms a ;
eqn = atan(abs(wd)/(a-sigma)) - (atan(abs(wd)/(-sigma+p1)))  - (pi-atan(abs(wd)/(sigma))) == -pi;
a = double(solve(eqn,a))


% Calculo de la ganacia del controlador
syms kd ;
eqn_1 = kd*(0.5)*(sqrt(wd^2+(a-sigma)^2))/((sqrt(wd^2+sigma^2))*(sqrt(wd^2+(p1-sigma)^2))) ==1 ;
kd = double(solve(eqn_1,kd))

%%estrucutra del PI
s = tf('s');
C1 = kd*(s+a)/s;

step(1*feedback(C1*G1,1));
title('Controlador esclavo')

Gf = (G2)*(feedback(C1*G1,1));
t1  = 0:0.2:90;
y1 = step(Gf,t1);
%%
% Identificacion del sistema en lazo cerrado
s = tf('s');
vector = [t1',y1];
% calculo de variables
dy = max(y1) - min(y1);
du = 1 ;
% du = max(u) - min(u);

dy25 = dy*0.25;
dy50 = dy*0.50;
dy75 = dy*0.75;
% -----------------------------
% interpolar
dt25 = CalcularPunto(vector,dy25)
dt50 = CalcularPunto(vector,dy50)
dt75 = CalcularPunto(vector,dy75)
k = dy / du;
tao =  0.9102*(dt75 - dt25)
tm = 1.2620*dt25 - 0.2620*dt75

Gf1 = k*exp(-tm*s)/(tao*s+1);
Gf3 = k/(tao*s+1);
figure
step(Gf1)
hold on
step(Gf)
title('identificacion')
legend('identificado','Feedback','SouthEast')

%%
% Diseño del controlador 2
% Criterios de Diseño
ts = 90; 
Mp = 0.05 ; %%overshoot

sigma1 = 4/ts;
wd1 = -pi*sigma1/(log(Mp));

% Calculo del punto a compensado
% Notas importantes p1 y p2 son los polos de la planta 
% si los polos son menores a 
p1 = 0.06816;
syms a1 ;
eqn = atan(abs(wd1)/(a1-sigma1)) - (atan(abs(wd1)/(-sigma1+p1)))  - (pi-atan(abs(wd1)/(sigma1))) == -pi;
a1 = double(solve(eqn,a1))


% Calculo de la ganacia del controlador
syms kd1 ;
eqn_1 = kd1*(0.5)*(sqrt(wd^2+(a-sigma)^2))/((sqrt(wd^2+sigma^2))*(sqrt(wd^2+(p1-sigma)^2))) ==1 ;
kd1 = double(solve(eqn_1,kd1))

s = tf('s');
C2 = kd1*(s+a1)/s;

figure
step(80*feedback(C2*Gf1,1))
title('Controlador 2')

s = tf('s');
P_presion = exp(-0.5*s);
P_temperatura = (exp(-5*s))/(1+5*s);