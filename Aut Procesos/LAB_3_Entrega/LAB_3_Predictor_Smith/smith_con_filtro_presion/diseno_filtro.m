clc
clear all
close all

load('datos_mod_2_completos.mat')
u = PRESION(60:200);
%filtro pasabajos
t = 0.6;
t = 0.2;
wn = 1/t;
fm = wn/2;
[B,A]=butter(2,0.83,'low');%filtramos se�al de entada
% [B,A]=butter(2,0.66,'low');%filtramos se�al de entada
% [B,A]=butter(2,0.2,'low');%filtramos se�al de entada filtro a probar
H1=filter(B,A,u);
figure
plot(H1)
figure
plot(H1)
grid on
hold on 
plot(u)
title(['Filtro Pasabajas'])
legend('Se�al Filtrada','Se�al Fuente','Location','SouthEast')
xlabel('Tiempo [s]') 
ylabel('Presi�n [bar]')