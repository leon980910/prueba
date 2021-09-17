clear all 
close all 
clc 

load('datos_02.mat')

F=FUJO(60:end,1);
B=BOMBA(60:963,1);
P=Perturbacion(60:963,1);
pre=PRESION(60:end,1);
T=tout(60:963,1)-5.8;


figure (1)
plot(T,F)

figure (2)
plot(T,P)

figure (3)
plot(T,B)
figure (4)
plot(T,pre)

