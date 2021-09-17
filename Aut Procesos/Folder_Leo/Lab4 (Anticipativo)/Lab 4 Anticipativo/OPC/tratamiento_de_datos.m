clear all 
close all 
clc 

load('datos_01.mat')

F=FUJO(34:1000,1);
B=BOMBA(34:1000,1);
P=Perturbacion(34:1000,1);
T=tout(34:1000,1);


figure (1)
plot(T,F)

figure (2)
plot(T,P)

figure (3)
plot(T,B)

