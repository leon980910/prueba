%%Modelo de Smith Flujo Juan Velandia - Harold Leon
clear all
close all
load('data_lineal.mat') 
l = 1.2620*dt25-0.2620*dt75
s= tf('s')
P = pade(k* exp(-l*s)/(tao*s+1),1)
Gs = k/(tao*s+1)
[N,D]= pade(l,1)
retardo = tf(N,D)
% parametros del PID tool
 kp = 13522.88
 ki = 9748
 kd =0
 pi = kp+ki/s+kd*s
C = feedback(pi,Gs)
ceq = feedback(C,-P)
step(feedback(ceq*P,1))