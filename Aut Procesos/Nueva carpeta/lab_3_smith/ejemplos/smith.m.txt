clear all
close all
clc
t=[1 2 3 4 5 6 7 12 22 32 42];
Tem =[0 0 2.2 3.9 5.3 6.3 7.1 9.2 9.9 9.9 10]
du =1
dy = 10
y25=dy*0.25
y50=dy*0.5
y75=dy*0.75
%interpolacion lineal
t25=interp1([2.2 3.9],[3 4],y25,'linear')
t50=interp1([3.9 5.3],[4 5],y50,'linear')
t75=interp1([7.1 9.2],[7 12],y75,'linear')
k= dy/du
tao = 0.9102*(t75-t25)
l = 1.2620*t25-0.2620*t75
s= tf('s')
P = pade(k* exp(-l*s)/(tao*s+1),1)
Gs = k/(tao*s+1)
[N,D]= pade(l,1)
retardo = tf(N,D)
 kp = 0.018013
 ki = 0.01396
 kd =0.0058106
 pi = kp+ki/s+kd*s
C = feedback(pi,Gs)
ceq = feedback(C,-P)
step(feedback(ceq*P,1))