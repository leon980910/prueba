clear all
close all
clc

s=tf('s')

r = exp(-105.5*s)
[Nr,Dr] = pade(105.5,1)
re = tf(Nr,Dr)
G = 1.48/((13.4*s+1)*(2*s+1))

kp = 1.68363538710378
ki = 0.113189043355807

C = kp*(1+ki/s)

PIz = c2d(C,0.01,'tustin')
rez = c2d(re,0.01,'tustin')
gz =c2d(G,0.01,'tustin')

Cs = (C/(1+C*G))

Ceq = minreal(Cs/(1-Cs*G*re))


[NN,DD] = tfdata(Ceq,'v')

H= ((C*G)/(1+C*G))*re
step(H,'+')
hold on
step(feedback(Ceq*G*re,1))


[N,D] = tfdata(G,'v')
