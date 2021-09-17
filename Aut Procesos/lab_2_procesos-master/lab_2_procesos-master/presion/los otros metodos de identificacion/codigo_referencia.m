
clear all
close all
clc
format long
tao=3
k=2
s=tf('s')

Gs = k/(tao * s +1)
[n,d] = pade(2,1)
Gt = tf(n,d)*Gs
Gz = c2d(Gt,7/10,'tustin')
step(Gz)
hold on

step(Gt)
[y,tt]=step(Gz,200);

u = ones(1,length(tt))
t(1)=-0.1006
t(2) =-0.1203
for k=3:length(tt)
    
    t(k) = -0.1006*u(k)+0.1083*u(k-1)+0.209*u(k-2) +1.273*t(k-1)-0.3809*t(k-2)
end
hold on
plot(tt,t,'r+')
t1=-0.1006
t2 =-0.1203
u0=1
u1=1
u2=1
for k=3:length(tt)
    
    t0 = -0.1006*u0+0.1083*u1+0.209*u2 +1.273*t1-0.3809*t2
    plot(tt(k),t0,'go')
    t2=t1
    t1=t0
    u2=u1
    u1=u0
    
end
hold on
plot(tt,t,'r+')   