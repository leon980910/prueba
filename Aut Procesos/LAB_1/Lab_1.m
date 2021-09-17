clear all
close all
clc

s = tf('s') %asociar la s en funcion de transferencia

Gp = 1.25*exp(-0.25*s)/((16*s+1)*(4*s+1)*(2*s+1)*(s+1))

t = Data(:,1)

