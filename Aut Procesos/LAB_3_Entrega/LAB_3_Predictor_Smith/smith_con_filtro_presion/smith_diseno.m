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

% step(2*feedback(PI*G,1))

load('data_lineal.mat') 
l = 1.2620*dt25-0.2620*dt75
s= tf('s')
P = pade(k* exp(-l*s)/(tao*s+1),1)
Gs = k/(tao*s+1);
[N,D]= pade(l,1);
retardo = tf(N,D);
% parametros del PID tool

pi = PI;
C = feedback(pi,Gs)
ceq = feedback(C,-P)
step(feedback(ceq*P,1))


%Discretizar
T = 0.6;
ceq_d = c2d(ceq,T)
P_d = c2d(P,T)
step(feedback(ceq_d*P_d,1))
%%ecuaciones de diferencia
[NGz,DGz]=tfdata(P_d,'v')
[NGcz,DGcz]=tfdata(ceq_d,'v')

U=ones(1,100);
% U(1)=0; U(2)=0; U(3)=0; U(4)=0; U(5)=0; U(6)=0; U(7)=0; U(8)=0; U(9)=0; U(10)=0;


U0=0; U1=0; U2=0; U3=0; U4=0; U5=0; U6=0; U7=0; U8=0; 
Y1=0; Y2=0; Y3=0; Y4=0; Y5=0; Y6=0;
E1=0; E2=0; E3=0; E4=0; E5=0; E6=0; E7=0; E8=0;
ref=2;
for k = 9:1000
%     Y(k) = NGz(2)*U(k)+NGz(3)*U(k-1)-...
%         DGz(2)*Y(k-1)-DGz(3)*Y(k-2);
%     E(k) = ref-Y(k);
%     U(k) = NGcz(1)*E(k)+NGcz(2)*E(k-1)+NGcz(3)*E(k-2)+NGcz(4)*E(k-3)+NGcz(5)*E(k-4)-...
%         DGcz(2)*U(k-1)-DGcz(3)*U(k-2)-DGcz(4)*U(k-3)-DGcz(5)*U(k-4);
    
    Y0 = NGz(2)*U0+NGz(3)*U1-...
        DGz(2)*Y1-DGz(3)*Y2;
    E0 = ref-Y0;
    U0 = NGcz(1)*E0+NGcz(2)*E1+NGcz(3)*E2+NGcz(4)*E3+NGcz(5)*E4-...
        DGcz(2)*U1-DGcz(3)*U2-DGcz(4)*U3-DGcz(5)*U4;
    
    salida(k)=Y0;
    salida_c(k) = U0;
    U4=U3; U3=U2; U2=U1; U1=U0;
    Y2=Y1; Y1=Y0;
    E4=E3; E3=E2; E2=E1; E1=E0;   
end
figure
plot(salida)
title(['Control Smith Ecuaciones de Diferencia'])
legend('Presion','Location','SouthEast')
xlabel('iteraciones') 
ylabel('Bar')


