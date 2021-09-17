clc 
clear all
close all

load('Datos.mat')
s = tf('s');
T = Data(:,1); %Temperatura
F = Data(:,2); %Flujo
P = Data(:,3); %Presion
u = Data(:,4); %Entrada

%% Tratamiento de Datos
%Al inicio
T(1:34)=[];
F(1:34)=[];
P(1:34)=[];
u(1:34)=[];
%Al final
T(99:133)=[];
F(99:133)=[];
P(99:133)=[];
u(99:133)=[];

for i=1:size(T,1)
    T(i)=T(i)-10;
end
%% Graficas de respuesta en el tiempo
figure
hold on
subplot(1,3,1)
plot(T,u)
ylabel('Entrada (Pulsos)')
xlabel('Tiempo (s)')
title('Tiempo vs Entrada')

subplot(1,3,2)
plot(T,F,'r')
xlabel('Tiempo (s)')
ylabel('Flujo (l/h)')
title('Flujo vs Tiempo')

subplot(1,3,3)
plot(T,P,'g')
xlabel('Tiempo (s)')
ylabel('Presion (Bar)')
title('Presion vs Tiempo')
hold off

%% SISTEMA PRESION VS TIEMPO
%% Interpolacion de Datos

%Valores de identificacion gradual
P25=0.25*max(P);
P50=0.50*max(P);
P75=0.75*max(P);
PP=[P25 P50 P75];

%Barrido de vectores para interpolar la variable independiente
countP=1;
while countP~=4
for i = 1:size(P,1)
    if(P(i)<PP(countP)||P(i)==PP(countP))
        AP(countP)=i;
    end
end
BP(countP)=AP(countP)+1;
countP=countP+1;
end
i=0;

%Formacion vector de valores de tiempo interpolados
for i = 1:3
x1=T(AP(i));
x2=T(BP(i));
y1=P(AP(i));
y2=P(BP(i));
TTP(i) = (((x2-x1)/(y2-y1))*(PP(i)-y1))+x1; %Interpolacion de los tiempos para los delta Y.
end 
TTP

%% ULTIMO METODO DE SINTONIZACION
%% Modelo de polo doble más tiempo muerto (PDMTM) para Flujo
Kp2P=(max(P)-min(P))/(max(u)-min(u));
tao2P=0.5776*(TTP(3)-TTP(1));
tm2P=(1.5552*TTP(1))-(0.5552*TTP(3));

num2P=Kp2P*exp(-tm2P*s);
den2P=((tao2P*s)+1)^2;
Gp2P=num2P/den2P

%  Error de predicción cuadrático
[ye2P,te2P]=step(Gp2P,T);
i=0;
for i= 1:size(P,1)
S2(i,1)=(P(i)/20000-ye2P(i))^2;
end
S22=sum(S2)

%Comparación
ye2P=20000*ye2P;
figure
plot(T,P)
% plot(T,F/(max(u)-min(u)))
hold on
plot(T,ye2P)
legend('Real','Estimado2')
title(['PDMTM con Error de predición = ',num2str(S22)])

%% Para discretizar
[N,D]=pade(tm2P,1);
ret = tf(N,D);
Gf = (Kp2P*ret)/((tao2P*s)+1);
% zpk(Gf)

% Discretizacion de la planta
Ts = 11.09/20;
Gz = c2d(Gf,Ts,'tustin')
figure
step(Gp2P)
hold on
step(Gf)
hold on
step(Gz)
title(['Comparacion Referente a la Planta'])
legend('Continuo','Aproximacion','Discreto')
%% Método de Sung, O, Lee, Lee y Yi

z = tm2P/tao2P %% <0.9 y <0.4

%Guia encontrada
% Kc=((-0.04+(0.333+(0.949*((tm2P/tao2P)^-0.983))))*z)/Kp2P;
% Ti=((2.055+(0.072*(tm2P/tao2P)))*z)*tao2P;
Kc=((-0.544+(0.308*(z))+((1.408*((z)^-0.832))*z))/Kp2P)/2;
Ti=(1.768+(0.329*(z)))*z*tao2P
A=((1-exp(((-(tm2P/tao2P)^1.060)*z)/0.870))*(0.55+(1.683*((tm2P/tao2P)^-1.090))));
Td=tao2P/A

% Control PID ideal en continuo
Gc = Kc*(1+(1/(Ti*s))+((Td*s)/(1+(tao2P*s))));
cl_sys = feedback(Gc*Gf,1);

% Control PID Ideal en Discreto
Gc2z = c2d(Gc,Ts,'tustin')
figure
step(feedback(Gc*Gp2P,1))
hold on
step(cl_sys)
hold on
step(feedback(Gc2z*Gz,1))

title(['Control PID ideal de Sung, O, Lee, Lee y Yi'])
legend('Continuo','Aproximado','Discreto')

% Ecuacion en diferencias PID 
U2=ones(1,140);
Y(1)=0;
Y(2)=0;
ref=1;
for k = 3:140
    Y(k) = -(2.438e-05)*U2(k)+(8.929e-06)*U2(k-1)+(3.331e-05)*U2(k-2)+1.405*Y(k-1)-0.4927*Y(k-2);
    E(k) = ref-Y(k);
    U2(k) = 8062*E(k)-(1.43e04)*E(k-1)+6326*E(k-2)+1.673*U2(k-1)-0.6731*U2(k-2);
end
Yult=Y;
figure
plot(Yult,'g')
hold on
step(feedback(Gc*Gp2P,1))
hold on
step(feedback(Gc2z*Gz,1))
title (['Comparación PID Ideal de Sung O, Lee (Ecu.Diferencias)'])
legend('Diferencias','Continuo','Discreto')

