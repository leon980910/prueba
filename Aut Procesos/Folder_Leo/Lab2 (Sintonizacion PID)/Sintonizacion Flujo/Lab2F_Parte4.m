%% Limpieza de Consola y reinicio
clear all
close all
clc
%% Carga y asignacion de datos
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

%% SISTEMA FLUJO VS TIEMPO
%% Interpolacion de Datos

%Valores de identificacion gradual
F25=0.25*max(F);
F50=0.50*max(F);
F75=0.75*max(F);
FF=[F25 F50 F75];

%Barrido de vectores para interpolar la variable independiente
countF=1;
while countF~=4
for i = 1:size(F,1)
    if(F(i)<FF(countF)||F(i)==FF(countF))
        AF(countF)=i;
    end
end
BF(countF)=AF(countF)+1;
countF=countF+1;
end
i=0;

%Formacion vector de valores de tiempo interpolados
for i = 1:3
x1=T(AF(i));
x2=T(BF(i));
y1=F(AF(i));
y2=F(BF(i));
TTF(i) = (((x2-x1)/(y2-y1))*(FF(i)-y1))+x1; %Interpolacion de los tiempos para los delta Y.
end 
TTF

%% Modelo de polo doble más tiempo muerto (PDMTM) para Flujo
Kp2F=(max(F)-min(F))/(max(u)-min(u));
tao2F=0.5776*(TTF(3)-TTF(1));
tm2F=(1.5552*TTF(1))-(0.5552*TTF(3));

num2F=Kp2F*exp(-tm2F*s);
den2F=((tao2F*s)+1)^2;
Gp2F=num2F/den2F

%  Error de predicción cuadrático
[ye2F,te2F]=step(Gp2F,T);
i=0;
for i= 1:size(F,1)
S2(i,1)=(F(i)/20000-ye2F(i))^2;
end
S22=sum(S2)

%Comparación
ye2F=20000*ye2F;
figure
plot(T,F)
% plot(T,F/(max(u)-min(u)))
hold on
plot(T,ye2F)
legend('Real','Estimado2')
title(['PDMTM con Error de predición = ',num2str(S22)])

%% Para discretizar
[N,D]=pade(tm2F,1);
ret = tf(N,D);
Gf = (Kp2F*ret)/((tao2F*s)+1);
% zpk(Gf)
figure
step(Gp2F)
hold on
step(Gf)
title(['Planta Continuo vs Discreto'])
legend('Continuo','Aproximacion')

% Discretizacion de la planta
Ts = 8.01/20;
Gz = c2d(Gf,Ts,'tustin')

%% Método de Sung, O, Lee, Lee y Yi

z = tm2F/tao2F %% <0.9 y <0.4

%Guia encontrada
Kc=((-0.04+(0.333+(0.949*((tm2F/tao2F)^-0.983))))*z)/Kp2F;
Ti=((2.055+(0.072*(tm2F/tao2F)))*z)*tao2F;
A=((1-exp(((-(tm2F/tao2F)^1.060)*z)/0.870))*(0.55+(1.683*((tm2F/tao2F)^-1.090))));
Td=tao2F/A

% Control PID ideal en continuo
Gc = Kc*(1+(1/(Ti*s))+((Td*s)/(1+(tao2F*s))));
cl_sys = feedback(Gc*Gf,1);

% Control PID Ideal en Discreto
Gc2z = c2d(Gc,Ts,'tustin')
figure
step(cl_sys)
hold on
step(feedback(Gc2z*Gz,1))
title(['Control PID ideal de Sung, O, Lee, Lee y Yi'])
legend('Continuo','Discreto')

% Ecuacion en diferencias PI 
U2=ones(1,45);
Y(1)=0;
Y(2)=0;
ref=1;
for k = 3:45
    Y(k) = -0.0007803*U2(k)+0.0006391*U2(k-1)+0.001419*U2(k-2)+1.443*Y(k-1)-0.4912*Y(k-2);
    E(k) = ref-Y(k);
    U2(k) = 88.9*E(k)-167.8*E(k-1)+79.52*E(k-2)+1.893*U2(k-1)-0.8935*U2(k-2);
end

figure
plot(Y,'g')
hold on
step(cl_sys)
hold on
step(feedback(Gc2z*Gz,1))
title (['Comparación PID Ideal de Sung, O, Lee, Lee y Yi)'])
legend('Diferencias','Continuo','Discreto')
