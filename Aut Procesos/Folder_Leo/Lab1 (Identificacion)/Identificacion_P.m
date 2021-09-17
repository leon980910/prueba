clc
clear all
close all

load('Datos.mat')
s = tf('s');
T = Data(:,1); %Temperatura
F = Data(:,2); %Flujo
P = Data(:,3); %Presion
u = Data(:,4); %Entrada

%% Graficas sin tratamiento de datos
% plot(T,F)
% figure
% plot(T,P)
% figure
% plot(T,u)

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

%% Modelo de primer orden más tiempo muerto (POMTM) para Presion
KpP=(max(P)-min(P))/(max(u)-min(u));
taoP=0.9102*(TTP(3)-TTP(1));
tmP=(1.2620*TTP(1))-(0.2620*TTP(3));

num1P=KpP*exp(-tmP*s);
den1P=(taoP*s)+1;
Gp1P=num1P/den1P %Funcion de transferencia estimada

%  Error de predicción cuadrático modelo identificado
% tnewF=0:0.4:40;
[yeP,teP]=step(Gp1P,T);
% yeF=yeF;
for i= 1:size(yeP,1)
S(i,1)=((P(i)/20000)-yeP(i))^2;
end
S21=sum(S)

%  Comparación
yeP=20000*yeP;
figure
plot(T,P)
hold on
plot(T,yeP)
legend('Real','Estimado1')
title(['POMTM con Error de predición = ',num2str(S21)])

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

%%  Modelo de segundo orden más tiempo muerto (SOMTM) para Flujo

%%% Metodo simplificado (SOMTMs)
Kp3P = (max(P)-min(P))/(max(u)-min(u));
tm3P=tm2P;
a=(TTP(2)-tm2P-(1.4362*tao2P))/((1.9844*tao2P)-TTP(2)+tm2P);
tao3P=(2*tao2P)/(1+a);
T1P=tao3P;
T2P=a*tao3P;

num3P=Kp3P*exp(-tm3P*s);
den3P=(tao3P*s+1)*(a*tao3P*s+1);
Gp3P=num3P/den3P

%  Error de predicción cuadrático
[ye3P,te3P]=step(Gp3P,T);
i=0;
for i= 1:size(P,1)
S3(i,1)=(P(i)/20000-ye3P(i))^2;
end
S23=sum(S3)

%Comparación
ye3P=20000*ye3P;
figure
plot(T,P)
% plot(T,F/(max(u)-min(u)))
hold on
plot(T,ye3P)
legend('Real','Estimado3')
title(['SOMTMs con Error de predición = ',num2str(S23)])

%% Método general (SOMTMg) para Flujo
Kp4P=(max(P)-min(P))/(max(u)-min(u));
a4=(-0.6240*TTP(1)+0.9866*TTP(2)-0.3626*TTP(3))/(0.3533*TTP(1)-0.7036*TTP(2)+0.3503*TTP(3));
tao4P=(TTP(3)-TTP(1))/(0.9866+(0.7036*a4));
T14P=tao4P;
T24P=a4*tao4P;
tm4P=TTP(3)-((1.3421+1.3455*a4)*tao4P);

num4P=Kp4P*exp(-tm4P*s);
den4P=(tao4P*s+1)*(a4*tao4P*s+1);
Gp4P=num4P/den4P

%  Error de predicción cuadrático
[ye4P,te4P]=step(Gp4P,T);
i=0;
for i= 1:size(P,1)
S4(i,1)=(P(i)/20000-ye4P(i))^2;
end
S24=sum(S4)

ye4P=20000*ye4P;
figure
plot(T,P)
% plot(T,F/(max(u)-min(u)))
hold on
plot(T,ye4P)
legend('Real','Estimado3')
title(['SOMTMg con Error de predición = ',num2str(S24)])
