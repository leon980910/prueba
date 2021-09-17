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
%% Modelo de primer orden más tiempo muerto (POMTM) para Flujo
KpF=(max(F)-min(F))/(max(u)-min(u));
taoF=0.9102*(TTF(3)-TTF(1));
tmF=(1.2620*TTF(1))-(0.2620*TTF(3));

num1F=KpF*exp(-tmF*s);
den1F=(taoF*s)+1;
Gp1F=num1F/den1F %Funcion de transferencia estimada

%  Error de predicción cuadrático modelo identificado
% tnewF=0:0.4:40;
[yeF,teF]=step(Gp1F,T);
for i= 1:size(yeF,1)
S(i,1)=(F(i)/20000-yeF(i))^2;
end
S21=sum(S)

%  Comparación
yeF=20000*yeF;
figure
plot(T,F)
% plot(T,F/max(u)-min(u))
hold on
plot(T,yeF)
legend('Real','Estimado1')
title(['POMTM con Error de predición = ',num2str(S21)])

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

%%  Modelo de segundo orden más tiempo muerto (SOMTM) para Flujo

%%% Metodo simplificado (SOMTMs)
Kp3F = (max(F)-min(F))/(max(u)-min(u));
tm3F=tm2F;
a=(TTF(2)-tm2F-(1.4362*tao2F))/((1.9844*tao2F)-TTF(2)+tm2F);
tao3F=(2*tao2F)/(1+a);
T1F=tao3F;
T2F=a*tao3F;

num3F=Kp3F*exp(-tm3F*s);
den3F=(tao3F*s+1)*(a*tao3F*s+1);
Gp3F=num3F/den3F

%  Error de predicción cuadrático
[ye3F,te3F]=step(Gp3F,T);
i=0;
for i= 1:size(F,1)
S3(i,1)=(F(i)/20000-ye3F(i))^2;
end
S23=sum(S3)

%Comparación
ye3F=20000*ye3F;
figure
plot(T,F)
% plot(T,F/(max(u)-min(u)))
hold on
plot(T,ye3F)
legend('Real','Estimado3')
title(['SOMTMs con Error de predición = ',num2str(S23)])

%% Método general (SOMTMg) para Flujo
Kp4F=(max(F)-min(F))/(max(u)-min(u));
a4=(-0.6240*TTF(1)+0.9866*TTF(2)-0.3626*TTF(3))/(0.3533*TTF(1)-0.7036*TTF(2)+0.3503*TTF(3));
tao4F=(TTF(3)-TTF(1))/(0.9866+(0.7036*a4));
T14F=tao4F;
T24F=a4*tao4F;
tm4F=TTF(3)-((1.3421+1.3455*a4)*tao4F);

num4F=Kp4F*exp(-tm4F*s);
den4F=(tao4F*s+1)*(a4*tao4F*s+1);
Gp4F=num4F/den4F

%  Error de predicción cuadrático
[ye4F,te4F]=step(Gp4F,T);
i=0;
for i= 1:size(F,1)
S4(i,1)=(F(i)/20000-ye4F(i))^2;
end
S24=sum(S4)

ye4F=20000*ye4F;
figure
plot(T,F)
% plot(T,F/(max(u)-min(u)))
hold on
plot(T,ye4F)
legend('Real','Estimado3')
title(['SOMTMg con Error de predición = ',num2str(S24)])