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

% pidtool(Gp1F)
%% Para discretizar
[N,D]=pade(tmF,1);
ret = tf(N,D);
Gf = (KpF*ret)/((taoF*s)+1);
% zpk(Gf)
figure
step(Gp1F)
hold on
step(Gf)
title(['Sistema real con aproximacion de tiempo muerto'])


%% Método de Ziegler y Nichols
% Factor de controlabilidad:
Fcont = tmF/taoF

% Control PI
Kp1 = 0.9*(taoF/(KpF*tmF));
Ti1 = (tmF/0.3);
Gc1 = Kp1*(1+(1/(Ti1*s)));
cl_sys1 = feedback(Gc1*Gp1F,1);
figure
step(cl_sys1)
title(['Control PI de Ziegler y Nichols'])

% Discretizacion PI ideal
Ts = 12.9/20;
Gz = c2d(Gf,Ts,'tustin')
Gcz = c2d(Gc1,Ts,'tustin')

figure
step(Gf)
hold on
step(Gz)
title(['Planta real en Continuo vs Discreto'])

figure
step(feedback(Gc1*Gf,1))
hold on
step(feedback(Gcz*Gz,1))
title(['Control PI continuo vs Discreto'])

% Ecuacion en diferencias PI 
U=ones(1,70);
y(1)=0;
y(2)=0;
ref=1;
for k = 3:70
    y(k) = -0.000961*U(k)+0.000488*U(k-1)+0.001449*U(k-2)+1.554*y(k-1)-0.5911*y(k-2);
    e(k) = ref-y(k);
    U(k) = 61.27*e(k)-57.66*e(k-1)+1*U(k-1);
end

figure
plot(y,'r')
hold on
step(feedback(Gc1*Gp1F,1))
hold on
step(feedback(Gcz*Gz,1))
title (['Continuo vs discreto vs Ecu. Diferencias (PI Ziegler y Nichols)'])
legend('Diferencias','Continuo','Discreto')

% Control PID ideal
Kp2 = (1.2*(taoF/(KpF*tmF)))/1.5;
Ti2 =  2*tmF;
Td2 = 0.5*tmF;
Gc2 = Kp2*(1+(1/(Ti2*s))+((Td2*s)/(1+(taoF*s))));
cl_sys = feedback(Gc2*Gp1F,1);
figure
step(cl_sys)
title(['Control PID ideal de Ziegler y Nichols'])

%Discretizacion PID Ideal
Gc2z = c2d(Gc2,Ts,'tustin')
figure
step(feedback(Gc2z*Gz,1))

% Ecuacion en diferencias PID
U2=ones(1,45);
Y(1)=0;
Y(2)=0;
for k = 3:45
    Y(k) = -0.000961*U2(k)+0.000488*U2(k-1)+0.001449*U2(k-2)+1.554*Y(k-1)-0.5911*Y(k-2);
    E(k) = ref-Y(k);
%     U2(k) = 104.6*E(k)-192.1*E(k-1)+88.38*E(k-2)+1.891*U2(k-1)-0.8913*U2(k-2);
    U2(k) = 69.73*E(k)-128.1*E(k-1)+58.92*E(k-2)+1.891*U2(k-1)-0.8913*U2(k-2);
%     U2(k) = 474.9*E(k)-775.1*E(k-1)+316.3*E(k-2)+1*U2(k-2);
end
Yziger=Y;
figure
plot(Y,'g')
hold on
step(feedback(Gc2*Gp1F,1))
hold on
step(feedback(Gc2z*Gz,1))
title (['Continuo vs discreto vs Ecu. Diferencias (PID Ziegler y Nichols)'])
legend('Diferencias','Continuo','Discreto')

%% Método de Cohen y Coon

% Control PI Cohen y Coon
Kcc1 = ((taoF/(KpF*tmF))*(0.9+(tmF/(12*taoF))));
Tic1 = (tmF*(30+((3*tmF)/taoF)))/(9+((20*tmF)/taoF));

Gcc1 = Kcc1*(1+(1/(Tic1*s)));
cl_sysc1 = feedback(Gcc1*Gp1F,1);
figure
step(cl_sysc1)
title(['Control PI de Cohen y Coon'])

% Discretizacion PI Cohen y Coon
Gcc1z = c2d(Gcc1,Ts,'tustin');
cl_sysc1z = feedback(Gcc1z*Gz,1);
figure
step(cl_sysc1)
hold on
step(cl_sysc1z)
title(['Control Discreto PI de Cohen y Coon'])

% Ecuacion en Diferencias Cohen y Coon
Uc1=ones(1,70);
Yc1(1)=0;
Yc1(2)=0;
for k = 3:70
    Yc1(k) = -0.000961*Uc1(k)+0.000488*Uc1(k-1)+0.001449*Uc1(k-2)+1.554*Yc1(k-1)-0.5911*Yc1(k-2);
    Ec1(k) = ref-Yc1(k);
    Uc1(k) = 66.66*Ec1(k)-58.53*Ec1(k-1)+Uc1(k-1);
end

figure
plot(Yc1,'g')
hold on
step(feedback(Gcc1*Gp1F,1))
hold on
step(feedback(Gcc1z*Gz,1))
title (['Comparacion PI Cohen y Coon'])
legend('Diferencias','Continuo','Discreto')

% Control PID ideal Cohen y Coon
Kcc = ((taoF/(KpF*tmF))*((4/3)+(tmF/(4*taoF))))/2;
Tic = tmF*((32+(6*(tmF/taoF)))/(13+(8*(tmF/taoF))));
Tdc = tmF*(4/(11+(2*(tmF/taoF))));

Gcc = Kcc*(1+(1/(Tic*s))+((Tdc*s)/(1+(taoF*s))));
cl_sysc = feedback(Gcc*Gp1F,1);
figure
step(cl_sysc)
title(['Control PID ideal de Cohen y Coon'])

%Dicretizacion PID ideal Cohen y Coon
Gcc2z=c2d(Gcc,Ts,'tustin')
figure
step(feedback(Gcc2z*Gz,1))
hold on
step(feedback(Gcc*Gf,1))
title(['Control Discreto PID Cohen y Coon'])

% Ecuacion en diferencias PID Cog
Uc2=ones(1,45);
Yc2(1)=0;
Yc2(2)=0;
for k = 3:45
    Yc2(k) = -0.000961*Uc2(k)+0.000488*Uc2(k-1)+0.001449*Uc2(k-2)+1.554*Yc2(k-1)-0.5911*Yc2(k-2);
    Ec2(k) = ref-Yc2(k);
    Uc2(k) = 59.81*Ec2(k)-109.2*Ec2(k-1)+49.89*Ec2(k-2)+1.891*Uc2(k-1)-0.8913*Uc2(k-2);
end

figure
plot(Yc2,'g')
hold on
step(feedback(Gcc*Gp1F,1))
hold on
step(feedback(Gcc2z*Gz,1))
title (['Comparacion PID Cohen y Coon'])
legend('Diferencias','Continuo','Discreto')

%% ULTIMO METODO DE SINTONIZACION
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
legend('Continuo','Aproximado')

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

% Ecuacion en diferencias PID 
U2=ones(1,45);
Y(1)=0;
Y(2)=0;
ref=1;
for k = 3:45
    Y(k) = -0.0007803*U2(k)+0.0006391*U2(k-1)+0.001419*U2(k-2)+1.443*Y(k-1)-0.4912*Y(k-2);
    E(k) = ref-Y(k);
    U2(k) = 88.9*E(k)-167.8*E(k-1)+79.52*E(k-2)+1.893*U2(k-1)-0.8935*U2(k-2);
end
Yult=Y;
figure
plot(Yult,'g')
hold on
step(cl_sys)
hold on
step(feedback(Gc2z*Gz,1))
title (['Comparación PID Ideal de Sung, O, Lee, Lee y Yi)'])
legend('Diferencias','Continuo','Discreto')

%% COMPARACIONES
%% Comparacion PI en Ecuaciones en diferencias

figure
plot(y,'r')
hold on
plot(Yc1,'g')

title(['Comparacion de controles PI'])
legend('Ziegler y Nichols','Cohen y Coon')

%% Comparacion PID ideal

figure
plot(Yziger,'b')
hold on
plot(Yc2,'r')
hold on
plot(Yult,'g')

title(['Comparacion de controles PID'])
legend('Ziegler y Nichols','Cohen y Coon','Sung O, Lee')