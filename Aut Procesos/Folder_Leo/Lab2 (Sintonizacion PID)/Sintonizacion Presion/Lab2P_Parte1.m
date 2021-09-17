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

%% Para discretizar

[N,D]=pade(tmP,1);
ret = tf(N,D);
Gf = (KpP*ret)/((taoP*s)+1);
% zpk(Gf)

% Discretizacion
Ts = 4.91/14;
Gz = c2d(Gf,Ts,'tustin');
figure
step(Gp1P)
hold on
step(Gf)
hold on
step(Gz,'g')
title(['Sistema real con aproximacion de tiempo muerto'])
legend('REAL','APROXIMADO','DISCRETIZACION')

%% Método de Ziegler y Nichols
% Factor de controlabilidad:
Fcont = tmP/taoP

% Control PI
Kp1 = 0.9*(taoP/(KpP*tmP));
Ti1 = (tmP/0.3);
Gc1 = Kp1*(1+(1/(Ti1*s)));
cl_sys1 = feedback(Gc1*Gp1P,1);

% Discretizacion PI ideal
Gcz = c2d(Gc1,Ts,'tustin')

figure
step(cl_sys1)
hold on
step(feedback(Gcz*Gz,1))
title(['Control PI de Ziegler y Nichols'])
legend('Continuo','Discreto')

% Ecuacion en diferencias PI 
U=ones(1,200);
y(1)=0;
y(2)=0;
ref=1;
for k = 3:200
    y(k) = ((-1.259e-05)*U(k))+((2.234e-06)*U(k-1))+((1.482e-05)*U(k-2))+((1.704)*y(k-1))-(0.7258*y(k-2));
    e(k) = ref-y(k);
    U(k) = 2323*e(k)-2267*e(k-1)+U(k-1);
end

figure
plot(y,'r')
hold on
step(feedback(Gc1*Gp1P,1),200)
hold on
step(feedback(Gcz*Gz,1),200)
title (['Continuo vs discreto vs Ecu. Diferencias (PI Ziegler y Nichols)'])
legend('Diferencias','Continuo','Discreto')

% Control PID ideal
Kp2 = (1.2*(taoP/(KpP*tmP)));
Ti2 = 2*tmP;
Td2 = 0.5*tmP;
Gc2 = Kp2*(1+(1/(Ti2*s))+((Td2*s)/(1+(taoP*s))));
cl_sys = feedback(Gc2*Gp1P,1);

% Discretizacion PID Ideal
Gc2z = c2d(Gc2,Ts,'tustin')
% X11=c2d((Kp2*(1+(1/(Ti2*s)))),Ts,'tustin')
% Z11=c2d(((Td2*s)/(1+(taoP*s))),Ts,'zoh')
% Gc2z = X11+Z11
figure
step(cl_sys)
hold on
step(feedback(Gc2z*Gz,1))
hold on
step(feedback(Gc2*Gf,1))
title(['Control PID ideal de Ziegler y Nichols'])
legend('Continuo','Discreto','Aproximacion')

% Ecuacion en diferencias PID
U2=ones(1,120);
Y(1)=0;
Y(2)=0;
for k = 3:120
    Y(k) = -(1.259e-05)*U2(k)+(2.234e-06)*U2(k-1)+(1.482e-05)*U2(k-2)+1.704*Y(k-1)-0.7258*Y(k-2);
    E(k) = ref-Y(k);
    U2(k) = 5852*E(k)-(1.112e04)*E(k-1)+5291*E(k-2)+1.855*U2(k-1)-0.8546*U2(k-2);
end
Yzieg=Y;
figure
plot(Y,'g')
hold on
step(feedback(Gc2*Gp1P,1))
hold on
step(feedback(Gc2z*Gz,1))
title (['Contro PID Ideal de Ziegler y Nichols)'])
legend('Diferencias','Continuo','Discreto')

%% Método de Cohen y Coon

% Control PI Cohen y Coon
Kcc1 = ((taoP/(KpP*tmP))*(0.9+(tmP/(12*taoP))));
Tic1 = (tmP*(30+((3*tmP)/taoP)))/(9+((20*tmP)/taoP));
Gcc1 = Kcc1*(1+(1/(Tic1*s)));
cl_sysc1 = feedback(Gcc1*Gf,1);

% Discretizacion PI Cohen y Coon
Gcc1z = c2d(Gcc1,Ts,'tustin');
cl_sysc1z = feedback(Gcc1z*Gz,1);
figure
step(cl_sysc1)
hold on
step(cl_sysc1z)
title(['Control Discreto PI de Cohen y Coon'])

% Ecuacion en Diferencias Cohen y Coon
Uc1=ones(1,200);
Yc1(1)=0;
Yc1(2)=0;
for k = 3:200
    Yc1(k) = -(1.259e-05)*Uc1(k)+(2.234e-06)*Uc1(k-1)+(1.482e-05)*Uc1(k-2)+1.704*Yc1(k-1)-0.7258*Yc1(k-2);
    Ec1(k) = ref-Yc1(k);
    Uc1(k) = 2850*Ec1(k)-2557*Ec1(k-1)+Uc1(k-1);
end

figure
plot(Yc1,'g')
hold on
step(feedback(Gcc1*Gp1P,1))
hold on
step(feedback(Gcc1z*Gz,1))
title (['Comparacion PI Cohen y Coon'])
legend('Diferencias','Continuo','Discreto')

% Control PID ideal Cohen y Coon
Kcc = ((taoP/(KpP*tmP))*((4/3)+(tmP/(4*taoP))));
Tic = tmP*((32+(6*(tmP/taoP)))/(13+(8*(tmP/taoP))));
Tdc = tmP*(4/(11+(2*(tmP/taoP))));

Gcc = Kcc*(1+(1/(Tic*s))+((Tdc*s)/(1+(taoP*s))));
cl_sysc = feedback(Gcc*Gp1P,1);

%Dicretizacion PID ideal Cohen y Coon
Gcc2z=c2d(Gcc,Ts,'tustin')
figure
step(feedback(Gcc2z*Gz,1))
hold on
step(feedback(Gcc*Gf,1))
title(['Control Discreto PID Cohen y Coon'])
legend ('Discreto','Continuo Aprox.')

% Ecuacion en diferencias PID Cohen y con
Uc2=ones(1,140);
Yc2(1)=0;
Yc2(2)=0;
for k = 3:140
    Yc2(k) = -(1.259e-05)*Uc2(k)+(2.234e-06)*Uc2(k-1)+(1.482e-05)*Uc2(k-2)+1.704*Yc2(k-1)-0.7258*Yc2(k-2);
    Ec2(k) = ref-Yc2(k);
    Uc2(k) = 6973*Ec2(k)-(1.301e4)*Ec2(k-1)+6072*Ec2(k-2)+1.855*Uc2(k-1)-0.8546*Uc2(k-2);
end

figure
plot(Yc2,'g')
hold on
step(feedback(Gcc*Gp1P,1))
hold on
step(feedback(Gcc2z*Gz,1))
title (['Comparacion PID Cohen y Coon'])
legend('Diferencias','Continuo','Discreto')

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

%% Método de Sung, O, Lee, Lee y Yi

z = tm2P/tao2P %% <0.9 y <0.4

%Guia encontrada
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
plot(Yzieg,'b')
hold on
plot(Yc2,'r')
hold on
plot(Yult,'g')

title(['Comparacion de controles PID'])
legend('Ziegler y Nichols','Cohen y Coon','Sung O, Lee')