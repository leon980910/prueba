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

%% Método de Kaya y Sheib
%% Integral  del  error  absoluto  (IAE)
a1= 0.98089; b1=-0.76167; c1=0.91032; d1=-1.05211; e1=0.59974; f1=0.89819;

% Control PID ideal (IAE)
KcAL = ((a1*((tmP/taoP)^b1))/KpP);
TiAL = (taoP/c1)*((tmP/taoP)^-d1);
TdAL = (e1*((tmP/taoP)^f1))*taoP;

GcAL = KcAL*(1+(1/(TiAL*s))+((TdAL*s)/(1+(taoP*f1*s))));

GcALz = c2d(GcAL,Ts,'tustin');
figure
step(feedback(GcAL*Gf,1))
hold on
step(feedback(GcALz*Gz,1))
title(['Control PID ideal Kaya (IAE)'])
legend('Continuo','Discreto')

U1=ones(1,200);
Y1(1)=0;
Y1(2)=0;
ref=1;
for k = 3:200
    Y1(k) = -(1.259e-05)*U1(k)+(2.234e-06)*U1(k-1)+(1.482e-05)*U1(k-2)+1.704*Y1(k-1)-0.7258*Y1(k-2);   
    E1(k) = ref-Y1(k); %Y(k)=Ref-Y(k)*Bettaz
    U1(k) = 6260*E1(k)-(1.182e04)*E1(k-1)+5598*E1(k-2)+1.839*U1(k-1)-0.8394*U1(k-2); %U/E = Alpha 
end

figure
plot(Y1,'g')
hold on
step(feedback(GcAL*Gp1P,1))
hold on
step(feedback(GcALz*Gz,1))
title (['Comparación PID Kaya(IAE)'])
legend('Diferencias','Continuo','Discreto')

% Control PID Industrial (IAE)
a1= 0.91000; b1=-0.79380; c1=1.01495; d1=-1.00403; e1=0.54140; f1=0.7848;

KcAL = ((a1*((tmP/taoP)^b1))/KpP);
TiAL = (taoP/c1)*((tmP/taoP)^-d1);
TdAL = (e1*((tmP/taoP)^f1))*taoP;

Alpha=KcAL*(1+1/(TiAL*s));
Betta=((1+(TdAL*s))/(1+(f1*TdAL*s)));
sys_1=feedback(Alpha*Gf,Betta);

Alphaz=c2d(Alpha,Ts,'tustin');
Bettaz=c2d(Betta,Ts,'zoh');
sys_1z=feedback(Alphaz*Gz,Bettaz);
figure
step(sys_1)
hold on
step(sys_1z)
title(['Control PID Industrial de Kaya (IAE)'])
legend('Continuo','Discreto')

U2=ones(1,160);
Y2(1)=0;
Y2(2)=0;
ref=1;
for k = 3:160
    Y2(k) = -(1.259e-05)*U2(k)+(2.234e-06)*U2(k-1)+(1.482e-05)*U2(k-2)+1.704*Y2(k-1)-0.7258*Y2(k-2);
    Ec2(k) = ref-((1.274*Y2(k)-1.105*Y2(k-1)+0.831*Y2(k-1))); %Y(k)=Ref-Y(k)*Bettaz
    U2(k) = 3038*Ec2(k)-2818*Ec2(k-1)+U2(k-1); %U/E = Alpha 
end

figure
plot(Y2,'g')
hold on
step(feedback(Alpha*Gp1P,Betta))
hold on
step(sys_1z)
title (['Comparacion PID Industrial Kaya (IAE)'])
legend('Diferencias','Continuo','Discreto')

%% Integral del error absoluto por el tiempo (ITAE) 
a1= 0.77902; b1=-1.06401; c1=1.14311; d1=-0.70949; e1=0.57137; f1=1.03826;

% Control PID ideal (ITAE)
KcTL = ((a1*((tmP/taoP)^b1))/KpP);
TiTL = (taoP/c1)*((tmP/taoP)^-d1);
TdTL = (e1*((tmP/taoP)^f1))*taoP;

GcTL = KcTL*(1+(1/(TiTL*s))+((TdTL*s)/(1+(taoP*f1*s))));

GcTLz = c2d(GcTL,Ts,'tustin');
figure
step(feedback(GcTL*Gf,1))
hold on
step(feedback(GcTLz*Gz,1))
title(['Control PID ideal Kaya (ITAE)'])
legend('Continuo','Discreto')

U1T=ones(1,200);
Y1T(1)=0;
Y1T(2)=0;
for k = 3:200
    Y1T(k) = -(1.259e-05)*U1T(k)+(2.234e-06)*U1T(k-1)+(1.482e-05)*U1T(k-2)+1.704*Y1T(k-1)-0.7258*Y1T(k-2);
    E1T(k) = ref-Y1T(k); %Y(k)=Ref-Y(k)*Bettaz
    U1T(k) = 3935*E1T(k)-7373*E1T(k-1)+3468*E1T(k-2)+1.86*U1T(k-1)-0.8596*U1T(k-2); %U/E = Alpha 
end

figure
plot(Y1T,'g')
hold on
step(feedback(GcTL*Gp1P,1))
hold on
step(feedback(GcTLz*Gz,1))
title (['Comparación PID Kaya (ITAE)'])
legend('Diferencias','Continuo','Discreto')

% Control PID Industrial (ITAE)
a1= 0.70580; b1=-0.88720; c1=1.03326; d1=-0.99138; e1=0.60006; f1=0.97100;

KcTL = ((a1*((tmP/taoP)^b1))/KpP);
TiTL = (taoP/c1)*((tmP/taoP)^-d1);
TdTL = (e1*((tmP/taoP)^f1))*taoP;

AlpT=KcTL*(1+1/(TiTL*s));
BettaT=((1+(TdTL*s))/(1+(f1*TdTL*s)));
sys_1T=feedback(AlpT*Gp1P,BettaT);

AlpTz=c2d(AlpT,Ts,'tustin');
BettaTz=c2d(BettaT,Ts,'zoh');
sys_1Tz=feedback(AlpTz*Gz,BettaTz);
figure
step(sys_1T)
hold on
step(sys_1Tz)
title(['Control PID Industrial de Kaya (ITAE)'])
legend('Continuo','Discreto')

U2T=ones(1,160);
Y2T(1)=0;
Y2T(2)=0;
for k = 3:160
    Y2T(k) = -(1.259e-05)*U2T(k)+(2.234e-06)*U2T(k-1)+(1.482e-05)*U2T(k-2)+1.704*Y2T(k-1)-0.7258*Y2T(k-2);
    E2T(k) = ref-((1.03*Y2T(k)-0.8964*Y2T(k-1)+0.8665*Y2T(k-1))); %Y(k)=Ref-Y(k)*Bettaz
    U2T(k) = 2012*E2T(k)-1797*E2T(k-1)+U2T(k-1); %U/E = Alpha 
end

figure
plot(Y2T,'g')
hold on
step(sys_1T)
hold on
step(sys_1Tz)
title (['Comparacion PID Industrial Kaya (ITAE)'])
legend('Diferencias','Continuo','Discreto')

%% Integral del error cuadrático (ISE).

a1= 1.11907; b1=-0.89711; c1=0.79870; d1=-0.95480; e1=0.54766; f1=0.87798;

% Control PID ideal (ISE)
KcSL = ((a1*((tmP/taoP)^b1))/KpP);
TiSL = (taoP/c1)*((tmP/taoP)^-d1);
TdSL = (e1*((tmP/taoP)^f1))*taoP;

GcSL = KcSL*(1+(1/(TiSL*s))+((TdSL*s)));

GcSLz = c2d(GcSL,Ts,'tustin');
figure
step(feedback(GcSL*Gf,1))
hold on
step(feedback(GcSLz*Gz,1))
title(['Control PID ideal Kaya (ISE)'])
legend('Continuo','Discreto')

U1S=ones(1,200);
Y1S(1)=0;
Y1S(2)=0;
for k = 3:200
    Y1S(k) = -(1.259e-05)*U1S(k)+(2.234e-06)*U1S(k-1)+(1.482e-05)*U1S(k-2)+1.704*Y1S(k-1)-0.7258*Y1S(k-2);
    E1S(k) = ref-Y1S(k); %Y(k)=Ref-Y(k)*Bettaz
    U1S(k) = 41020*E1S(k)-(75540)*E1S(k-1)+34920*E1S(k-2)+U1S(k-2); %U/E = Alpha 
end

figure
plot(Y1S,'g')
hold on
step(feedback(GcSL*Gf,1))
hold on
step(feedback(GcSLz*Gz,1))
title (['Comparación PID ideal Kaya (ISE)'])
legend('Diferencias','Continuo','Discreto')

% Control PID Industrial (ISE)
a1= 1.11470; b1=-0.89920; c1=0.93240; d1=-0.87530; e1=0.56508; f1=0.91107;

KcSL = ((a1*((tmP/taoP)^b1))/KpP);
TiSL = (taoP/c1)*((tmP/taoP)^-d1);
TdSL = (e1*((tmP/taoP)^f1))*taoP;

AlpS=KcSL*(1+1/(TiSL*s));
BettaS=((1+(TdSL*s))/(1+(f1*TdSL*s)));
sys_1S=feedback(AlpS*Gf,BettaS);

AlpSz=c2d(AlpS,Ts,'tustin');
BettaSz=c2d(BettaS,Ts,'zoh');
sys_1Sz=feedback(AlpSz*Gz,BettaSz);
figure
step(sys_1S)
hold on
step(sys_1Sz)
title(['Control PID Industrial de Kaya (ISE)'])
legend('Continuo','Discreto')

U2S=ones(1,160);
Y2S(1)=0;
Y2S(2)=0;
for k = 3:160
    Y2S(k) = -(1.259e-05)*U2S(k)+(2.234e-06)*U2S(k-1)+(1.482e-05)*U2S(k-2)+1.704*Y2S(k-1)-0.7258*Y2S(k-2);
    E2S(k) = ref-((1.098*Y2S(k)-0.9355*Y2S(k-1)+0.8378*Y2S(k-1))); %Y(k)=Ref-Y(k)*Bettaz
    U2S(k) = 3154*E2S(k)-2950*E2S(k-1)+U2S(k-1); %U/E = Alpha 
end

figure
plot(Y2S,'g')
hold on
step(sys_1S)
hold on
step(sys_1Sz)
title (['Comparacion PID Industrial Kaya (ISE)'])
legend('Diferencias','Continuo','Discreto')

%% COMPARACIONES DE PID IDEAL
figure
plot(Y1,'b')
hold on
plot(Y1T,'r')
hold on
plot(Y1S,'g')
title (['Comparacion PID Ideal de Kaya (Ecu.Diferencias)'])
legend('IAE','ITAE','ISE')

%% COMPARACIONES DE PID INDUSTRIAL
figure
plot(Y2,'b')
hold on
plot(Y2T,'r')
hold on
plot(Y2S,'g')
title (['Comparacion PID Industriales de Kaya (Ecu.Diferencias)'])
legend('IAE','ITAE','ISE')