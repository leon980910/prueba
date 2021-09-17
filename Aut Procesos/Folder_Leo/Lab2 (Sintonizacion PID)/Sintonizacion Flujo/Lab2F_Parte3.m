clc 
clear all
close all

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

% pidtool(Gp1F)
%% Para discretizar
[N,D]=pade(tmF,1);
ret = tf(N,D);
Gf = (KpF*ret)/((taoF*s)+1);
% zpk(Gf)

% Discretizacion de la planta
Ts = 12.9/20;
Gz = c2d(Gf,Ts,'tustin')

%% Método de Kaya y Sheib

%% Integral  del  error  absoluto  (IAE)
a1= 0.98089; b1=-0.76167; c1=0.91032; d1=-1.05211; e1=0.59974; f1=0.89819;

% Control PID ideal (IAE)
KcAL = ((a1*((tmF/taoF)^b1))/KpF);
TiAL = (1/(c1+(d1*(tmF/taoF))))*taoF;
TdAL = (e1*((tmF/taoF)^f1))*taoF;

GcAL = KcAL*(1+(1/(TiAL*s))+((TdAL*s)/(1+(taoF*f1*s))));

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
    Y1(k) = -0.000961*U1(k)+0.000488*U1(k-1)+0.001449*U1(k-2)+1.554*Y1(k-1)-0.5911*Y1(k-2); %Y/U=Gz
    E1(k) = ref-Y1(k); %Y(k)=Ref-Y(k)*Bettaz
    U1(k) = 79.03*E1(k)-149.1*E1(k-1)+70.3*E1(k-2)+1.88*U1(k-1)-0.8797*U1(k-2); %U/E = Alpha 
end

figure
plot(Y1,'g')
hold on
step(feedback(GcAL*Gp1F,1))
hold on
step(feedback(GcALz*Gz,1))
title (['Comparación PID Kaya(IAE)'])
legend('Diferencias','Continuo','Discreto')

% Control PID Industrial (IAE)
a1= 0.91000; b1=-0.79380; c1=1.01495; d1=-1.00403; e1=0.54140; f1=0.7848;

KcAL = ((a1*((tmF/taoF)^b1))/KpF);
TiAL = (1/(c1+(d1*(tmF/taoF))))*taoF;
TdAL = (e1*((tmF/taoF)^f1))*taoF;

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
    Y2(k) = -0.000961*U2(k)+0.000488*U2(k-1)+0.001449*U2(k-2)+1.554*Y2(k-1)-0.5911*Y2(k-2); %Y/U=Gz
    Ec2(k) = ref-((1.274*Y2(k)-0.93*Y2(k-1)+0.6558*Y2(k-1))); %Y(k)=Ref-Y(k)*Bettaz
    U2(k) = 54.88*Ec2(k)-52.14*Ec2(k-1)+U2(k-1); %U/E = Alpha 
end

figure
plot(Y2,'g')
hold on
step(feedback(Alpha*Gp1F,Betta))
hold on
step(sys_1z)
title (['Comparacion PID Industrial Kaya (IAE)'])
legend('Diferencias','Continuo','Discreto')

%% Integral del error absoluto por el tiempo (ITAE) 
a1= 0.77902; b1=-1.06401; c1=1.14311; d1=-0.70949; e1=0.57137; f1=1.03826;

% Control PID ideal (ITAE)
KcTL = ((a1*((tmF/taoF)^b1))/KpF);
TiTL = (1/(c1+(d1*(tmF/taoF))))*taoF;
TdTL = (e1*((tmF/taoF)^f1))*taoF;

GcTL = KcTL*(1+(1/(TiTL*s))+((TdTL*s)/(1+(taoF*f1*s))));

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
    Y1T(k) = -0.000961*U1T(k)+0.000488*U1T(k-1)+0.001449*U1T(k-2)+1.554*Y1T(k-1)-0.5911*Y1T(k-2); %Y/U=Gz
    E1T(k) = ref-Y1T(k); %Y(k)=Ref-Y(k)*Bettaz
    U1T(k) = 71.11*E1T(k)-131.8*E1T(k-1)+61.2*E1T(k-2)+1.895*U1T(k-1)-0.8951*U1T(k-2); %U/E = Alpha 
end

figure
plot(Y1T,'g')
hold on
step(feedback(GcTL*Gp1F,1))
hold on
step(feedback(GcTLz*Gz,1))
title (['Comparación PID Kaya (ITAE)'])
legend('Diferencias','Continuo','Discreto')

% Control PID Industrial (ITAE)
a1= 0.70580; b1=-0.88720; c1=1.03326; d1=-0.99138; e1=0.60006; f1=0.97100;

KcTL = ((a1*((tmF/taoF)^b1))/KpF);
TiTL = (1/(c1+(d1*(tmF/taoF))))*taoF;
TdTL = (e1*((tmF/taoF)^f1))*taoF;

AlpT=KcTL*(1+1/(TiTL*s));
BettaT=((1+(TdTL*s))/(1+(f1*TdTL*s)));
sys_1T=feedback(AlpT*Gp1F,BettaT);

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
    Y2T(k) = -0.000961*U2T(k)+0.000488*U2T(k-1)+0.001449*U2T(k-2)+1.554*Y2T(k-1)-0.5911*Y2T(k-2); %Y/U=Gz
    E2T(k) = ref-((1.03*Y2T(k)-0.7403*Y2T(k-1)+0.7104*Y2T(k-1))); %Y(k)=Ref-Y(k)*Bettaz
    U2T(k) = 44.94*E2T(k)-42.57*E2T(k-1)+U2T(k-1); %U/E = Alpha 
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
KcSL = ((a1*((tmF/taoF)^b1))/KpF);
TiSL = (1/(c1+(d1*(tmF/taoF))))*taoF;
TdSL = (e1*((tmF/taoF)^f1))*taoF;

GcSL = KcSL*(1+(1/(TiSL*s))+((TdSL*s)/(1+(taoF*f1*s))));

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
    Y1S(k) = -0.000961*U1S(k)+0.000488*U1S(k-1)+0.001449*U1S(k-2)+1.554*Y1S(k-1)-0.5911*Y1S(k-2); %Y/U=Gz
    E1S(k) = ref-Y1S(k); %Y(k)=Ref-Y(k)*Bettaz
    U1S(k) = 95.64*E1S(k)-180.5*E1S(k-1)+85.13*E1S(k-2)+1.877*U1S(k-1)-0.8771*U1S(k-2); %U/E = Alpha 
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

KcSL = ((a1*((tmF/taoF)^b1))/KpF);
TiSL = (1/(c1+(d1*(tmF/taoF))))*taoF;
TdSL = (e1*((tmF/taoF)^f1))*taoF;

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
    Y2S(k) = -0.000961*U2S(k)+0.000488*U2S(k-1)+0.001449*U2S(k-2)+1.554*Y2S(k-1)-0.5911*Y2S(k-2); %Y/U=Gz
    E2S(k) = ref-((1.098*Y2S(k)-0.7856*Y2S(k-1)+0.688*Y2S(k-1))); %Y(k)=Ref-Y(k)*Bettaz
    U2S(k) = 71.31*E2S(k)-67.83*E2S(k-1)+U2S(k-1); %U/E = Alpha 
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
