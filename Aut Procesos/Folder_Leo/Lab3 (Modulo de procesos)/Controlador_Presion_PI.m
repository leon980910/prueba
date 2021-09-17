%% Limpieza de Consola y reinicio
clear all
close all
clc
%% Carga y asignacion de datos
A = load('DATOS1.mat')
s = tf('s');
T = A.tout; %Temperatura
F = A.FUJO; %Flujo
P = A.PRESION; %Presion
u = A.BOMBA; %Entrada

%% Graficas sin tratamiento de datos
% plot(T,F)
% figure
% plot(T,P)
% figure
% plot(T,u)

%% Tratamiento de Datos
%Al inicio

T(1:52)=[];
F(1:52)=[];
P(1:52)=[];
u(1:52)=[];
%Al final

T(166:218)=[];
F(166:218)=[];
P(166:218)=[];
u(166:218)=[];

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
F25=0.25*max(P);
F50=0.50*max(P);
F75=0.75*max(P);
FF=[F25 F50 F75];

%Barrido de vectores para interpolar la variable independiente
countF=1;
while countF~=4
for i = 1:size(P,1)
    if(P(i)<FF(countF)||P(i)==FF(countF))
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
y1=P(AF(i));
y2=P(BF(i));
TTF(i) = (((x2-x1)/(y2-y1))*(FF(i)-y1))+x1; %Interpolacion de los tiempos para los delta Y.
end 
TTF
%% Modelo de primer orden más tiempo muerto (POMTM) para Flujo
KpF=(max(P)-min(P))/(max(u)-min(u));
taoF=0.9102*(TTF(3)-TTF(1));
tmF=(1.2620*TTF(1))-(0.2620*TTF(3));

num1F=KpF*exp(-tmF*s);
den1F=(taoF*s)+1;
Gp1F=num1F/den1F %Funcion de transferencia estimada

%  Error de predicción cuadrático modelo identificado
% tnewF=0:0.4:40;
[yeF,teF]=step(Gp1F,T);
for i= 1:size(yeF,1)
S(i,1)=(P(i)/(max(u)-min(u))-yeF(i))^2;
end
S21=sum(S)

%  Comparación
yeF=(max(u)-min(u))*yeF;
figure
plot(T,P)
% plot(T,F/max(u)-min(u))
hold on
plot(T,yeF)
legend('Real','Estimado1')
title(['POMTM con Error de predición = ',num2str(S21)])
xlabel('Tiempo(s)')
ylabel('Flujo(l/h)')
% pidtool(Gp1F)

%% Controlador PI metodo 2

% Planta sin tiempo muerto
num2F=KpF;
den2F=(taoF*s)+1;
Gp2F=num2F/den2F; %Funcion de transferencia sin tiempo muerto
Gp = tf([9.217e-5],[1 0.71839]); %Adecuacion de la tf
[num,den]=tfdata(Gp,'v');

error = 0.03
K = (den(2)-(error*den(2)))/(error*num(2))
Ki=10*K;
% error=0.00001
% Ki = 1/(K*error*(num(2)/den(2)))
% PI = (K*(s+10))/s
% K=8125;
% Ki=1.752e4;
%PI =  K + (Ki/s)
%PI = (K*(s+(Ki/K)))/s
%PI=(1085*(s+0.7184))/s
PI=(1390*s+2281.33)/(s);
figure
step(feedback(PI*Gp,1),T)
title('Controlador PI para sistema de Presión')

[N1,D1]=tfdata(PI,'v')
sim('PI_sim')
%% CONTROLADOR SMITH
[N,D]=pade(tmF,5);
ret = tf(N,D);
Gr = (KpF*ret)/((taoF*s)+1);
CL = feedback(PI,Gp)
Ceq = feedback(CL,-Gr)

figure
step(feedback(Ceq*Gr,1),100)
title('Controlador SMITH')
% figure;step(Gp1F);hold on;step(Gr)

%% Ecuaciones en diferencias}

% Discretizacion de la planta
% Ts = 3.06/10
Ts=0.7;
% Ts= 10.4/25;
% Ts=0.01;
Gz = c2d(Gr,Ts,'tustin')
Gcz = c2d(Ceq,Ts,'tustin')

figure
step(feedback(Gcz*Gz,1),100)
title('Controlador SMITH Discreto')

[NGz,DGz]=tfdata(Gz,'v')
[NGcz,DGcz]=tfdata(Gcz,'v')

U=ones(1,100);
% U(1)=0; U(2)=0; U(3)=0; U(4)=0; U(5)=0; U(6)=0; U(7)=0; U(8)=0; U(9)=0; U(10)=0;
Y(1)=0; Y(2)=0; Y(3)=0; Y(4)=0; Y(5)=0; Y(6)=0; Y(7)=0; Y(8)=0; Y(9)=0;
Aux=[0];

U0=0; U1=0; U2=0; U3=0; U4=0; U5=0; U6=0; U7=0; U8=0; 
Y1=0; Y2=0; Y3=0; Y4=0; Y5=0; Y6=0;
E1=0; E2=0; E3=0; E4=0; E5=0; E6=0; E7=0; E8=0;
ref=2.5;
for k = 9:100
    Y(k) = NGz(1)*U(k)+NGz(2)*U(k-1)+NGz(3)*U(k-2)+NGz(4)*U(k-3)+NGz(5)*U(k-4)+NGz(6)*U(k-5)+NGz(7)*U(k-6)-...
        DGz(2)*Y(k-1)-DGz(3)*Y(k-2)-DGz(4)*Y(k-3)-DGz(5)*Y(k-4)-DGz(6)*Y(k-5)-DGz(7)*Y(k-6);
    E(k) = ref-Y(k);
    U(k) = NGcz(1)*E(k)+NGcz(2)*E(k-1)+NGcz(3)*E(k-2)+NGcz(4)*E(k-3)+NGcz(5)*E(k-4)+NGcz(6)*E(k-5)+NGcz(7)*E(k-6)+NGcz(8)*E(k-7)+NGcz(9)*E(k-8)-...
        DGcz(2)*U(k-1)-DGcz(3)*U(k-2)-DGcz(4)*U(k-3)-DGcz(5)*U(k-4)-DGcz(6)*U(k-5)-DGcz(7)*U(k-6)-DGcz(8)*U(k-7)-DGcz(9)*U(k-8);
    
    Y0 = NGz(1)*U0+NGz(2)*U1+NGz(3)*U2+NGz(4)*U3+NGz(5)*U4+NGz(6)*U5+NGz(7)*U6-...
        DGz(2)*Y1-DGz(3)*Y2-DGz(4)*Y3-DGz(5)*Y4-DGz(6)*Y5-DGz(7)*Y6;
    E0 = ref-Y0;
    U0 = NGcz(1)*E0+NGcz(2)*E1+NGcz(3)*E2+NGcz(4)*E3+NGcz(5)*E4+NGcz(6)*E5+NGcz(7)*E6+NGcz(8)*E7+NGcz(9)*E8-...
        DGcz(2)*U1-DGcz(3)*U2-DGcz(4)*U3-DGcz(5)*U4-DGcz(6)*U5-DGcz(7)*U6-DGcz(8)*U7-DGcz(9)*U8;
    Aux(k)=Y0;
%     
    U8=U7; U7=U6; U6=U5; U5=U4; U4=U3; U3=U2; U2=U1; U1=U0;
    Y6=Y5; Y5=Y4; Y4=Y3; Y3=Y2; Y2=Y1; Y1=Y0;
    E8=E7; E7=E6; E6=E5; E5=E4; E4=E3; E3=E2; E2=E1; E1=E0;   
end

figure
% plot(Y,'g')
% hold on
plot(Aux,'r')
%% Controlador PI metodo 3

figure
plot(U)
