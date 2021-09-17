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

% for i=1:size(T,1)
%     T(i)=T(i)-10;
% end
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
S(i,1)=(F(i)/(max(u)-min(u))-yeF(i))^2;
end
S21=sum(S)

%  Comparación
yeF=(max(u)-min(u))*yeF;
figure
plot(T,F)
% plot(T,F/max(u)-min(u))
hold on
plot(T,yeF)
legend('Real','Estimado1')
title(['POMTM con Error de predición = ',num2str(S21)])
xlabel('Tiempo(s)')
ylabel('Flujo(l/h)')
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
hold off
title(['Sistema real con aproximacion de tiempo muerto'])
legend('real','estimado')

%% Discretización de la planta.

Ts = 11.8/10
Gz = c2d(Gf,Ts,'tustin')
figure
step(Gf)
hold on
step(Gz)
hold off
%% Control PID por asignacion de polos.

% Descomposicion de la planta Gp1F
[num,den]=tfdata(Gf,'v');
a1=num(2)
a2=num(3)
b1=den(1)
b2=den(2)
b3=den(3)
% Dinamica del sistema
H = tf(num,den); 
figure 
pzmap(H)
grid on

% Condiciones de diseño
os=10
ts=2
zita=-log(os/100)/sqrt(pi^2+log(os/100)^2)
sigma=4/ts
wn=sigma/zita
wd=wn*sqrt(1-(zita^2))
theta=acos(zita)*180/pi 

Pd=[-zita*wn+1i*wn*sqrt(1-zita^2), -zita*wn-1i*wn*sqrt(1-zita^2)];
P1=(-zita*wn)*20;
% P2=0.142;
Pd1=[Pd P1];
Pds=poly(Pd1)

alpha=1; %Agresividad del filro pasa bajos

% Parametros del controlador PID
Kc = ((a1*b2)+(a2*(b1-Pds(1)))-(a1*Pds(2)))/(a1^2)
Ti1 = (a2*Kc)/Pds(4);
Ti2 = (a1*Kc)/(b3+(a2*Kc)-Pds(3)) %
Td = (b1-Pds(1))/(Kc*a1)

% ((a2*Kc)/Ti2)*2.193

%Renombrar parametros del PID con el filtro derivativo
d2= (Kc*Td)
d1= (Kc)
d0= Kc/Ti2;

%Ecuacion de control PID
ECpid=Kc*(1+(1/(Ti1*s))+(Td*s))

%Simulacion en lazo cerrado
H=((ECpid*Gf)/(1+ECpid*Gf));

figure
step(H)
title('Sistema en lazo cerrado PID')
figure
pzmap(H)
title ('Polos y ceros en lazo cerrado PID')