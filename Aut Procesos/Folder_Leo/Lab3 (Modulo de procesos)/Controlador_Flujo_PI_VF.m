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
[N,D]=pade(tmF,2);
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

%% Control PID por asignacion de polos.

% Planta sin tiempo muerto
num2F=KpF;
den2F=(taoF*s)+1;
Gp2F=num2F/den2F; %Funcion de transferencia sin tiempo muerto
Gp = tf([0.005718],[1 0.1869]); %Adecuacion de la tf
[num,den]=tfdata(Gp,'v');

% Error en estado estacionario sin compensacion
K = num(2)/den(2);
essp = (1/(1+K))*100

% Controlador PI por rlocus
Gc=(34*(s+0.5))/s
figure
step(feedback(Gc*Gp,1),100)
title('Controlador PI metodo 1')

% Controlador PI metodo 2
error = 0.05
K = (den(2)-(error*den(2)))/(error*num(2))
error=0.0001
Ki = 1/(K*error*(num(2)/den(2)))
% PI = (K*(s+0.8))/s
PI =  K + (Ki/s)
% PI = (K*(s+(Ki/K)))/s

figure
step(feedback(PI*Gp,1),T)
title('Controlador PI para Flujo')


[N1,D1]=tfdata(PI,'v')
sim('PI_sim')

%% CONTROLADOR SMITH
[N,D]=pade(tmF,2);
ret = tf(N,D);
Gr = (KpF*ret)/((taoF*s)+1);
CL = feedback(PI,Gp)
Ceq = feedback(CL,-Gr)

figure
step(feedback(Ceq*Gr,1),T)
title('Controlador SMITH')

%% Ecuaciones en diferencias}

% Discretizacion de la planta
% Ts = 11.7/10
Ts= 0.01;
Gz = c2d(Gr,Ts,'tustin')
Gcz = c2d(Ceq,Ts,'tustin')

figure
step(feedback(Gcz*Gz,1),100)
title('Controlador SMITH Discreto')

[NGz,DGz]=tfdata(Gz,'v')
[NGcz,DGcz]=tfdata(Gcz,'v')
%%
U=ones(1,100);
% U(1)=0; U(2)=0; U(3)=0; U(4)=0; U(5)=0; U(6)=0; U(7)=0; U(8)=0; U(9)=0; U(10)=0;
Y(1)=0; Y(2)=0; Y(3)=0; Y(4)=0; Y(5)=0;
Aux=[0];

U0=0; U1=0; U2=0; U3=0; U4=0; U5=0;
Y1=0; Y2=0; Y3=0;
E1=0; E2=0; E3=0; E4=0; E5=0;
ref=1;
for k = 6:100
    Y(k) = NGz(1)*U(k)+NGz(2)*U(k-1)+NGz(3)*U(k-2)+NGz(4)*U(k-3)-DGz(2)*Y(k-1)-DGz(3)*Y(k-2)-DGz(4)*Y(k-3);
    E(k) = ref-Y(k);
    U(k) = NGcz(1)*E(k)+NGcz(2)*E(k-1)+NGcz(3)*E(k-2)+NGcz(4)*E(k-3)+NGcz(5)*E(k-4)+NGcz(6)*E(k-5)-...
        DGcz(2)*U(k-1)-DGcz(3)*U(k-2)-DGcz(4)*U(k-3)-DGcz(5)*U(k-4)-DGcz(6)*U(k-5);
    
    Y0 = NGz(1)*U0+NGz(2)*U1+NGz(3)*U2+NGz(4)*U3-DGz(2)*Y1-DGz(3)*Y2-DGz(4)*Y3;
    E0 = ref-Y0;
    U0 = NGcz(1)*E0+NGcz(2)*E1+NGcz(3)*E2+NGcz(4)*E3+NGcz(5)*E4+NGcz(6)*E5-...
        DGcz(2)*U1-DGcz(3)*U2-DGcz(4)*U3-DGcz(5)*U4-DGcz(6)*U5;
    Aux(k)=Y0;
    
    U5=U4; U4=U3; U3=U2; U2=U1; U1=U0;
    Y3=Y2; Y2=Y1; Y1=Y0;
    E5=E4; E4=E3; E3=E2; E2=E1; E1=E0;   
end
% k=1:1:100
figure
plot(Y,'g')
hold on
plot(Aux,'r')
%% Controlador PI metodo 3

Aasa=5;
% cero=0.6;
% sys = ((K*num(2))*(s+cero))/(s^2+(den(2)*s))
% figure
% rlocus(sys)
% title('Lugar geometrico de la planta sin control') 
% 
% os=10;
% zita=-log(os/100)/sqrt(pi^2+log(os/100)^2)
% sgrid(zita,0)
% 
% Kp = 9.78
% Ki = cero*Kp
% PI = (Kp*(s+(Ki/Kp)))/s
% % PI =  K + (Ki/s)
% % PI = (K*(s+(Ki/K)))/s
% figure
% step(feedback(PI*Gp,1),100)
% title('Controlador PI metodo 3')
% 
% %% PID AUTOMATICO
% 
% Kp=34.08
% Ki=19.13
% PI=Kp+(Ki/s)
% % figure
% % step(feedback(PI*Gp2F,1))

sim('PI_sim')