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
step(feedback(PI*Gp,1),10)
title('Controlador PI metodo 2')


[N1,D1]=tfdata(PI,'v')

%% CONTROLADOR SMITH
Gr = pade(Gp1F,1)
CL = feedback(PI,Gp)
Ceq = feedback(CL,-Gr)

figure
step(feedback(Ceq*Gr,1),100)
title('Controlador SMITH')

%% Ecuaciones en diferencias}

% Discretizacion de la planta
Ts = 35/25
Gz = c2d(Gr,Ts,'tustin')

[N2,D2]=ss2tf(Ceq.A,Ceq.B,Ceq.C,Ceq.D)
Ceq1=tf(N2,D2)
Gcz = c2d(Ceq1,Ts,'tustin')
% Gpz= c2d(Gp,Ts,'tustin')
% PIz = c2d(PI,Ts,'tustin')
% % Gp2Fz = c2d(Gp2F,Ts,'tustin')
% CLz = feedback(PIz,Gpz)
% 
% Ceqz = feedback(CLz,-Gz)

figure
step(feedback(Gcz*Gz,1))
title('Controlador SMITH Discreto')

U2=ones(1,45);
Y(1)=0;
Y(2)=0;
ref=1;
for k = 3:45
    Y(k) = -0.03058*U1(k)+0.00434*U2(k-1)+0.001419*U2(k-2)+1.443*Y(k-1)-0.4912*Y(k-2);
    E(k) = ref-Y(k);
    U2(k) = 88.9*E(k)-167.8*E(k-1)+79.52*E(k-2)+1.893*U2(k-1)-0.8935*U2(k-2);
end

% Controlador PI metodo 3
% cero=1;
% sys = ((K*num(2))*(s+cero))/(s^2+(den(2)*s))
% figure
% rlocus(sys)
% title('Lugar geometrico de la planta sin control') 
% 
% os=10;
% zita=-log(os/100)/sqrt(pi^2+log(os/100)^2)
% sgrid(zita,0)
% 
% Kp = 274
% Ki = cero*Kp
% PI = (Kp*(s+(Ki/Kp)))/s
% % PI =  K + (Ki/s)
% % PI = (K*(s+(Ki/K)))/s
% figure
% step(feedback(PI*Gp2F,1),10)
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