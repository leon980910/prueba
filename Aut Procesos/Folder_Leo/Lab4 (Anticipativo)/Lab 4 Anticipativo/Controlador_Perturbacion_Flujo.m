%% Limpieza de Consola y reinicio
clear all
close all
clc
%% Carga y asignacion de datos
A = load('DATOS3.mat')
s = tf('s');
T = A.tout; %Temperatura
F = A.FUJO; %Flujo
P = A.PRESION; %Presion
u = A.VALVULA; %Entrada

%% Tratamiento de Datos
%Al inicio

T(1:287)=[];
F(1:287)=[];
P(1:287)=[];
u(1:287)=[];
%Al final

% T(173:312)=[];
% F(173:312)=[];
% P(173:312)=[];
% u(173:312)=[];

for i=1:size(T,1)
    T(i)=T(i)-57.2;
    F(i)=F(i)-452.51;
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
F25=0.25*(max(F)-min(F));
F50=0.50*(max(F)-min(F));
F75=0.75*(max(F)-min(F));
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
tmFd=10.4;
[Nd,Dd]=pade(tmFd,5);
retd = tf(Nd,Dd);
Grd = (KpF*retd)/((taoF*s)+1);
[Nfd,Dfd]=tfdata(Grd,'v')
%  Error de predicción cuadrático modelo identificado
% tnewF=0:0.4:40;
[yeF,teF]=step(Gp1F,T);
for i= 1:size(yeF,1)
S(i,1)=(F(i)/(max(u)-min(u))-yeF(i))^2;
end
S21=sum(S)

%  Comparación
yeF=((max(u)-min(u))*yeF);
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

%% Controlador

Gff=-((5.21*s+0.0125)/(8.124*s+0.01702))*exp(-7.06*s)
tmF=7.06
[Nff,Dff]=pade(tmF,5);
Grff = tf(Nff,Dff);




figure
step(feedback(Grff*Grd,1))

save('DatosDF')

%% Discretizacion Perturbacion y Controlador

Ts=0.7;

Gdd=c2d(Grd,Ts,'tustin');
Gffd=c2d(Grff,Ts,'tustin');


[Ndz,Ddz]=tfdata(Gdd,'v')
[Nffd,Dffd]=tfdata(Gffd,'v')


U=ones(1,100);
% U(1)=0; U(2)=0; U(3)=0; U(4)=0; U(5)=0; U(6)=0; U(7)=0; U(8)=0; U(9)=0; U(10)=0;
Y(1)=0; Y(2)=0; Y(3)=0; Y(4)=0; Y(5)=0;
Aux=[0];

U0=0; U1=0; U2=0; U3=0; U4=0; U5=0;
Y1=0; Y2=0; Y3=0;
E1=0; E2=0; E3=0; E4=0; E5=0;
ref=2.0;
for k = 6:100
    Y(k) = Ndz(1)*U(k)+Ndz(2)*U(k-1)+Ndz(3)*U(k-2)+Ndz(4)*U(k-3)-Ddz(2)*Y(k-1)-Ddz(3)*Y(k-2)-Ddz(4)*Y(k-3);
    E(k) = ref-Y(k);
    U(k) = Nffd(1)*E(k)+Nffd(2)*E(k-1)+Nffd(3)*E(k-2)+Nffd(4)*E(k-3)+Nffd(5)*E(k-4)+Nffd(6)*E(k-5)-...
        Dffd(2)*U(k-1)-Dffd(3)*U(k-2)-Dffd(4)*U(k-3)-Dffd(5)*U(k-4)-Dffd(6)*U(k-5);
    
    Y0 = Ndz(1)*U0+Ndz(2)*U1+Ndz(3)*U2+Ndz(4)*U3-Ddz(2)*Y1-Ddz(3)*Y2-Ddz(4)*Y3;
    E0 = ref-Y0;
    U0 = Nffd(1)*E0+Nffd(2)*E1+Nffd(3)*E2+Nffd(4)*E3+Nffd(5)*E4+Nffd(6)*E5-...
        Dffd(2)*U1-Dffd(3)*U2-Dffd(4)*U3-Dffd(5)*U4-Dffd(6)*U5;
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





figure
plot(U)

