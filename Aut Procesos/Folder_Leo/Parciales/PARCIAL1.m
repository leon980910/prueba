clc
clear all
close all
s = tf('s');

t = [0:0.5:9]'
t(20)=10
T = [26;26;26.01;27.35;27.80;28.20;28.60;28.90;29.30;29.65;29.85;30.15;30.30;30.40;30.50;30.60;30.65;30.70;30.70;30.70]

%% Interpolacion de Datos

%Valores de identificacion gradual
T25=(0.25*(max(T)-min(T)))+26;
T50=(0.50*(max(T)-min(T)))+26;
T75=(0.75*(max(T)-min(T)))+26;
TT=[T25 T50 T75];

%Barrido de vectores para interpolar la variable independiente
countT=1;
while countT~=4
for i = 1:size(T,1)
    if(T(i)<TT(countT)||T(i)==TT(countT))
        AT(countT)=i;
    end
end
BT(countT)=AT(countT)+1;
countT=countT+1;
end
i=0;

%Formacion vector de valores de tiempo interpolados
for i = 1:3
x1=t(AT(i));
x2=t(BT(i));
y1=T(AT(i));
y2=T(BT(i));
tt(i) = (((x2-x1)/(y2-y1))*(TT(i)-y1))+x1; %Interpolacion de los tiempos para los delta Y.
end 
tt

%% Punto A
%% Modelo de primer orden más tiempo muerto (POMTM) para Flujo
KpF=(max(T)-min(T))/(57.74*0.8);
taoF=0.9102*(tt(3)-tt(1));
tmF=(1.2620*tt(1))-(0.2620*tt(3));

num1F=KpF*exp(-tmF*s);
den1F=(taoF*s)+1;
Gp1F=num1F/den1F %Funcion de transferencia estimada

%  Error de predicción cuadrático modelo identificado
% tnewF=0:0.4:40;
[yeT,teT]=step(Gp1F);
%  Comparación
yeT=57.74*yeT+26;
figure
plot(teT,yeT)
title(['Comportamiento del sistema identificado'])

%% Punto B
% Factor de controlabilidad:
Fcont = tmF/taoF

% Control PI
Kp1 = (0.9*(taoF/(KpF*tmF)))/1.5;
Ti1 = (tmF/0.3);
Gc1 = Kp1*(1+(1/(Ti1*s)));
cl_sys1 = feedback(Gc1*Gp1F,1);
figure
step(cl_sys1)
hold on
step(Gp1F)
title(['Control PI de Ziegler y Nichols'])

% Control PID ideal
Kp2 = (1.2*(taoF/(KpF*tmF)))/1.5;
Ti2 =  2*tmF;
Td2 = 0.5*tmF;
Gc2 = Kp2*(1+(1/(Ti2*s))+((Td2*s)));
cl_sys = feedback(Gc2*Gp1F,1);
figure
step(cl_sys)
hold on
step(Gp1F)
title(['Control PID ideal de Ziegler y Nichols'])

%% Punto C
 [N,D]=pade(tmF,1);
ret = tf(N,D);
Gf = (KpF*ret)/((taoF*s)+1)
num = [-0.1017 0.3]
den = [2.627 8.746 2.948]

step(Gf)
[A,B,C,D]=tf2ss(num,den)

Polos=roots(den) %Polos del sistema.

%Condiciones de diseño
os=10;
ts=7;
zita=-log(os/100)/sqrt(pi^2+log(os/100)^2);
sigma=4/ts;
wn=sigma/zita;
wd=wn*sqrt(1-(zita^2));
theta=acos(zita)*180/pi; 

Pd=[-zita*wn+1i*wn*sqrt(1-zita^2), -zita*wn-1i*wn*sqrt(1-zita^2)];
% P1=(-zita*wn)*20;

Pd1=[Pd];
Pds=poly(Pd1)

Ac=A^2+Pds(2)*A+Pds(3)*[1 0 ; 0 1;];

Cm = [B (A*B)];
Rankc=rank(Cm)

K = [0 1]*Cm^-1*Ac

% Kx = acker(A,B,[Pd(1);Pd(2);P1]) Para comprobar
g = 1/(C*inv(B*K-A)*B)

Aclz=A-(B*K);

[num1,den1]=ss2tf(Aclz,B,C,D);

Td = tf(num1,den1)
roots(den1)

figure
step(g*Td)
hold on
step(cl_sys1)
title('Comparacion ')

