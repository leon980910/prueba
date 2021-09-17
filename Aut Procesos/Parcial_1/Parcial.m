%%Parcial de Procesos
% Datos
clc
clear all
close all
t = [1 2 3 4 5 6 7 12 22 32 42]';
t = t*60;
y = [0 0 2.2 3.9 5.3 6.3 7.1 9.2 9.9 9.9 10]';
s = tf('s');
vector = [t,y];
dy = max(y)-min(y);
du = 60;
u  = [du du du du du du du du du du du]';
dy25 = dy*0.25+min(y);
dy50 = dy*0.50+min(y);
dy75 = dy*0.75+min(y);
dt25 = CalcularPunto(vector,dy25);
dt50 = CalcularPunto(vector,dy50);
dt75 = CalcularPunto(vector,dy75);
% -----------------------------
% Modelo de primer orden más tiempo muerto (POMTM)
% interpolación
kp = dy / du;
tao =  0.9102*(dt75 - dt25);
tm =  1.2620*dt25 - 0.2620*dt75
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
Gp1 = (kp*exp(-tm*s))/(tao*s+1)
yi = lsim(Gp1,u,t)+26;
    plot(t,yi)
    hold on
    plot(t,y,'r')
    legend('Orden 1','Planta','Location','SouthEast')
% ------------------------------------------------------------------------
% Calculo del indice de aproximacion
% definir un tiempo para ambas señales
data_error = y-yi
error = sum(((y-yi).^2))
tam_dato = size(t);
error_promedio = error/tam_dato(1)
    title(['Identificación orden 1 con tiempo muerto Error: ',num2str(error_promedio),' '])
    legend('Orden 1 ','Planta','Location','SouthEast')
    xlabel('Tiempo[s]') 
    ylabel('Caudal')
infoy=stepinfo(feedback(Gp1,1));
TrId1=infoy.RiseTime;
TId1=TrId1/20;
GIpId1=pade(Gp1,1);
GIzId1=c2d(GIpId1,TId1,'tustin');
    figure
    step(feedback(GIpId1,1))
    hold on
    step(feedback(GIzId1,1))
    title(['Planta real en Continuo vs Discreto'])
    legend('Continuo','Discreto','Location','SouthEast')
    ylabel('Amplitude') 
    xlabel('Time [seconds]')
% ------------------------------------------------------------------------
%% Punto B
% Factor de controlabilidad:
Fcont = tm/tao;

%%PID Ziegler y Nichols
% Ecuaciones de Sintonizacion
Kc=(1.2*tao)/(kp*tm);
Ti=2*tm;
Td=0.5*tm;

PID=Kc*(1+(1/(Ti*s))+((Td*s)/(1+tao*s)));
infoy=stepinfo(feedback(Gp1,1));
Tr=infoy.RiseTime;
T=Tr/15;

GIp=pade(Gp1,1);
GIz=c2d(GIp,T,'tustin');
PIDz=c2d(PID,T,'tustin');
    figure
    step(feedback(PID*GIp,1))
    hold on
    step(feedback(PIDz*GIz,1))
    title(['Control PID de Ziegler y Nichols'])
    legend('Continuo','Discreto','Location','SouthEast')

% ------------------------------------------------------------------------
%% Punto C
Gp1F = tf(pade(Gp1,1));
num = [-0.1017 0.005]
den = [157.6/60 8.746 2.948]

figure
step(Gp1F);
[A,B,C,D]=tf2ss(num,den)

Polos = roots(den) %Polos del sistema

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
    step(feedback(PID*GIp,1))
    title('Comparacion ')
    legend('Ackerman','Nichols','Location','SouthEast')
figure
step(g*Td)