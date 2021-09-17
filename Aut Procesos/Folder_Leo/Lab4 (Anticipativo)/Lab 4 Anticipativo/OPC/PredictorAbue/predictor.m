clear all
close all
clc
load Datos

t = tout(:,1)%tiempo
f = FUJO(:,1)%flujo
p = PRESION(:,1)%presion
u = BOMBA(:,1)%entrada
n = NIVEL(:,1)%entrada

%Al inicio
t(1:200)=[];
f(1:200)=[];
p(1:200)=[];
u(1:200)=[];
%Al final
t(264:464)=[];
f(264:464)=[];
p(264:464)=[];
u(264:464)=[];

% figure 
% plot(t,f)
% hold on
  figure
  plot(t,p)
  hold on
  % ------(DELTA DE U Y DELTA DE Y)(DELTA DE Y*0.25,0.5,0.75)----FLUJO---///
 DuF=(u(220)-6000);% tomamos el dato de la fila 220 columna 1 que es la entrada escalon. le restamos 6000 para el delta de u 
 DyF=f(220);% tomamos el dato de la fila 220 que seria el flujo.

Dy25F=DyF*0.25;
Dy50F=DyF*0.5;
Dy75F=DyF*0.75;
 
%------INTERPOLACION LINEAL---FLUJO---
 t25F=((Dy25F-f(99))*(t(100)-t(99))/(f(100)-f(99)))+t(99);
 t50F=((Dy50F-f(110))*(t(111)-t(110))/(f(111)-f(110)))+t(110);
 t75F=((Dy75F-f(127))*(t(128)-t(127))/(f(128)-f(127)))+t(127);
 
 %------Modelo de primer orden mas tiempo muerto--FLUJO--
%Gp1(s)= (kp*exp(-tm*s))/T*s+1
 s=tf('s');
 kpF=DyF/DuF;
 taoF=0.9102*(t75F-t25F);
 tmF=(1.2620*t25F)-(0.2620*t75F);
 POMTM= (kpF*exp(-tmF*s))/(taoF*s+1);
 
 % ------(DELTA DE U Y DELTA DE Y)(DELTA DE Y*0.25,0.5,0.75)----PRESIÓN---///
 DuP=(u(220)-6000);% tomamos el dato de la fila 220 columna 1 que es la entrada escalon. le restamos 6000 para el delta de u 
 DyP=p(220);% tomamos el dato de la fila 220 que seria el flujo.

Dy25P=DyP*0.25;
Dy50P=DyP*0.5;
Dy75P=DyP*0.75;
 
%------INTERPOLACION LINEAL---PRESIÓN---
 t25P=((Dy25P-p(97))*(t(98)-t(97))/(p(98)-p(97)))+t(97);
 t50P=((Dy50P-p(100))*(t(101)-t(100))/(p(101)-p(100)))+t(100);
 t75P=((Dy75P-p(107))*(t(108)-t(107))/(p(108)-p(107)))+t(107);
 
 %------Modelo de primer orden mas tiempo muerto--PRESIÓN--////////////////
%Gp1(s)= (kp*exp(-tm*s))/T*s+1
 s=tf('s');
 kpP=DyP/DuP;
 taoP=0.9102*(t75P-t25P);
 tmP=(1.2620*t25P)-(0.2620*t75P);
 POMTMp= (kpP*exp(-tmP*s))/(taoP*s+1);
 
 %------Grafica---de---comparación---Identificación---PRESIÓN--////////////
 [y,x]=step((POMTMp)*21000,t);
 plot(x,y,'r')
 title ('SISTEMA-(PRESIÓN)')
 legend('Real','Identificación')
 hold off
 
 %------Grafica---de---comparación---Identificación---FLUJO-----///////////
 figure 
 plot(t,f)
 hold on
 [y,x]=step((POMTM)*21000,t);
 plot(x,y,'r')
 title ('SISTEMA-(FLUJO)')
 legend('Real','Identificación')
 hold off
 
%------Diseño----LGR----Controlador--PD-------PID---------FLUJO--//////////

ts=100
Mp=10
sigma=4/ts; % para 2% error
wd=(-pi*sigma)/log(Mp);%calculamos wd 
[num,den]=tfdata(pade(POMTM,1),'v');
[z,p,k]=tf2zp(num,den);% sacamos los ceros polos y ganancias del sistema

rs= -(atand(abs(wd)/(abs(p(1))-sigma)))+(-180+(atand(abs(wd)/(sigma-abs(p(2))))))-(-180+(atand(abs(wd)/(sigma+z))));%hacemos el calculo de la constante
a=(abs(wd)/tand(-180-rs))+sigma;%hacemos el calculo del cero
%calculamos el Kd la constante derivativa 
Kd=(sqrt((abs(wd)^2)+(abs(p(1))-sigma)^2)*sqrt((abs(wd)^2)+(sigma-abs(p(2)))^2)/(k*sqrt((abs(wd)^2)+(a-sigma)^2)*sqrt((abs(wd)^2)+(sigma+z)^2)));

PD=Kd*(s+a)
PID=Kd*((s+a)*(s+0.037))/s
PI=Kd*((s+0.037)/s)
figure
step(feedback(pade(POMTM,1)*PI,1))
figure
step(feedback(pade(POMTM,1)*PD,1))
hold on 
step(feedback(pade(POMTM,1)*PID,1))
hold on
step(feedback(pade(POMTM,1),1))
hold off
title ('DISEÑO-LGR-(FLUJO)')
legend('PD','PID','PLANTA')


%------Diseño----LGR----Controlador--PD-------PID---------PRESIÓN--////////

ts=100
Mp=10
sigma=4/ts; % para 2% error
wd=(-pi*sigma)/log(Mp);%calculamos wd 
[num,den]=tfdata(pade(POMTMp,1),'v');
[zP,pP,kP]=tf2zp(num,den);% sacamos los ceros polos y ganancias del sistema

rs= -(atand(abs(wd)/(abs(pP(1))-sigma)))+(-180+(atand(abs(wd)/(sigma-abs(pP(2))))))-(-180+(atand(abs(wd)/(sigma+zP))));%hacemos el calculo de la constante
a=(abs(wd)/tand(-180-rs))+sigma;%hacemos el calculo del cero
%calculamos el Kd la constante derivativa 
Kd=(sqrt((abs(wd)^2)+(abs(pP(1))-sigma)^2)*sqrt((abs(wd)^2)+(sigma-abs(pP(2)))^2)/(kP*sqrt((abs(wd)^2)+(a-sigma)^2)*sqrt((abs(wd)^2)+(sigma+zP)^2)));

PD=Kd*(s+a)
PID=Kd*((s+a)*(s+0.037))/s
PI=Kd*((s+0.037)/s)
figure
step(feedback(pade(POMTMp,1)*PI,1))
figure
step(feedback(pade(POMTMp,1)*PD,1))
hold on 
step(feedback(pade(POMTMp,1)*PID,1))
hold on
step(feedback(pade(POMTMp,1),1))
hold off
title ('DISEÑO-LGR-(PRESIÓN)')
legend('PD','PID','PLANTA')

figure
step(feedback((pade(POMTMp,1)),1))
hold off