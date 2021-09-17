clc
clear all
close all

s = tf('s');

num = (1.25*exp(-0.25*s));
den = (16*s+1)*(4*s+1)*(2*s+1)*(s+1);
Gs=num/den

[y,t]=step(Gs)

y25=0.25*1.25 %0.25*max(y)
y50=0.50*1.25
y75=0.75*1.25
yy=[y25 y50 y75]'

% Interpolacion de datos:
count=1;
while count~=4
for i = 1:size(y,1)
    if(y(i)<yy(count)|| y(i)==yy(count))
        A(count)=i;
    end
end
B(count)=A(count)+1;
count=count+1;
end
i=0;

for i = 1:3
x1=t(A(i))
x2=t(B(i))
y1=y(A(i))
y2=y(B(i))
xx(i) = (((x2-x1)/(y2-y1))*(yy(i)-y1))+x1
end 

%% Modelo de primer orden más tiempo muerto (POMTM)
Kp=1.25/1;
tao=0.9102*(xx(3)-xx(1));
tm=(1.2620*xx(1))-(0.2620*xx(3));

num1=Kp*exp(-tm*s);
den1=(tao*s)+1;
Gp1=num1/den1

%  Error de predicción cuadrático
tnew=0:0.18:180;
[ye,te]=step(Gp1,tnew);
[y,t]=step(Gs,tnew)
for i= 1:size(ye,1)
S(i,1)=(ye(i)-y(i))^2;
end
S21=sum(S)

%  Comparación
figure
step(Gs)
hold on
step (Gp1)
legend('Real','Estimado1')
title(['POMTM con Error de predición = ',num2str(S21)])
%% Modelo de polo doble más tiempo muerto (PDMTM)
Kp2=1.25/1;
tao2=0.5776*(xx(3)-xx(1));
tm2=(1.5552*xx(1))-(0.5552*xx(3));

num2=Kp2*exp(-tm2*s);
den2=((tao2*s)+1)^2;
Gp2=num2/den2

%  Error de predicción cuadrático
[ye2,te2]=step(Gp2,tnew);
i=0;
for i= 1:size(y,1)
S2(i,1)=(y(i)-ye2(i))^2;
end
S22=sum(S2)

%Comparación
figure
step(Gs)
hold on
step (Gp2)
legend('Real','Estimado2')
title(['PDMTM con Error de predición = ',num2str(S22)])

%%  Modelo de segundo orden más tiempo muerto (SOMTM)

%%% Metodo simplificado (SOMTMs)
Kp3 = 1.25/1;
tm3=tm2;
a=(xx(2)-tm2-(1.4362*tao2))/((1.9844*tao2)-xx(2)+tm2)
tao3=(2*tao2)/(1+a)
T1=tao3;
T2=a*tao3;

num3=Kp3*exp(-tm3*s);
den3=(tao3*s+1)*(a*tao3*s+1);
Gp3=num3/den3

%  Error de predicción cuadrático
[ye3,te3]=step(Gp3,tnew);
i=0;
for i= 1:size(y,1)
S3(i,1)=(y(i)-ye3(i))^2;
end
S23=sum(S3)

%Comparación
figure
step(Gs)
hold on
step (Gp3)
legend('Real','Estimado3')
title(['SOMTMs con Error de predición = ',num2str(S23)])

%% Método general (SOMTMg)
Kp4=1.25/1
a4=(-0.6240*xx(1)+0.9866*xx(2)-0.3626*xx(3))/(0.3533*xx(1)-0.7036*xx(2)+0.3503*xx(3))
tao4=(xx(3)-xx(1))/(0.9866+(0.7036*a4))
T14=tao4;
T24=a4*tao4;
tm4=xx(3)-((1.3421+1.3455*a4)*tao4)

num4=Kp4*exp(-tm4*s);
den4=(tao4*s+1)*(a4*tao4*s+1);
Gp4=num4/den4

%  Error de predicción cuadrático
[ye4,te4]=step(Gp4,tnew);
i=0;
for i= 1:size(y,1)
S4(i,1)=(y(i)-ye4(i))^2;
end
S24=sum(S4)

%Comparación
figure
step(Gs)
hold on
step (Gp4)
legend('Real','Estimado3')
title(['SOMTMg con Error de predición = ',num2str(S24)])
