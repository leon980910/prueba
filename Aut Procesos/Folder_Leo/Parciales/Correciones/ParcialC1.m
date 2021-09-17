clc
clear all
close all
s = tf('s');

Gp = ((0.1017)*exp(-0.678*s))/((2.627*s)+1)

Vs1 = ((0.1017))/((2.627*s)+1)
%% Diseño de controlador PI

tss = 5;
Ov = 0.10;

Z = abs(log(Ov))/sqrt((pi^2)+((log(Ov))^2));
Wn = 3/(Z*tss);

Pr = -Z*Wn; %Parte real 
Pi =(Wn*sqrt(1-(Z^2))); %Parte imaginaria
%Pd = Pr+j*Pi   %Polo deseado

%% Comprobar si el polo deseado pertenece a LGR
a = Pr+i*Pi;
Ang = ((0.1017)/(2.627*a+1)) *(1/a);
[theta,rho]=cart2pol(real(Ang),imag(Ang));
theta = theta*(180/pi)

% theta no es multiplo de 180°, se agrega un cero de compensación.
if theta<180
theta1 = 180-theta  %Angulo del cero aportado al LGR
cero = Pi/tand(theta1); %Ubicacion del cero.
end
C1 = (s+cero)/s %Controlador sin ganancia

%Ganancia Kp del controlador PI
Kp = 1/(abs(((0.1017)/((2.627*a)+1))*((a+cero)/a)))

%% Controlador diseñado
Cs = Kp*(s+cero)/s;

figure
step(Vs1)
hold on
step(feedback(Cs*Vs1,1))
title(['Planta sin control vs con control'])
legend('Sin control','Con control')
ylabel('Viscosidad')


%% Mitigar perturbaciones = Controlador Smith

tmF = 0.678;
KpF=0.1017;
taoF=2.627;

[N1,D1]=pade(tmF,2);
ret = tf(N1,D1);
Gr = (KpF*ret)/((taoF*s)+1);
CL = feedback(Cs,Vs1)
Ceq = feedback(CL,-Gr)

figure
step(feedback((Ceq*Gr),1))
% hold on
% step(-(Gper*0.3),3)
title('Controlador SMITH')

%% Punto B
% Factor de controlabilidad:
Fcont = tmF/taoF

% Control PI
Kp1 = (0.9*(taoF/(KpF*tmF)))/1.5;
Ti1 = (tmF/0.3);
Gc1 = Kp1*(1+(1/(Ti1*s)));
cl_sys1 = feedback(Gc1*Gp,1);
figure
step(cl_sys1)
hold on
step(Gp)
title(['Control PI de Ziegler y Nichols'])

figure
step(cl_sys1)
hold on
step(feedback((Ceq*Gr),1))
legend('Sintonizacion','Diseño LGR con SMith')
title(['Comparación entre LGR y Sintonizacion'])