clc 
clear all
close all
s = tf('s');

Gp = ((9)*exp(-20*s))/((40*s)+1);

V = ((0.1)*exp(-1*s))/((3*s)+1);

Greal = Gp*V
Vs1 = ((0.9))/((120*s^2)+(43*s)+1);
%% Diseño de controlador PI

tss = 90;
Ov = 0.05;

Z = abs(log(Ov))/sqrt((pi^2)+((log(Ov))^2));
Wn = 3/(Z*tss);

Pr = -Z*Wn; %Parte real 
Pi =(Wn*sqrt(1-(Z^2))); %Parte imaginaria
%Pd = Pr+j*Pi   %Polo deseado

%% Comprobar si el polo deseado pertenece a LGR
a = Pr+i*Pi;
Ang = ((0.9)/((120*a^2)+(43*a)+1)) *(1/a);
[theta,rho]=cart2pol(real(Ang),imag(Ang));
theta = theta*(180/pi)

% theta no es multiplo de 180°, se agrega un cero de compensación.
if theta<180
theta1 = 180-theta  %Angulo del cero aportado al LGR
cero = Pi/tand(theta1); %Ubicacion del cero.
end
C1 = (s+cero)/s %Controlador sin ganancia

%Ganancia Kp del controlador PI
Kp = 1/(abs(((0.9)/((120*a^2)+(43*a)+1))*((a+cero)/a)))

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

Gr = pade(Greal,2);
CL = feedback(Cs,Vs1)
Ceq = feedback(CL,-Gr)

figure
step(feedback((Ceq*Gr),1))
% hold on
% step(-(Gper*0.3),3)
title('Controlador SMITH')
