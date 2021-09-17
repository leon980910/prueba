%% PARCIAL TEORICO

clc
clear all
close all
s = tf('s');

%% Planta del sistema lazo abierto
Vs = (10)*exp(-0.5*s)/(5*s+6); 
Vs1 = (10)/(5*s+6) 

%% Diseño de controlador PI

tss = 3;
Ov = 0.10;

Z = abs(log(Ov))/sqrt((pi^2)+((log(Ov))^2));
Wn = 3/(Z*tss);

Pr = -Z*Wn; %Parte real 
Pi =(Wn*sqrt(1-(Z^2))); %Parte imaginaria
%Pd = Pr+j*Pi   %Polo deseado

%% Comprobar si el polo deseado pertenece a LGR
a = Pr+i*Pi;
Ang = ((10)/(5*a+6)) *(1/a);
[theta,rho]=cart2pol(real(Ang),imag(Ang));
theta = theta*(180/pi)

% theta no es multiplo de 180°, se agrega un cero de compensación.
if theta<180
theta1 = 180-theta  %Angulo del cero aportado al LGR
cero = Pi/tand(theta1); %Ubicacion del cero.
end
C1 = (s+cero)/s %Controlador sin ganancia

%Ganancia Kp del controlador PI
Kp = 1/(abs(((10)/((5*a)+6))*((a+cero)/a)))

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

tmF = 0.5;
KpF=10;
taoF=5;

[N1,D1]=pade(tmF,2);
ret = tf(N1,D1);
Gr = (KpF*ret)/((taoF*s)+6);
CL = feedback(Cs,Vs1)
Ceq = feedback(CL,-Gr)

Gper = 5/(5*s+6);

figure
step(feedback((Ceq*Gr)-(Gper*0.3),1))
% hold on
% step(-(Gper*0.3),3)
title('Controlador SMITH para flujo')

%% Discretización

Ts = 1.83/10;
Gz = c2d(Gr,Ts,'tustin')
Gcz = c2d(Ceq,Ts,'tustin')

figure
step(feedback(Gcz*Gz,1))
title('Controlador SMITH Discreto')

[NGz,DGz]=tfdata(Gz,'v')
[NGcz,DGcz]=tfdata(Gcz,'v')

%% Implementacion PLC

U=zeros(1,100);
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
    Uaux(k)=U0;
    U5=U4; U4=U3; U3=U2; U2=U1; U1=U0;
    Y3=Y2; Y2=Y1; Y1=Y0;
    E5=E4; E4=E3; E3=E2; E2=E1; E1=E0;   
    
end

figure
plot(Y,'r')
hold on
plot(Aux,'b')
title('Implementacion en PLC')
xlabel('tiempo(s)')
ylabel('vsicosidad')
legend('Ciclo for','Constantes')

%% Error en estado estacionario
