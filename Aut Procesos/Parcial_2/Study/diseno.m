% diseno PI
clc
clear all
close all
%%
% sistema
s = tf('s');
G = 10/(5*s+6);
Gp1 = G*exp(-0.5*s);
ts = 5; 
Mp = 0.1 ; %%overshoot 10% %0.01 1%
sigma = 4/ts;
wd = -pi*sigma/(log(Mp));
%%
% Calculo del punto a compensado
% Notas importantes p1 son negativos arriba positivos
% si los polos son menores a 
p1 = 1.2;
syms a ;
eqn = atan(abs(wd)/(a-sigma)) - (pi-atan(abs(wd)/(sigma-p1)))  - (pi-atan(abs(wd)/(sigma))) == -pi;
a = double(solve(eqn,a))


% Calculo de la ganacia del controlador
syms kd ;
eqn_1 = kd*(2)*(sqrt(wd*wd+(a-sigma)*(a-sigma)))/((sqrt(wd*wd+sigma*sigma))*(sqrt(wd*wd+(p1-sigma)*(p1-sigma)))) == 0.8 ;
kd = double(solve(eqn_1,kd))

%%estrucutra del PI
s = tf('s');
PI = kd*(s+a)/s;
step(G)
hold on
con_control = feedback(PI*G,1);
step(con_control)

%%Diseño del controlador SMITH
G = 10/(5*s+6);
Gp1 = G*exp(-0.5*s);
pi = PI;
P = pade(Gp1,1);
C = feedback(pi,G);
Ceq = feedback(C,-P);
SMITH = feedback(Ceq*G,1);
step(SMITH)
title('Controlador SMITH')
%%
%Discretizar
T = 0.2;%cambio del tiempo de muestreo
ceq_d = c2d(ceq,T)
Gz = c2d(P,T)
Gcz =c2d(G,T)
G_d =c2d(G,T)

%%ecuaciones de diferencia
[NGz,DGz]=tfdata(Gz,'v')
[NGcz,DGcz]=tfdata(ceq_d,'v')

U=ones(1,100);
U0=0; U1=0; U2=0; U3=0; U4=0; U5=0; U6=0; U7=0; U8=0; 
Y1=0; Y2=0; Y3=0; Y4=0; Y5=0; Y6=0;
E1=0; E2=0; E3=0; E4=0; E5=0; E6=0; E7=0; E8=0;
ref=1;
for k = 9:200
%     Y(k) = NGz(2)*U(k)+NGz(3)*U(k-1)-...
%         DGz(2)*Y(k-1)-DGz(3)*Y(k-2);
%     E(k) = ref-Y(k);
%     U(k) = NGcz(1)*E(k)+NGcz(2)*E(k-1)+NGcz(3)*E(k-2)+NGcz(4)*E(k-3)+NGcz(5)*E(k-4)-...
%         DGcz(2)*U(k-1)-DGcz(3)*U(k-2)-DGcz(4)*U(k-3)-DGcz(5)*U(k-4);
    
    Y0 = NGz(2)*U0+NGz(3)*U1-...
        DGz(2)*Y1-DGz(3)*Y2;
    E0 = ref-Y0;
    U0 = NGcz(1)*E0+NGcz(2)*E1+NGcz(3)*E2+NGcz(4)*E3+NGcz(5)*E4-...
        DGcz(2)*U1-DGcz(3)*U2-DGcz(4)*U3-DGcz(5)*U4;
    
    salida(k)=Y0;
    s_c(k) = U0;
    
    U4=U3; U3=U2; U2=U1; U1=U0;
    Y2=Y1; Y1=Y0;
    E4=E3; E3=E2; E2=E1; E1=E0;   
end
figure
plot(salida)
title(['Control Smith Ecuaciones de Diferencia'])
figure
plot(s_c)
title(['Salida del controlador'])


% Calculo de los errores
clear s
syms s

G_i=10*exp(-0.5*s)/(5*s+6);
Error_Planta=double(limit((s/(1+G_i))*(1/s)))

P=(0.7333*s + 1.303)/(s);
Error_PID=double(limit((s/(1+(G_i*P)))*(1/s^2)))





