clc
clear all
close all 
% Data
clear all
close all
clc

load Datos

ti = Data(:,1);
fi = Data(:,2);
pi = Data(:,3);
ui = Data(:,4);
% plot(t,f)
% figure
% plot(t,p)
% figure
% plot(t,u)
%recorte de data
cut_data_i = 34;
cut_data_f = 6;
t = ti(1:end-cut_data_i-cut_data_f);
y = fi(cut_data_i:end-cut_data_f-1);
u = ui(cut_data_i:end-cut_data_f-1);
%--    
% ________________________________________________________
% Diseño
s = tf('s');
vector = [t,y];
% calculo de variables
dy = max(y) - min(y);
du = max(u);
dy25 = dy*0.25;
dy50 = dy*0.50;
dy75 = dy*0.75;
% -----------------------------
% interpolar
dt25 = CalcularPunto(vector,dy25)
dt50 = CalcularPunto(vector,dy50)
dt75 = CalcularPunto(vector,dy75)
kp = dy / du;
tao1 = 0.5776*(dt75-dt25);%t'
tm1 = 1.5552*dt25-0.5552*dt75;
tm2 = tm1;
a = (dt50-tm1-1.4362*tao1)/(1.9844*tao1-dt50+tm1);
tao2 = 2*tao1/(1+a);%t´´
t1 = tao2;
t2 = tao2*a;
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
Gp1 = (kp*exp(-tm2*s))/((tao2*s + 1)*(a*tao2*s + 1)) ;
yi = lsim(Gp1,u,t);
plot(t,yi)
hold on
plot(t,y,'r')
legend('Simplificado','Planta','Location','SouthEast')
% ------------------------------------------------------------------------
% Calculo del indice de aproximacion
% definir un tiempo para ambas señales
data_error = y-yi
error = sum(((y-yi).^2))