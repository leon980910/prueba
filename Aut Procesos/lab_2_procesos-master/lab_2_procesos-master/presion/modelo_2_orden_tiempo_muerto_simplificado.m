clc
clear all
close all
% Cargar los Datos
load('data_text.mat')
% plot(t,f)
% figure
% plot(t,p)
% figure
% plot(t,u)
%recorte de data
cut_data_i = 137+5;
cut_data_f = 150;
t = TIME(1:end-cut_data_i-cut_data_f)';
y = PRESION(cut_data_i:end-cut_data_f-1);
u = BOMBA(cut_data_i:end-cut_data_f-1);
% ___________________________________________________
% Diseño
s = tf('s');
vector = [t,y];
% calculo de variables
dy = average(y(30:end)) - min(average(y(1:30)));
du = max(u)-min(u);
dy25 = dy*0.25;
dy50 = dy*0.50;
dy75 = dy*0.75;
% -----------------------------
% interpolar
dt25 = CalcularPunto(vector,dy25);
dt50 = CalcularPunto(vector,dy50);
dt75 = CalcularPunto(vector,dy75);
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
yi = step(du*Gp1,t);
plot(t,yi)
hold on
plot(t,y,'r')
legend('Simplificado','Planta','Location','SouthEast')
% ------------------------------------------------------------------------
% Calculo del indice de aproximacion
% definir un tiempo para ambas señales
data_error = y-yi;
error = sum(((y-yi).^2));
num_data = size(t);
error_medio = error / num_data(1) ;

title(['Identificación orden 2 simplificado con t muerto Error: ',num2str(error_medio),' '])
legend('Orden 2 simplificado ','Planta','Location','SouthEast')
xlabel('Tiempo[s]') 
ylabel('Presión [bar]')
