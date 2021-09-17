 clc
clear all
close all 
% Data
clear all
close all
clc

% load('Datos_Modulo_2.mat');
% load('datos_mod_2_completos.mat');
load('datos_filtro.mat')

% ti = TIEMPO;
ti = [0:0.2:149.8]';
fi = FUJO;
pi = PRESION;
ui = BOMBA;

%recorte de data
cut_data_i = 324;
cut_data_f = -450+750;
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
du = max(u) ;
% du = max(u) - min(u);

dy25 = dy*0.25;
dy50 = dy*0.50;
dy75 = dy*0.75;
% -----------------------------
% interpolar
dt25 = CalcularPunto(vector,dy25)
dt50 = CalcularPunto(vector,dy50)
dt75 = CalcularPunto(vector,dy75)
k = dy / du;
tao =  0.9102*(dt75 - dt25)
tm = 1.2620*dt25 - 0.2620*dt75
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
Gp1 = k*exp(-tm*s)/(tao*s+1);
G = k/(tao*s+1);
yi = lsim(Gp1,u,t);
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
ylabel('Caudal [L/s]')
save data_lineal.mat
save planta.mat G Gp1 tm

