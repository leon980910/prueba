clc
clear all
close all
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
cut_data_i = 48;
cut_data_f = 6;
t = ti(1:end-cut_data_i-cut_data_f);
y = pi(cut_data_i:end-cut_data_f-1);
u = ui(cut_data_i:end-cut_data_f-1);
%--    
% ________________________________________________________
% Dise?o
s = tf('s');
vector = [t,y];
% calculo de variables
% dy = average(y);
dy = max(y); %Mejor Opcion
du = max(u);
% du = max(u) - min(u);

dy25 = dy*0.25;
dy50 = dy*0.50;
dy75 = dy*0.75;
% -----------------------------
% interpolar
dt25 = CalcularPunto(vector,dy25)
dt50 = CalcularPunto(vector,dy50)
dt75 = CalcularPunto(vector,dy75)
kp = dy / du;
tao =  0.9102*(dt75 - dt25)
tm = 1.2620*dt25 - 0.2620*dt75
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
Gp1 = kp*exp(-tm*s)/(tao*s+1)
yi = lsim(Gp1,u,t);
plot(t,yi)
hold on
plot(t,y,'r')
legend('Orden 1','Planta','Location','SouthEast')
% ------------------------------------------------------------------------
% Calculo del indice de aproximacion
% definir un tiempo para ambas se?ales
data_error = y-yi
error = sum(((y-yi).^2))


