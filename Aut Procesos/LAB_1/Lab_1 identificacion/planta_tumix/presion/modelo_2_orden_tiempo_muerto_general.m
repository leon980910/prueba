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
cut_data_i = 40;
cut_data_f = 6;
t = ti(1:end-cut_data_i-cut_data_f);
y = pi(cut_data_i:end-cut_data_f-1);
u = ui(cut_data_i:end-cut_data_f-1);
%--    
% ________________________________________________________
% Diseño
s = tf('s');
vector = [t,y];
% calculo de variables
% dy = max(y) - min(y);
dy = average(y(30:end));
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
a = (-0.6240*dt25 + 0.9866*dt50-0.3626*dt75)/(0.3533*dt25-0.7036*dt50 + 0.3503*dt75);
tao2 =   (dt75-dt25)/(0.9866 + 0.7036*a);%???
t1 = tao2;
t2 = tao2*a;
tm2 = dt75-(1.3421 + 1.3455*a)*tao2;


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
Gp1 = (kp*exp(-tm2*s))/((tao2*s + 1)*(a*tao2*s + 1)) ;

yi = lsim(Gp1,u,t);
plot(t,yi)
hold on
plot(t,y,'r')
legend('General','Planta','Location','SouthEast')
% ------------------------------------------------------------------------
% Calculo del indice de aproximacion
% definir un tiempo para ambas señales
data_error = y-yi
error = sum(((y-yi).^2))
