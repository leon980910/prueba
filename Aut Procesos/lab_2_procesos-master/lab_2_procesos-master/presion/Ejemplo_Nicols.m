% Ejemplo de Nicols 1 m�todo
% Tomado de 
% https://intranet.ceautomatica.es/old/actividades/jornadas/XXIX/pdf/205.pdf
ng = [1];
dg = [1 6 5 0];
G = tf(ng,dg);
clf,rlocus(G);
[km,pole] = rlocfind(G)
wm = max(imag(pole));
kp = 0.6*km;
kd = (kp*pi)/(4*wm);
ki = (kp*wm)/pi; 
nk = [kd kp ki];
dk = [1 0];
gc = tf(nk,dk)
gd = series(G,gc)
GT = feedback(gd,1)
step(GT,'r') %Escalon rojo para la funci�n de transferencia con el regulador
hold on
GS=feedback(G,1)
step(GS,'g') %Escalon Verde para la funcion de transferencia sola
title(['Nicole PID sintonizaci�n'])
legend('PID','Plant','Location','SouthEast')