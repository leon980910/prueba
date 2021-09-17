function [KP a_f Cs]=rootlocusflujo(G,ts,mp1)
%controlador por raices
KP = [];
a_f = [];
s = tf('s');
G1=zpk(G);
[N,D] = tfdata(G,'v');
%parametros de diseño
mp  = mp1;
tss = ts;

theta = 4/tss;
wd = (-pi*theta)/(log(mp/100));

%polo de diseño 
pd = -theta + wd*j;



%diseño de polos en polo de diseños
dp = 1/(s + theta - wd*j);
%plotear polos
figure
pzmap(G*1/s*dp)

%polo de la funcion trans
[p,z] = pzmap(G);

%polo del controlador 
Kp1 = 1;
a1 = 1;
Cs1 = Kp1 * (s + a1)/s; 
Cp = pole(Cs1);

% sumatoria de polos

angulo_Cs1 = 180 - (radtodeg(atan(abs(imag(pd))/abs(real(pd)))));
angulo_Gsf1 =180-(radtodeg(atan(abs(imag(pd))/(abs(pd(1,1))-abs(real(p))))));


angula_a = 180 - (angulo_Cs1 + angulo_Gsf1 );

a_f = abs(imag(pd))/(tan(degtorad(angula_a)))+ (abs(real(pd)));

%hallar k


vectp2 = [imag(pd) (real(pd)-p(1,1))];
vectz1 = [ imag(pd) (a_f - abs(real(pd)))];


prodepolos = norm(vectp2) ; 
prodezeros = norm(vectz1) ;

KP = prodepolos / (prodezeros * N(1,2));
%Kpr = 0.12416;



figure
Cs = KP * ( s + a_f)/s;
step(feedback(Cs*G,1))
end