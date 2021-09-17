clc
clear all
close all
s = tf('s');

Gp = (9*exp(-20*s))/(1 + 40*s);

Gd1 = (20*exp(-3*s))/(1 + 30*s);
Gd2 = (exp(-20*s))/(1 + 60*s);

figure
step(Gp)

ts=105;
mp1 =6.34;
G=(9)/(1 + 40*s);

[KP a_f Cs]=rootlocusflujo(G,ts,mp1)


[N6,D6]=tfdata(Cs,'v')
figure
step(feedback(Cs*G,1))
title('Controlador PI para planta sin tiempo muerto')

%% CONTROLADOR SMITH
tmF=20;
KpF=9;
taoF=40;
Gps=(9)/(1 + 40*s);
[N,D]=pade(tmF,2);
ret = tf(N,D);
Gr = (KpF*ret)/((taoF*s)+1);
CL = feedback(Cs,Gps)
Ceq = feedback(CL,-Gr)

figure
step(feedback(Ceq*Gr,1),200)
title('Controlador SMITH')

%% Perturbaciones

Gff1=-((20+40*s)/(9+30*s)) %-3
Gff2 = -((1+40*s)/(9+60*s)) %-20

%Planta
[Ngr,Dgr]=tfdata(Gr,'v');

%Smith
[Nsmith,Dsmith]=tfdata(Ceq,'v');

% Control perturbacion:
[Nf1,Df1]=tfdata(Gff1,'v')
[Nf2,Df2]=tfdata(Gff2,'v')

% Valvula:
Gv = (0.1*exp(-1*s))/(1+3*s)
[N,D]=pade(1,2);
ret = tf(N,D);
Gvr = (0.1*ret)/((3*s)+1);
[Nf3,Df3]=tfdata(Gvr,'v')


%Perturbacion 1:
Gd1 = (20*exp(-3*s))/(1 + 30*s);
[N,D]=pade(3,2);
ret = tf(N,D);
Gd1s = (20*ret)/((30*s)+1);
[Nf4,Df4]=tfdata(Gd1s,'v')

%Perturbacion 2:
Gd2 = (exp(-20*s))/(1 + 60*s);
[N,D]=pade(20,2);
ret = tf(N,D);
Gd2s = (ret)/((60*s)+1);
[Nf5,Df5]=tfdata(Gd2s,'v')



