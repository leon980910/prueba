x=[0.97 1.12 2.92 3.00 3.33 3.97 6.10 8.39 8.56 9.44];
y=[2.58 0.43 0.06 5.74 7.44 8.07 6.37 2.51 1.44 0.52];
xx=[1.0 2.0 3.5 5.5 8.0];
yy=interp1(x,y,xx,'linear');
disp([xx' yy'])
hold on
plot(x,y,'-bo','markersize',3,'markerfacecolor','b')
plot(xx,yy,'ro','markersize',4,'markerfacecolor','r')
xlabel('x')
ylabel('y')
grid on
title('Interpolación lineal');
hold off