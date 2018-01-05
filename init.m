close all

x=linspace(0,pi,1000);
figure
plot(x, cos(2*x)); hold on
plot(x, (abs(cos(x)).^0.8))