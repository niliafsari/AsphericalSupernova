m=4.825996760367126832e-3*2;
e=9.199700591776255711e-4*2;
r=0.5;
rho_c=m/r^3;
v_c=r/sqrt(e/m);
t=1;
name=['/home/nilou/poly_gamma133.csv'] ;
data= csvread(name)
dens1=data(:,2);
vel1=data(:,4);
pres1=data(:,3);
radius=data(:,1);



a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
close all
% 
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])
y = logspace(-16,2,50);
rad=ones(1,50)*0.75;
plot(radius/0.5,log10(((radius/0.5).^2).*dens1),'r');hold on
plot(radius/0.5,log10(vel1),'m');hold on
plot(radius/0.5,log10(pres1),'g');hold on


xlabel('r/R_*')
ylabel('log r^3 \rho/\rho_* ')
legend('\rho r^3','v','p' )
plot(rad,log10(y),'k');hold on
axis([0 2 -16 2])
