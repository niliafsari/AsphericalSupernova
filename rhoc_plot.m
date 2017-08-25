mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


ptot=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048)/rtot;

xx=zeros(2048,2048);
yy=zeros(2048,2048);

thet=zeros(2048,2048);


for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
       thet(i,j)=theta(j);
       rad(i,j)=radius(i);
    end
end
time=load('/home/nilou/Data/processeddata/timesteps.mat');
time=time.time1(:)/ttot;
cone=[0 0 0 0 0 0 0.02 0.02 0.1 0.1 0.1 0.15 0.17 0.17 0.2 0.2 0.22 0.23 0.23 0.24 0.25 0.27 0.3 0.32 ...
    0.33 0.34 0.34 0.34 0.34 0.34 0.34 0.35 0.36 0.36 0.36 0.37 0.37 0.38 0.39  0.4 0.4 0.41 0.42 0.42 0.42];
% 
% dens_c=zeros(1536,5);
% p_c=zeros(1536,5);
% j=1;
% for t=30:10:70
%     t
% 
%     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
%     dens= csvread(name)/rhotot;
%     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
%     pres= csvread(name)/ptot;
%     dens_c(:,j)=dens(513:2048,2048);
%     p_c(:,j)=pres(513:2048,2048);
%     j=j+1;
% end
% 
% close all
% a=get(gcf,'Position');
% x0=15;
% y0=15;
% width=350;
% height=300;
% myFigure = figure('PaperPositionMode','auto','Color','w');
% set(myFigure,'units','points','position',[x0,y0,width,height])
% 
% set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
% 
% 
% plot(asin(1./xx(513:2048,2048))*180/pi,log10(p_c(:,2)),'-k','LineWidth',2), hold on
% plot(asin(1./xx(513:2048,2048))*180/pi,log10(p_c(:,3)),'--k','LineWidth',2), hold on
% plot(asin(1./xx(513:2048,2048))*180/pi,log10(p_c(:,4)),':k','LineWidth',2), hold on
% plot(asin(1./xx(513:2048,2048))*180/pi,log10(p_c(:,5)),'-.k','LineWidth',2)
% 
% legend(['t/t_*=' num2str(time(40),3)],['t/t_*=' num2str(time(50),3)],['t/t_*=' num2str(time(60),3)],...
%     ['t/t_*=' num2str(time(70),3)])
% ylabel('Log (p_c/p_*)'); 
% xlabel('\theta_c [deg]');
% set(gca,'LineWidth',1.5,'FontSize',11);
% name=['/home/nilou/Data/plot/p_c.pdf'];    
% print('-dpdf',name) 
% export_fig(name, '-pdf')

% close all
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')

plot(asin(1./xx(513:2048,2048))*180/pi,log10(dens_c(:,2)),'-k','LineWidth',2), hold on
plot(asin(1./xx(513:2048,2048))*180/pi,log10(dens_c(:,3)),'--k','LineWidth',2), hold on
plot(asin(1./xx(513:2048,2048))*180/pi,log10(dens_c(:,4)),':k','LineWidth',2), hold on
plot(asin(1./xx(513:2048,2048))*180/pi,log10(dens_c(:,5)),'-.k','LineWidth',2)

legend(['t/t_*=' num2str(time(40),3)],['t/t_*=' num2str(time(50),3)],['t/t_*=' num2str(time(60),3)],...
    ['t/t_*=' num2str(time(70),3)])
ylabel('Log (\rho_c/\rho_*)'); 
xlabel('\theta_c [deg]'); 
set(gca,'LineWidth',1.5,'FontSize',11);
name=['/home/nilou/Data/plot/rho_c.pdf'];    
print('-dpdf',name) 
export_fig(name, '-pdf')