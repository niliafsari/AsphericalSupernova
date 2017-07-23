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

time=load('/home/nilou/Data/timesteps.mat');
cone=[ 0.02 0.02 0.1 0.1 0.1 0.15 0.17 0.17 0.2 0.2 0.22 0.23 0.23 0.24 0.25 0.27 0.3 0.32 ...
    0.33 0.34 0.34 0.34 0.34 0.34 0.34 0.35 0.36 0.36 0.36 0.37 0.37 0.38 0.39  0.4 0.4 0.41 0.42 0.42 0.42];

close all
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])


time_cone=time.time1(33:71)/ttot;


f= fit(time_cone,cone','poly3');
ff=f(time_cone);
plot(time_cone,ff,'LineWidth',2)
xlabel('t/t_*'); 
ylabel('h_{c}/R_*');



%legend('\Theta=0','\Theta=\pi/2','Total','Spherical','Location','northeast');
set(gca,'LineWidth',2,'FontSize',12);
name=['/home/nilou/Data/cone_hight.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')

