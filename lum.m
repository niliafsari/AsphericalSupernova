mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=15*msun;
e=1e51;
r=49*rsun;
kappa=0.34;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

time=load('/home/nilou/Data/timesteps.mat');
  
pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;
times=(time.time1*tconv);
load('luminosity.mat', 'luminosity')
close all
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])
loglog(times(1:50),luminosity,'b','Linewidth',1),hold on,
xlabel('time (sec)'); 
ylabel('L_{tot} (erg/s)');
%title('Time of Maximum Compression Rate');
%legend('v_\phi','fitted','Location','northeast');
set(gca,'LineWidth',1.5,'FontSize',10);
name=['/home/nilou/Data/totlum'];
print('-dpdf',name) 
export_fig(name, '-pdf')