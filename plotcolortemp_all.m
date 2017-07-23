close all
mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=5*msun;
e=1e51;
r=0.2*rsun;
Z=0.005;
a=7.566e-15;
X=0.7;
kappa_T=0.2*(1+X);
kappa=kappa_T;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;
%load('/home/nilou/Data/processeddata/ic/colortemp_phi0_v0.mat','t_color_phi0')

t_color=load('/home/nilou/Data/processeddata/ic/colortemp_tot_v0.mat','t_color_tot');
t_color_ic=t_color.t_color_tot;
t_color=load('/home/nilou/Data/processeddata/RSG/colortemp_tot_v3.mat','t_color_tot');
t_color_rsg=t_color.t_color_tot;
t_color=load('/home/nilou/Data/processeddata/BSG/colortemp_tot.mat','t_color_tot');
t_color_bsg=t_color.t_color_tot;
a=get(gcf,'Position');
x0=15;
y0=15;
width=650;
height=650;


myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');

set(myFigure,'units','points','position',[x0,y0,width,height])   


time=load('/home/nilou/Data/timesteps.mat');t_color_bsg
time_bsg=(time.time1*tconv);


plot(time.time1(23:32),log10(t_color_ic(1,23:32)/11604.52),'k','LineWidth',2);hold on
plot(time.time1(23:42),log10(t_color_bsg(1,23:42)/11604.52),'b','LineWidth',2);hold on
plot(time.time1(24:47),log10(t_color_rsg(1,24:47)/11604.52),'r','LineWidth',2);hold on

xlabel('t/t_*'); 
ylabel('Log (T_{color} [eV])');

%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortemp_all_v0.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')