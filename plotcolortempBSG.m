close all
mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=15*msun;
e=1e51;
r=49*rsun;
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

% t_color=load('/home/nilou/Data/processeddata/BSG/colortemp_tot_v3.mat','t_color_tot');
% t_color_bsg=t_color.t_color_bsg;

load('/home/nilou/Data/processeddata/BSG/colortemp_tot.mat','t_color_bsg');
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height]) 


time=load('/home/nilou/Data/timesteps.mat');
time_bsg=(time.time1*tconv);

t_sph=time_bsg(23:41)-time_bsg(23)+113.59;
T_c=zeros(41-23+2,1);
t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
t1=13*60*(m/(10*msun))^-0.24*(r/(20*rsun))^0.94*(e/1e51)^0.29;
t2=20*60*(m/(10*msun))^-0.77*(r/(20*rsun))^0.62*(e/1e51)^0.99;

T_c(t_sph<t_s)=150*(m/(10*msun))^-1*(r/(20*rsun))^-0.1 *(e/1e51)^1.2* (t_sph(t_sph<t_s)/60).^-0.39;
T_c(t_sph>t_s & t_sph<t1)=70*(m/(15*msun))^-0.9*(r/(20*rsun))^-0.7 *(e/1e51)^1.1* (t_sph(t_sph>t_s & t_sph<t1)/t_s).^-2.1;
T_c(t_sph>t1 & t_sph<t2)=15*(m/(10*msun))^0.05*(r/(20*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/(15*60)).^-0.4;
T_c(t_sph>t2)=7*(m/(10*msun))^-0.11*(r/(20*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;

plot(log10(time_bsg(23:42)),log10(t_color_bsg(1,23:42)/11604.52),'k','LineWidth',2); hold on
plot(log10(time_bsg(23:42)),log10(T_c),'--k','LineWidth',2);
xlabel('Log (t [sec])'); 
ylabel('Log (T_{color} [eV])');
legend('Total','Spherical')
%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortempBSG_f.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')