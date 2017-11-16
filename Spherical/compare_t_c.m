mtot=0.0048187313*2;
etot=9.1997e-4*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=5*msun;
e=1e51;
r=0.2*rsun;
kappa=0.2;
c=3e10;
a=7.566e-15;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

time=load('/home/nilou/bgq/output/timesteps.mat')
time_ic=time.time1*tconv;

t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
t1=200*(m/(15*msun))^-0.24*(r/(5*rsun))^0.94*(e/1e51)^0.29;
t2=400*(m/(15*msun))^-0.77*(r/(5*rsun))^0.62*(e/1e51)^0.99;

t_sph=time_ic(116:143)-time_ic(116)
T_c(t_sph<t_s)=10^3*(m/(15*msun))^-1.5*(r/(5*rsun))^-0.2 *(e/1e51)^1.4* t_sph(t_sph<t_s).^-0.4;
T_c(t_sph>t_s & t_sph<t1)=140*(m/(15*msun))^-1.2*(r/(5*rsun))^-0.9 *(e/1e51)^1.7* (t_sph(t_sph>t_s & t_sph<t1)/t_s).^-2.2;
T_c(t_sph>t1 & t_sph<t2)=40*(m/(15*msun))^0.05*(r/(5*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/60).^-0.4;
T_c(t_sph>t2)=5*(m/(15*msun))^-0.11*(r/(5*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;

close all
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

plot(log10(time_ic(116:143)-time_ic(116)+time_ic(66)),log10(T_c),'LineWidth',2); hold on


load(['/home/nilou/t_c.mat'],'modeltc','output','x1','y1')
time=load('/home/nilou/Data/processeddata/timesteps_1024.mat');
time_ic_1=time.time1*tconv

t_c=modeltc([170:180 180:2:230]);
plot(log10(time_ic(116:115+length(t_c))-time_ic(116)+time_ic(66)),t_c,'LineWidth',2); hold on

xlabel('Log (t [sec])'); 
ylabel('Log (T_{color} [eV])');
legend('Aspherical Simulation','Spherical N&S2010')
%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortempic_compare.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')