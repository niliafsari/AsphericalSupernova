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

load('/home/nilou/Data/processeddata/ic/colortemp_tot_v1.mat','t_color_bsg');
%load('/home/nilou/Data/processeddata/ic/colortemp_phi0_v0.mat','t_color_phi0')

%load('/home/nilou/Data/processeddata/BSG/colortemp_tot.mat','t_color_bsg');

%t_color_ic=t_color.t_color_tot;
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

time=load('/home/nilou/Data/processeddata/timesteps.mat');
time_ic1=(time.time1*tconv);

% T_c=zeros(10,1);
% t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
% t1=200*(m/(15*msun))^-0.24*(r/(5*rsun))^0.94*(e/1e51)^0.29;
% t2=400*(m/(15*msun))^-0.77*(r/(5*rsun))^0.62*(e/1e51)^0.99;
% 
% t_sph=time_ic(23:32)-time_ic(23)+0.46;
% T_c(t_sph<t_s)=10^3*(m/(15*msun))^-1.5*(r/(5*rsun))^-0.2 *(e/1e51)^1.4* t_sph(t_sph<t_s).^-0.4;
% T_c(t_sph>t_s & t_sph<t1)=140*(m/(15*msun))^-1.2*(r/(5*rsun))^-0.9 *(e/1e51)^1.7* (t_sph(t_sph>t_s & t_sph<t1)/t_s).^-2.2;
% T_c(t_sph>t1 & t_sph<t2)=40*(m/(15*msun))^0.05*(r/(5*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/60).^-0.4;
% T_c(t_sph>t2)=5*(m/(15*msun))^-0.11*(r/(5*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;
% 
% 


load('/home/nilou/Data/processeddata/ic/colortemp_ic_spherical.mat','t_loc_all','time_ic','etha_all');


time_ic([113 114 115])=[];
t_loc_all([113 114 115])=[];

plot(log10(time_ic1(23:32)),log10(t_color_bsg(1,23:32)*5/11604.52),'k','LineWidth',2); hold on

f= fit(log10(time_ic(105:140)-time_ic(105)+time_ic(66)),log10(t_loc_all(105:140)/11604.52),'rat34');
t_t=f(log10(time_ic(105:137)-time_ic(105)+time_ic(66)));
plot(log10(time_ic(105:137)-time_ic(105)+time_ic(66)),t_t,'--k','LineWidth',1.5),hold on,

% plot(log10(time_ic(105:140)-time_ic(105)+time_ic(66)),log10(t_loc_all(105:140)/11604.52),'--r','LineWidth',1.5),hold on,


%plot(log10(time_ic(105:140)-time_ic(105)+time_ic(66)),log10(t_loc_all(105:140)/11604.52),'--k','LineWidth',2);
xlabel('Log (t [sec])'); 
ylabel('Log (T_{color} [eV])');
legend('Total','Spherical')
%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortempic_1024.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')