close all
% mtot=0.0048187313;
% etot=0.0130949665;
% rtot=0.5;
% msun=1.989e33;
% rsun=6.955e10;
% m=14*msun;
% e=1e51;
% r=400*rsun;
% Z=0.005;
% a=7.566e-15;
% X=0.7;
% kappa_T=0.2*(1+X);
% kappa=kappa_T;
% c=3e10;
% rconv=r/rtot;
% econv=e/etot;
% mconv=m/mtot;
% 
% pconv=econv/rconv^3;
% rhoconv=mconv/rconv^3;
% vconv=sqrt(econv/mconv);
% tconv=rconv/vconv;

mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=14*msun;
e=1e51;
r=400*rsun;
kappa=0.34;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

cita=0;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end



load([path '/processeddata/RSG/colortem_1024_tot_v2.mat'],'t_tot')
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])
time=load('/home/nilou/Data/processeddata/timesteps_1024.mat');
time_rsg=(time.time1*tconv);
% 

f= fit(log10(time_rsg(173:1:420)),log10(t_tot(173:420)/(1.5*11604.52))','poly3');
t_t=f(log10(time_rsg(173:1:420)'));
plot(log10(time_rsg(173:1:420)),t_t,'-k','LineWidth',1.5),hold on,

time_rsg(420)/tconv




mtot=9.651994e-3;
etot=2.997521e-2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=14*msun;
e=1e51;
r=400*rsun;
kappa=0.34;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;


load('/home/nilou/Data/Spherical/processeddata/RSG/colortemp_RSG_spherical.mat','t_loc_all','time_rsg','etha_all');
time_rsg1=time_rsg;

t_sph=time_rsg1(141:270)-time_rsg1(141)+185;
T_c=zeros(length(t_sph),1);
t_s=14*3600* (m/(15*msun))^0.43*(r/(500*rsun))^1.26*(e/1e51)^-0.56;
T_c(t_sph<t_s)=10*(m/(15*msun))^-0.22*(r/(500*rsun))^0.12 *(e/1e51)^0.23* (t_sph(t_sph<t_s)/3600).^-0.36;
T_c(t_sph>t_s )=3*(m/(15*msun))^-0.13*(r/(500*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t_s)/(24*3600)).^-0.56;


f= fit(log10(time_rsg1(141:270)-time_rsg1(109)+time_rsg1(65)),log10(t_loc_all(141:270)/(1.5*11604.52)),'poly5');
t_t=f(log10(time_rsg1(141:270)-time_rsg1(109)+time_rsg1(65)));
%plot(log10(time_rsg1(141:270)-time_rsg1(109)+25000),(log10(T_c)+0.5*t_t)/1.5,'--k','LineWidth',1.5),hold on,


T_c_RSG=(log10(T_c)+0.3*t_t)/1.3;
time_RSG_log=log10(time_rsg1(141:270)-time_rsg1(109)+25000);
T_c_RSG(72:87)=[];
time_RSG_log(72:87)=[];
plot(time_RSG_log,T_c_RSG,'--k','LineWidth',1.5),hold on,
plot(log10(time_rsg1(141:270)-time_rsg1(109)+25000),log10(T_c),'-.k','LineWidth',2);

save('/home/nilou/Data/Spherical/processeddata/RSG/T_c_RSG_spherical.mat','T_c_RSG','time_RSG_log','time_rsg1');


xlabel('Log (t [sec])'); 
ylabel('Log (T_{color} [eV])');
legend('Total','Spherical','N&S 2010')
%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortempRSG_1024_v1.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')