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

load('/home/nilou/Data/processeddata/RSG/colortemp_RSG_spherical.mat','t_loc_all','time_rsg','etha_all');
time_rsg1=time_rsg;

% time=load('/home/nilou/Data/processeddata/timesteps.mat');
% time_rsg1=(time.time1*tconv);

time=load('/home/nilou/Data/processeddata/timesteps_1024.mat');
time_rsg=(time.time1*tconv);
% 
% t_sph=time_rsg1(23:50)-time_rsg1(23)+927.3;
% T_c=zeros(50-23+1,1);
% t_s=14*3600* (m/(15*msun))^0.43*(r/(500*rsun))^1.26*(e/1e51)^-0.56;
% T_c(t_sph<t_s)=10*(m/(15*msun))^-0.22*(r/(500*rsun))^0.12 *(e/1e51)^0.23* (t_sph(t_sph<t_s)/3600).^-0.36;
% T_c(t_sph>t_s )=3*(m/(15*msun))^-0.13*(r/(500*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t_s)/(24*3600)).^-0.56;

f= fit(log10(time_rsg(173:1:420)),log10(t_tot(173:420)/11604.52)','poly3');
t_t=f(log10(time_rsg(173:1:420)'));
plot(log10(time_rsg(173:1:420)),t_t,'-k','LineWidth',1.5),hold on,

f= fit(log10(time_rsg1(109:210)-time_rsg1(109)+time_rsg1(65)+927.3),log10(t_loc_all(109:210)/11604.52),'poly5');
t_t=f(log10(time_rsg1(109:210)-time_rsg1(109)+time_rsg1(65)+927.3));
plot(log10(time_rsg1(109:210)-time_rsg1(109)+time_rsg1(65)+927.3),t_t,'--k','LineWidth',1.5),hold on,

%plot(log10(time_rsg1(109:210)-time_rsg1(109)+time_rsg1(65)+927.3),log10(t_loc_all(109:210)/11604.52),'--k','LineWidth',2);

% plot(log10(time_rsg(175:420)),log10(t_tot(175:420)/11604.52),'k','LineWidth',2); hold on
% plot(log10(time_rsg1(23:50)),log10(T_c),'--k','LineWidth',2);
xlabel('Log (t [sec])'); 
ylabel('Log (T_{color} [eV])');
legend('Total','Spherical')
%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortempRSG_1024.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')