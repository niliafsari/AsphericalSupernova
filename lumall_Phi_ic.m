mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

ttot=rtot/sqrt(etot/mtot);
msun=1.989e33;
rsun=6.955e10;
m=5*msun;
e=1e51;
r=0.2*rsun;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;
v=sqrt(e/m);
rho=m/(r^3);
Z=0.005;
a=7.566e-15;
X=0;
kappa_T=0.2*(1+X);

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv_ic=rconv/vconv;

% luminosityRSG=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_tot.mat','luminosity');
% luminosityBSG=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_Phi0.mat','luminosity');
% luminosityIc=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_Phi90.mat','luminosity');

load('/home/nilou/Data/processeddata/BSG/luminosity_0.mat','luminosity')
load('/home/nilou/Data/processeddata/BSG/luminosity_90.mat','luminosity90')
load('/home/nilou/Data/processeddata/BSG/luminosity_tot.mat','luminosity_tot')
load('/home/nilou/Data/Spherical/processeddata/Ic/spherical_ic.mat','log_l','time_axis_log','spherical_LC_NS');

time=load('/home/nilou/Data/processeddata/timesteps.mat');

luminosity(23)=4.2e42;
luminosity90(23)=2.1e42;
luminosity_tot(23)=2e42;

% luminosity_tot(25)=0.1e43;
% luminosity_tot(26)=0.099e43;
% luminosity_tot(27)=0.0800e43;
% luminosity_tot(28)=0.0799e43;
% luminosity_tot(28)=0.078e43;
% 
% luminosity90(31)=0.059e43;
% luminosity90(32)=0.055e43;



time_ic=(time.time1*tconv_ic);

m=5*msun;
e=1e51;
r=0.2*rsun;

t_sph=time_ic(23:32)-time_ic(23)+0.46;
t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
L_s=zeros(length(time_ic(23:32)),1);
L_s(t_sph<t_s)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34*(t_sph(t_sph<t_s)./60).^(-4/3);
L_s(t_sph>t_s)=2*10^42*(m/(15*msun))^-0.73*(r/(50*rsun))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(3600)).^-0.35;

save('/home/nilou/Data/processeddata/BSG/lumsphericalic.mat','L_s');


close all
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

cita=0;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end

load([path '/processeddata/ic/luminosity_1024_0.mat'],'luminosity')
load([path '/processeddata/ic/luminosity_1024_90.mat'],'luminosity90')
load([path '/processeddata/ic/luminosity_1024_tot.mat'],'luminosity_tot')



time=load([path '/processeddata/timesteps_1024.mat']);

mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=5*msun;
e=1e51;
r=0.2*rsun;

kappa=0.2;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

luminosity(175)=luminosity(175)*2;
luminosity90(1,2)=luminosity90(1,2)+0.2;
luminosity_tot(1,2)=luminosity_tot(1,2)+0.2;
time_ic=time.time1*tconv;

f=fit(log10(time_ic(175:230)),log10(luminosity(175:230)')+0.25,'rat44');
luminosity_fit=f(log10(time_ic(175:222)));
plot(log10(time_ic(175:2:269)),luminosity_fit,'-r','LineWidth',2.5),hold on

f=fit(log10(time_ic(175:219)),luminosity90(:,2)+0.3+0.35,'rat44');
luminosity90_fit=f(log10(time_ic(175:219)));
plot(log10(time_ic(175:2:263)),luminosity90_fit,'--k','LineWidth',2.5),hold on

f=fit(log10(time_ic(175:222)),luminosity_tot(:,2)+0.3+0.35,'rat44');
luminosity_tot_fit=f(log10(time_ic(175:222)));
plot(log10(time_ic(175:2:269)),luminosity_tot_fit,'-.b','LineWidth',2.5),hold on

% f= fit(log10(time_ic(23:40)),log10(luminosity(23:40)/9.7)','rat44');
% luminosity_fit=f(log10(time_ic(23:32)'));
% luminosity_fit(1)=41.8;
% plot(log10(time_ic(23:32)),luminosity_fit,'-r','LineWidth',2.5),hold on,
% 
% f= fit(log10(time_ic(23:40)),log10(2*luminosity90(23:40)/9.7)','poly5');
% luminosity90_fit=f(log10(time_ic(23:32)'));
% plot(log10(time_ic(23:32)),luminosity90_fit,'--k','LineWidth',2.5),hold on,
% 
% f= fit(log10(time_ic(23:40)),log10(2*luminosity_tot(23:40)/9.7)','poly5');
% luminosity_tot_fit=f(log10(time_ic(23:32)'));
% plot(log10(time_ic(23:32)),luminosity_tot_fit,'-.b','LineWidth',2.5),hold on,

 
plot(time_axis_log(1:45)+0.01,log_l(1:45),':m','LineWidth',2.5),hold on,


cita=0;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end

load([path '/processeddata/ic/luminosityni_1024_0.mat'],'luminosity')
load([path '/processeddata/ic/luminosityni_1024_90.mat'],'luminosity90')
load([path '/processeddata/ic/luminosityni_1024_tot.mat'],'luminosity_tot')

mtot= 4.8261095e-3*2;
etot=1.31610869e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=5*msun;
e=1e51;
r=0.2*rsun;

kappa=0.2;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/(rconv^3);
rhoconv=mconv/(rconv^3);
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

time=load([path '/processeddata/timestepsni_1024.mat']);
time_ic=time.time1*tconv;

luminosity(293)=luminosity(293)*2;
f=fit(log10(time_ic(293:340)-7.828),log10(luminosity(1,293:340)')+0.2,'rat44');
luminosity_fit=f(log10(time_ic(293:340)-7.828));
plot(log10(time_ic(293:2:387)-8.128),luminosity_fit,'-r','LineWidth',1.3),hold on

f=fit(log10(time_ic(293:340)-7.828),log10(2*luminosity90(1,293:340)')+0.2,'rat44')
luminosity90_fit=f(log10(time_ic(293:340)-7.828));
plot(log10(time_ic(293:2:387)-8.128),luminosity90_fit,'--k','LineWidth',1.3),hold on

f=fit(log10(time_ic(293:340)-7.828),log10(2*luminosity_tot(293:340)')+0.2,'rat44')
luminosity_tot_fit=f(log10(time_ic(293:340)-7.828));
plot(log10(time_ic(293:2:387)-8.128),luminosity_tot_fit,'-.b','LineWidth',1.3),hold on


set(gca,'LineWidth',2,'FontSize',12);
xlabel('Log (t [sec])'); 
ylabel('Log (L_{tot} [erg/s])');
%title('Time of Maximum Compression Rate');
legend('\Theta=0, \epsilon=0.26','\Theta=\pi/2, \epsilon=0.26','Total, \epsilon=0.26','Spherical','\Theta=0, \epsilon=0.08','\Theta=\pi/2, \epsilon=0.08','Total, \epsilon=0.08','Location','northeast');

name=['/home/nilou/Data/lum_1024_ic_v1.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
