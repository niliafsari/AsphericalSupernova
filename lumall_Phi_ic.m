mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
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

time=load('/home/nilou/Data/timesteps.mat');

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

f= fit(log10(time_ic(23:40)),log10(luminosity(23:40)/9.7)','rat44');
luminosity_fit=f(log10(time_ic(23:32)'));
plot(log10(time_ic(23:32)),luminosity_fit,'-r','LineWidth',1.5),hold on,

f= fit(log10(time_ic(23:40)),log10(2*luminosity90(23:40)/9.7)','rat44');
luminosity90_fit=f(log10(time_ic(23:32)'));
plot(log10(time_ic(23:32)),luminosity90_fit,'--k','LineWidth',1.5),hold on,

f= fit(log10(time_ic(23:40)),log10(2*luminosity_tot(23:40)/9.7)','rat44');
luminosity_tot_fit=f(log10(time_ic(23:32)'));
plot(log10(time_ic(23:32)),luminosity_tot_fit,'-.b','LineWidth',1.5),hold on,

 
plot(log10(time_ic(23:32)),log10(L_s)',':m','LineWidth',1.5),hold on,


xlabel('Log (t [sec])'); 
ylabel('Log (L_{tot} [erg/s])');

legend('\Theta=0','\Theta=\pi/2','Total','Spherical','Location','northeast');
set(gca,'LineWidth',2,'FontSize',12);
name=['/home/nilou/Data/lum_ic.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
