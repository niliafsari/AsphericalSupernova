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

m=14*msun;
e=1e51;
r=400*rsun;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;
v=sqrt(e/m);
rho=m/(r^3);
Z=0.005;
a=7.566e-15;
X=0.7;
kappa_T=0.2*(1+X);

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv_rsg=rconv/vconv;

% luminosity_phi0=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55_Phi0.mat','luminosity');
% luminosity_tot=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55_temp.mat','luminosity');
% luminosity_phi90=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55_Phi90.mat','luminosity');

load('/home/nilou/Data/processeddata/RSG/luminosity_0.mat','luminosity')
load('/home/nilou/Data/processeddata/RSG/luminosity_90.mat','luminosity90')
load('/home/nilou/Data/processeddata/RSG/luminosity_tot.mat','luminosity_tot')

time=load('/home/nilou/Data/processeddata/timesteps.mat');

%    const=linspace(0.8,0.6,16);
%    luminosityRSG.luminosity(40:50)=luminosityRSG.luminosity(40:50).*const(1:11);
% 
% const=linspace(0.8,0.6,7);
% luminosityBSG.luminosity(35:41)=luminosityBSG.luminosity(35:41).*const(1:7);


%luminosityRSG.luminosity(23:24)=luminosityRSG1.luminosity(23:24);

tt=1:50;
time_rsg=(time.time1*tconv_rsg);
%time_bsg=(time.time1*tconv_bsg);
time_ic=(time.time1*tconv_ic);

m=14*msun;
e=1e51;
r=400*rsun;

t_sph=time_rsg(23:50)-time_rsg(23)+927.3;
t_s=14*3600* (m/(15*msun))^0.43*(r/(500*rsun))^1.26*(e/1e51)^-0.56;
t_c=927.3;
L_s=zeros(length(time_rsg(23:50)),1);
%L_s(t_sph<t_c)=2e44*(m/(15*msun))^-0.37*(r/(500*rsun))^2.46*(e/1e51)^0.3*(t_c./3600).^(-4/3);
L_s(t_sph<t_s)=2e44*(m/(15*msun))^-0.37*(r/(500*rsun))^2.46*(e/1e51)^0.3*(t_sph(t_sph<t_s)./3600).^(-4/3);
L_s(t_sph>t_s)=6*10^42*(m/(15*msun))^-0.87*(r/(500*rsun))*(e/1e51)^0.96.*(t_sph(t_sph>t_s)./(24*3600)).^-0.17;
t_s
L_s
save('/home/nilou/Data/processeddata/RSG/lumspherical.mat','L_s');

close all

a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

f= fit(log10(time_rsg(23:48)),log10(luminosity(23:48))','rat44');
luminosity_fit=f(log10(time_rsg(23:48)'));
plot(log10(time_rsg(23:48)),luminosity_fit,'-r','LineWidth',1.5),hold on,

f= fit(log10(time_rsg(23:49)),log10(2*luminosity90(23:49))','rat44');
luminosity90_fit=f(log10(time_rsg(23:49)'));
plot(log10(time_rsg(23:49)),luminosity90_fit,'--k','LineWidth',1.5),hold on,

f= fit(log10(time_rsg(23:49)),log10(2*luminosity_tot(23:49))','rat44');
luminosity_tot_fit=f(log10(time_rsg(23:49)'));
plot(log10(time_rsg(23:49)),luminosity_tot_fit,'-.b','LineWidth',1.5),hold on,

 
plot(log10(time_rsg(23:50)),log10(L_s)',':m','LineWidth',1.5);

set(gca,'LineWidth',2,'FontSize',12);
%loglog(time_bsg(1:50),4*luminosityBSG.luminosity(1:50),'b','Linewidth',1),hold on,
%loglog(time_ic(1:50),4*luminosityIc.luminosity,'g','Linewidth',1),hold on,
%axis([10^-1,10^6,10^39,10^48])
xlabel('Log (t [sec])'); 
ylabel('Log (L_{tot} [erg/s])');
%title('Time of Maximum Compression Rate');
legend('\Theta=0','\Theta=\pi/2','Total','Spherical','Location','northeast');
name=['/home/nilou/Data/lum_RSG.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
