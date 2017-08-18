mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
ttot=rtot/sqrt(etot/mtot);
msun=1.989e33;
rsun=6.955e10;





Z=0.005;
a=7.566e-15;
X=0.7;
kappa_T=0.2*(1+X);


m=15*msun;
e=1e51;
r=49*rsun;



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
tconv_bsg=rconv/vconv;

% load('/home/nilou/Data/processeddata/BSG/luminosity_0.mat','luminosity');
% load('/home/nilou/Data/processeddata/BSG/luminosity_tot.mat','luminosity_tot');
% load('/home/nilou/Data/processeddata/BSG/luminosity_90.mat','luminosity90');
 load('/home/nilou/Data/processeddata/BSG/luminosityV1_0.mat','luminosity')
 load('/home/nilou/Data/processeddata/BSG/luminosityV1_90.mat','luminosity90')
 load('/home/nilou/Data/processeddata/BSG/luminosityV1_tot.mat','luminosity_tot')
time=load('/home/nilou/Data/processeddata/timesteps.mat');

%  const=linspace(0.9,0.5,16);
%  luminosityRSG.luminosity(35:50)=luminosityRSG.luminosity(35:50).*const(1:16);
% 
% const=linspace(0.8,0.6,7);
% luminosityBSG.luminosity(35:41)=luminosityBSG.luminosity(35:41).*const(1:7);


%luminosityRSG.luminosity(23:24)=luminosityRSG1.luminosity(23:24);
%close all
tt=1:50;
% time_rsg=(time.time1*tconv_rsg);
 time_bsg=(time.time1*tconv_bsg);
%time_ic=(time.time1*tconv_ic);


m=15*msun;
e=1e51;
r=49*rsun;

t_sph=time_bsg(23:41)-time_bsg(23)+113.59;
t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
L_s=zeros(length(time_bsg(23:41)),1);
L_s(t_sph<t_s)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34*(t_sph(t_sph<t_s)./60).^(-4/3);
L_s(t_sph>t_s)=2*10^42*(m/(15*msun))^-0.73*(r/(50*rsun))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(3600)).^-0.35;
%L_s(t_sph>t_s)=3.3*10^42*(m/(msun))^-0.74*(r/(1e12))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(24*3600)).^-0.34;

save('/home/nilou/Data/processeddata/BSG/sphericallum.mat','L_s');

% t_s
% L_s
close all

a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])


f= fit(log10(time_bsg(23:40)),log10(1*luminosity(23:40))','rat44');
luminosity_fit=f(log10(time_bsg(23:40)'));
plot(log10(time_bsg(23:40)),log10(1*luminosity(23:40)),'-r','LineWidth',1.5),hold on,

f= fit(log10(time_bsg(23:40)),log10(2*luminosity90(23:40))','rat44');
luminosity90_fit=f(log10(time_bsg(23:40)'));
plot(log10(time_bsg(23:40)),log10(2*luminosity90(23:40)),'--k','LineWidth',1.5),hold on,

f= fit(log10(time_bsg(23:40)),log10(2*luminosity_tot(23:40))','rat44');
luminosity_tot_fit=f(log10(time_bsg(23:40)'));
plot(log10(time_bsg(23:40)),log10(2*luminosity_tot(23:40)),'-.b','LineWidth',1.5),hold on,

 
plot(log10(time_bsg(23:41)),log10(L_s)',':m','LineWidth',1.5),hold on,


set(gca,'LineWidth',2,'FontSize',12);
%loglog(time_bsg(1:50),4*luminosityBSG.luminosity(1:50),'b','Linewidth',1),hold on,
%loglog(time_ic(1:50),4*luminosityIc.luminosity,'g','Linewidth',1),hold on,
%axis([10^-1,10^6,10^39,10^48])
xlabel('Log (t [sec])'); 
ylabel('Log (L_{tot} [erg/s])');
%title('Time of Maximum Compression Rate');
legend('\Theta=0','\Theta=\pi/2','Total','Spherical','Location','northeast');
name=['/home/nilou/Data/lum_BSGV1.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
