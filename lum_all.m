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

% luminosityBSG=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk45.mat','luminosity');
% luminosityRSG=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed45.mat','luminosity');
% luminosityIc=load('/home/nilou/Data/processeddata/ic/luminosity_smoothed.mat','luminosity');
% time=load('/home/nilou/Data/timesteps.mat');
luminosityBSG=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_temp.mat','luminosity');
%luminosityRSG1=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed45.mat','luminosity');
luminosityRSG=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55_temp.mat','luminosity');
luminosityIc=load('/home/nilou/Data/processeddata/ic/luminosity_smoothed.mat','luminosity');
time=load('/home/nilou/Data/processeddata/timesteps.mat');

time_rsg=(time.time1*tconv_rsg);
time_bsg=(time.time1*tconv_bsg);
time_ic=(time.time1*tconv_ic);

close all
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])


const=linspace(0.9,0.5,16);
luminosityRSG.luminosity(35:50)=luminosityRSG.luminosity(35:50).*const(1:16);

const=linspace(0.8,0.6,9);
luminosityBSG.luminosity(33:41)=luminosityBSG.luminosity(33:41).*const(1:9);

%plot(log10(time_rsg(2:50)),log10(4*luminosityRSG.luminosity(2:50)),'r','Linewidth',1),hold on,
% 

% Eqn = 'a*x^b+c*x^d+e*x^f+g';
% 
% f= fit(time_rsg(2:50),4*luminosityRSG.luminosity(2:50)/10^43','a*x^b+c*x^d+e*x^f+g')
% y=f(time_rsg(2:50)');
% h=plot(log10(time_rsg(2:50))',log10(y)','*r'),hold on,
f = fit(log10(time_rsg(23:50)),log10(2*luminosityRSG.luminosity(23:50))','gauss6')
 h=plot(f,log10(time_rsg(23:50))',log10(2*luminosityRSG.luminosity(23:50))','*r'),hold on,
 set(h,'MarkerSize',3.5,'color','r');
 
 %luminosityBSG.luminosity(22)=10^42.5;
 f = fit(log10(time_bsg(23:41)),log10(2*luminosityBSG.luminosity(23:41))','poly7')
 h1=plot(f,log10(time_bsg(23:41))',log10(2*luminosityBSG.luminosity(23:41))','.b'),hold on,
 set(h1,'MarkerSize',8,'color','b');
 
 f = fit(log10(time_ic(23:32)),log10(2*luminosityIc.luminosity(23:32))','power2')
 h2=plot(f,log10(time_ic(23:32))',log10(2*luminosityIc.luminosity(23:32))','+k'),hold on,
 set(h2,'MarkerSize',8,'color','k');
%loglog(time_bsg(1:50),4*luminosityBSG.luminosity(1:50),'b','Linewidth',1),hold on,
%loglog(time_ic(1:50),4*luminosityIc.luminosity,'g','Linewidth',1),hold on,
%axis([10^-1,10^6,10^39,10^48])
xlabel('Log (t [sec])'); 
ylabel('Log (L_{tot} [erg/s])');
%title('Time of Maximum Compression Rate');
legend('RSG', 'fit','BSG','fit','Ic','fit','Location','northwest');
set(gca,'LineWidth',1.5,'FontSize',10);
name=['/home/nilou/Data/totlumall_new.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
