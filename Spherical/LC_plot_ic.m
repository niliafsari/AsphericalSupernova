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

%close all
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;

myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

% plot(log10(time_rsg(109:267)-time_rsg(109)+time_rsg(61)),log_l,'LineWidth',2), hold on
% plot(log10(time_rsg(109:267)-time_rsg(109)+time_rsg(61)),log10(L_s(1:267-108)),'-.','LineWidth',2)
% 
% time_axis_log=log10(time_rsg(109:267)-time_rsg(109)+time_rsg(61)+927.3);
% spherical_LC_NS=log10(L_s(1:267-108));
load('/home/nilou/Data/processeddata/ic/spherical_ic.mat','log_l','time_axis_log','spherical_LC_NS');
% find(luminosity>10^42.3)
% time.time1(luminosity>10^42.3)=[];
% luminosity90(luminosity>10^42.3)=[];
% luminosity_tot(luminosity>10^42.3)=[];
% luminosity(luminosity>10^42.3)=[];

% time_rsg=time.time1*tconv;
% f= fit(log10(time_rsg(161:7:203)),log10(luminosity(1,161:7:203))','poly5');
% luminosity_fit=f(log10(time_rsg(161:7:203)'));
% plot(log10(time_rsg(161:7:203)),luminosity_fit,'-r','LineWidth',1.5),hold on,
% 
% f= fit(log10(time_rsg(161:7:203)),log10(2*luminosity90(1,161:7:203))','poly5');
% luminosity90_fit=f(log10(time_rsg(161:7:203)'));
% %luminosity90_fit(40:50)= 42.8542;
% plot(log10(time_rsg(161:7:203)),luminosity90_fit,'--k','LineWidth',1.5),hold on,
% 
% f= fit(log10(time_rsg(161:7:203)),log10(2*luminosity_tot(1,161:7:203))','poly5');
% luminosity_tot_fit=f(log10(time_rsg(161:7:203)'));
% %luminosity_tot_fit(40:50)= 42.8742;
% plot(log10(time_rsg(161:7:203)),luminosity_tot_fit,'-.b','LineWidth',1.5),hold on,
% plot(time_axis_log(1:80),log_l(1:80),':m','LineWidth',1.5), hold on

plot(log10(time.time1(161:203)*tconv),log10(luminosity(1,161:203)),'-r','LineWidth',2),hold on
plot(log10(time.time1(161:203)*tconv),log10(2*luminosity90(1,161:203)),'-b','LineWidth',2),hold on
plot(log10(time.time1(161:203)*tconv),log10(2*luminosity_tot(1,161:203)),'-k','LineWidth',2),hold on

set(gca,'LineWidth',2,'FontSize',12);
xlabel('Log (t [sec])'); 
ylabel('Log (L_{tot} [erg/s])');
%title('Time of Maximum Compression Rate');
legend('\Theta=0','\Theta=\pi/2','Total','Spherical','Location','best');
name=['/home/nilou/Data/lum_1024_ic.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')

