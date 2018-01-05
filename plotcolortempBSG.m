close all
mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=15*msun;
e=1e51;
r=49*rsun;

kappa=0.34;
c=3e10;
h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;
Z=0.005;
a=7.566e-15;
X=0.7;

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
% t_color=load('/home/nilou/Data/processeddata/BSG/colortemp_tot_v3.mat','t_color_tot');
% t_color_bsg=t_color.t_color_bsg;
load('/home/nilou/Data/Spherical/processeddata/BSG/colortemp_BSG_spherical.mat','t_loc_all','time_rsg','etha_all');

load([path '/processeddata/BSG/colortem_1024_tot.mat'],'t_tot')
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height]) 


time=load('/home/nilou/Data/processeddata/timesteps_1024.mat');
time_bsg=(time.time1*tconv);



% t1=13*60*(m/(10*msun))^-0.24*(r/(20*rsun))^0.94*(e/1e51)^0.29;
% t2=20*60*(m/(10*msun))^-0.77*(r/(20*rsun))^0.62*(e/1e51)^0.99;
% 
% T_c(t_sph<t_s)=150*(m/(10*msun))^-1*(r/(20*rsun))^-0.1 *(e/1e51)^1.2* (t_sph(t_sph<t_s)/60).^-0.39;
% T_c(t_sph>t_s & t_sph<t1)=70*(m/(15*msun))^-0.9*(r/(20*rsun))^-0.7 *(e/1e51)^1.1* (t_sph(t_sph>t_s & t_sph<t1)/t_s).^-2.1;
% T_c(t_sph>t1 & t_sph<t2)=15*(m/(10*msun))^0.05*(r/(20*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/(15*60)).^-0.4;
% T_c(t_sph>t2)=7*(m/(10*msun))^-0.11*(r/(20*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;

ind=[173 194   201   202   203   204   205   206   208   209 210   213   214   215]
time_bsg(ind)=[];
t_tot(ind)=[];

f= fit(log10(time_bsg(172:352-14)),log10(t_tot(172:352-14)/11604.52)','poly5');
t_t=f(log10(time_bsg(172:352-14)'));
plot(log10(time_bsg(172:352-14)),t_t,'-k','LineWidth',1.5); hold on 

a=7.566e-15;
mtot=9.651994e-3;
etot=2.997521e-2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=15*msun;
e=1e51;
r=49*rsun;
kappa=0.34;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

t_sph=time_rsg(146:230)-time_rsg(145);
T_c=zeros(length(t_sph),1);
t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
T_c(t_sph<t_s)=50*(m/(25*msun))^-0.19*(r/(70*rsun))^0.06 *(e/1e51)^0.22* (t_sph(t_sph<t_s)/60).^-(16/45);
T_c(t_sph>t_s)=10*(m/(25*msun))^-0.11*(r/(70*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t_s)/3600).^-0.61;



% plot(log10(time_bsg(23:42)),log10(t_color_bsg(1,23:42)/11604.52),'k','LineWidth',2); hold on
% 
f= fit(log10(time_rsg(146:1:230)-4.3021e+03),log10(t_loc_all(146:230)/(11604.52)),'smoothingspline');
t_t=f(log10(time_rsg(146:1:230)-4.3021e+03));
T_c_BSG=(log10(T_c)+0.3*t_t)/1.3;
time_BSG_log=log10(time_rsg(146:1:230)-4.3021e+03);
T_c_BSG(18:27)=[];
time_BSG_log(18:27)=[];
plot(time_BSG_log,T_c_BSG,'--k','LineWidth',1.5),hold on,
plot(log10(time_rsg(146:1:230)-4.3021e+03),log10(T_c),'-.k','LineWidth',1.5);

save('/home/nilou/Data/Spherical/processeddata/BSG/T_c_BSG_spherical.mat','T_c_BSG','time_BSG_log','time_rsg');

xlabel('Log (t [sec])'); 
ylabel('Log (T_{color} [eV])');
legend('Total','Spherical','N&S 2010')
%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortempBSG_1024_v1.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')