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
load('/home/nilou/Data/processeddata/BSG/colortemp_BSG_spherical.mat','t_loc_all','time_rsg','etha_all');

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

% t_sph=time_bsg(23:41)-time_bsg(23)+113.59;
% T_c=zeros(41-23+2,1);
% t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
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

f= fit(log10(time_rsg(114:180)-time_rsg(114)+time_rsg(67)),log10(t_loc_all(114:180)/11604.52),'poly5');
t_t=f(log10(time_rsg(114:180)-time_rsg(114)+time_rsg(67)));
plot(log10(time_rsg(114:180)-time_rsg(114)+time_rsg(67)),t_t,'--k','LineWidth',1.5),hold on,


% plot(log10(time_bsg(23:42)),log10(t_color_bsg(1,23:42)/11604.52),'k','LineWidth',2); hold on
% plot(log10(time_bsg(23:42)),log10(T_c),'--k','LineWidth',2);

xlabel('Log (t [sec])'); 
ylabel('Log (T_{color} [eV])');
legend('Total','Spherical')
%title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/colortempBSG_1024.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')