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

h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;

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
load('/home/nilou/Data/processeddata/RSG/spherical_RSG.mat','log_l','time_axis_log','spherical_LC_NS');
load([path '/processeddata/RSG/luminosity_1024_0.mat'],'luminosity')
load([path '/processeddata/RSG/luminosity_1024_90.mat'],'luminosity90')
load([path '/processeddata/RSG/luminosity_1024_tot.mat'],'luminosity_tot')
time=load([path '/processeddata/timesteps_1024.mat']);
all=[0 inf];
X=[4.83e16 4.83e19 ];
FUV=c./([1806 1340]*10^-8);
NUV=c./([3007 1693]*10^-8);
U=c./([4028 3048]*10^-8);
G=c./([5549 3783]*10^-8);
R= c./([6989 5415]*10^-8);
I= c./([8389 6689]*10^-8);
Z= c./([10833 7960]*10^-8);

bands=[ G; FUV ; X];

load([path '/processeddata/RSG/colortem_1024_tot_v2.mat'],'t_tot')
load('/home/nilou/Data/Spherical/processeddata/RSG/colortemp_RSG_spherical.mat','t_loc_all','time_rsg','etha_all');
time_rsg1=time_rsg;
time_rsg=(time.time1*tconv);

f= fit(log10(time_rsg(173:1:420)),log10(luminosity(1,173:1:420))','rat55');
luminosity_fit=10.^f(log10(time_rsg(173:1:420)'));
%plot(log10(time_rsg(173:1:420)),luminosity_fit,'-r','LineWidth',1.5),hold on,

f= fit(log10(time_rsg(173:1:420)),log10(2*luminosity90(1,173:1:420))','rat55');
luminosity90_fit=10.^f(log10(time_rsg(173:1:420)'));
%luminosity90_fit(40:50)= 42.8542;
%plot(log10(time_rsg(173:1:420)),luminosity90_fit,'--k','LineWidth',1.5),hold on,

ff= fit(log10(time_rsg(173:1:420)),log10(2*luminosity_tot(1,173:1:420))','rat55');
luminosity_tot_fit=10.^ff(log10(time_rsg(173:1:420)'));
luminosity_tot_fit(1)= 10.^43;
%plot(log10(time_rsg(173:1:420)),luminosity_tot_fit,'-.b','LineWidth',1.5),hold on,

f=fit(log10(time_rsg1(141:270)-time_rsg1(109)+25000),log_l(1:130),'poly6');
L_s=10.^f(log10(time_rsg(173:1:420)));

f= fit(log10(time_rsg(173:1:420)),log10(t_tot(173:420)/(1.5*11604.52))','poly3');
t_color_rsg=(10.^f(log10(time_rsg(173:1:420)')))*11604.52;
%plot(log10(time_rsg(173:1:420)),t_t,'-k','LineWidth',1.5),hold on,

load('/home/nilou/Data/Spherical/processeddata/RSG/T_c_RSG_spherical.mat','T_c_RSG','time_RSG_log','time_rsg1');

f= fit(time_RSG_log,T_c_RSG,'smoothingspline');
T_c=(10.^f(log10(time_rsg(173:1:420)')))*11604.52;
%plot(log10(time_rsg1(109:210)-time_rsg1(109)+time_rsg1(65)+927.3),t_t,'--k','LineWidth',1.5),hold on,

% clear T_c
% t_sph=time_rsg1(141:270)-time_rsg1(140);
% t_s=14*3600* (m/(15*msun))^0.43*(r/(500*rsun))^1.26*(e/1e51)^-0.56;
% T_c(t_sph<t_s)=10*(m/(15*msun))^-0.22*(r/(500*rsun))^0.12 *(e/1e51)^0.23* (t_sph(t_sph<t_s)/3600).^-0.36;
% T_c(t_sph>t_s )=3*(m/(15*msun))^-0.13*(r/(500*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t_s)/(24*3600)).^-0.56;
% T_c=T_c*11604.52;
% 
% f= fit(log10(time_rsg1(141:270)-time_rsg1(109)+25000),log10(T_c'/(11604.52)),'smoothingspline');
% T_c=(10.^f(log10(time_rsg(173:1:420)')))*11604.52;

bandL=zeros(length(bands),248,4);
factor=zeros(1,248);
factor1=zeros(1,248);

for i=1:length(bands)
    if (i==length(bands))
        bandL(i,:,1)=luminosity_tot_fit;
        bandL(i,:,2)=luminosity_fit;
        bandL(i,:,3)=luminosity90_fit;
        bandL(i,:,4)=L_s;
    else 
        for j=1:248
            T=t_color_rsg(j);
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*T))-1);
            factor(j)= integral(fun,bands(i,1),bands(i,2));
            Tp=T_c(j);
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*Tp))-1);
            factor1(j)= integral(fun,bands(i,1),bands(i,2));
        end
        bandL(i,:,1)=((pi*luminosity_tot_fit')./(sigma*t_color_rsg'.^4)).*factor;
        bandL(i,:,2)=((pi*luminosity_fit')./(sigma*t_color_rsg'.^4)).*factor;
        bandL(i,:,3)=((pi*luminosity90_fit')./(sigma*t_color_rsg'.^4)).*factor;
        bandL(i,:,4)=((pi*L_s')./(sigma*T_c'.^4)).*factor1;
    end
end

set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
bandM=71.197425-2.5*log10(bandL*1e-7);
close all

x0=15;
y0=15;
width=350;
height=950;
myFigure = figure('PaperPositionMode','auto','Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

ax1=subplot(3,1,1)
plot(log10(time_rsg(173:1:420)),reshape(bandM(1,:,:),[248 4]),'LineWidth',1.5)

%s=plotyy(log10(time_bsg(173:2:352-14)),reshape(bandM(1,:,:),[83 4]),log10(time_bsg(173:2:352-14)),reshape(bandM(1,:,:),[83 4]))
%axis([3.7 4.15 33 45])
str={'M_{Tot}','M_{obs}(\Theta==0)', 'M_{obs}(\Theta==\pi/2)', 'Spherical'};


xlabel('Log(t [sec])'); 
text(5.04,-10,'g^\prime','HorizontalAlignment','left','FontSize',11);
ax1.Position=[0.210 0.6533 0.6000 0.2717]
% s(2).YDir='reverse';
% s(2).YTickLabel=num2cell(round(log10((10.^((71.197425-(-4:-1:-10))/2.5))/1e-7),1));
set(gca,'Ydir','reverse')
set(gca,'LineWidth',1.5,'FontSize',12);
ax1_1=subplot(3,1,2)
plot(log10(time_rsg(173:1:420)),reshape(bandM(2,:,:),[248 4]),'LineWidth',1.5)
%ylabel('Absolute Magnitude [mag]');
ax1_1.Position=[0.210 0.6533 0.6000 0.2717]
set(gca,'LineWidth',1.5,'FontSize',11);
set(gca,'Ydir','reverse')
text(5.04,-13.5,'FUV','HorizontalAlignment','left','FontSize',11);
% columnlegend(1,str,'location', 'northeast');
% legend('boxoff')
%axis([3.7 4.15 34 45])
set(gca,'LineWidth',1.5,'FontSize',12);
ax1_2=subplot(3,1,3)
plot(log10(time_rsg(173:1:420)),reshape(bandM(3,:,:),[248 4]),'LineWidth',1.5)
%axis([3.7 4.15 32 43])
ax1_2.Position=[0.210 0.6533 0.6000 0.2717]
text(5.04,-20,'X-ray','HorizontalAlignment','left','FontSize',11);
set(gca,'Ydir','reverse')
set(gca,'LineWidth',1.5,'FontSize',12);

% subplot(5,1,4)
% plot(log10(time_bsg(23:42)),reshape(bandM(4,:,:),[20 4]),'LineWidth',1.5)
%axis([3.7 4.1 -15 -24])
% text(3.8,2,'FUV','HorizontalAlignment','left','FontSize',11);
% set(gca,'Ydir','reverse')
% subplot(5,1,5)
% plot(log10(time_bsg(23:42)),reshape(bandM(5,:,:),[20 4]),'LineWidth',1.5)
% %axis([3.7 4.15 34 46])
% text(3.8,-20,'X-ray','HorizontalAlignment','left','FontSize',11);
% set(gca,'Ydir','reverse')samexaxis('abc','xmt','on','ytac','join','yld',1)

samexaxis('abc','xmt','on','ytac','join','yld',1)
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',... 
    'YAxisLocation','right',...
    'Color','none');
ax2.YDir='reverse';
ax2.XTick=[]
z=linspace(0,1,9);
ax2.YTick=z(2:2:8);
ax1.YTick
ax2.YTickLabel=num2cell(round(log10((10.^((71.197425-(-13:1:-10))/2.5))/1e-7),1));

set(gca,'LineWidth',1.5,'FontSize',12);

ax1_pos = ax1_1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',... 
    'YAxisLocation','right',...
    'Color','none');
ax2.YDir='reverse';
ax2.XTick=[]
z=linspace(0,1,9);
ax2.YTick=z(2:2:8);
ax2.YLabel.String='Log(L[erg/s])'
ax1_1.YTick
ax2.YTickLabel=num2cell(round(log10((10.^((71.197425-(-16:1:-13))/2.5))/1e-7),1));

 

% ax2.XTick=ax1.XTick/100
% r=1./sin(ax1.XTick*pi/180)
% ax2.XTickLabel=num2cell(round(r,2))
% ax2.XLabel.String='x/R_*'


set(gca,'LineWidth',1.5,'FontSize',12);

ax1_pos = ax1_2.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',... 
    'YAxisLocation','right',...
    'Color','none');
ax2.YDir='reverse';
ax2.XTick=[]
z=linspace(0,1,5);
ax2.YTick=z(1:4);
ax1_2.YTick=-24:2:-18;
ax2.YTickLabel=num2cell(round(log10((10.^((71.197425-(-24:2:-18))/2.5))/1e-7),1));

set(gca,'LineWidth',1.5,'FontSize',12);
h=legend(ax1,'Total','\Theta=0','\Theta=\pi/2','Spherical')
h.Position=[ 0.3298    0.82    0.2670    0.0950]
 name=['/home/nilou/Data/plot/BandM_RSG_1024.pdf'];    
 print('-dpdf',name) 
 export_fig(name, '-pdf')