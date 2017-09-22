mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=15*msun;
e=1e51;
r=49*rsun;
kappa=0.34;

h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;


X=0.7;
kappa_T=0.2*(1+X);
kappa=kappa_T;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv_bsg=rconv/vconv;

all=[0 inf];
X=[4.83e16 4.83e18];
FUV=c./([1806 1340]*10^-8);
NUV=c./([3007 1693]*10^-8);
U=c./([4028 3048]*10^-8);
G=c./([5549 3783]*10^-8);
R= c./([6989 5415]*10^-8);
I= c./([8389 6689]*10^-8);
Z= c./([10833 7960]*10^-8);

bands=[ G ; FUV; X];

cita=0;
if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end

load([path '/processeddata/BSG/luminosity_1024_0.mat'],'luminosity')
load([path '/processeddata/BSG/luminosity_1024_90.mat'],'luminosity90')
load([path '/processeddata/BSG/luminosity_1024_tot.mat'],'luminosity_tot')
load('/home/nilou/Data/processeddata/BSG/spherical_BSG.mat','log_l','time_axis_log','spherical_LC_NS');

time=load([path '/processeddata/timesteps_1024.mat']);

time_bsg=(time.time1*tconv_bsg);

ind=find(luminosity_tot>10^42.3 )
time.time1(ind(12:23))=[];
luminosity90(ind(12:23))=[];
luminosity_tot(ind(12:23))=[];
luminosity(ind(12:23))=[];

f= fit(log10(time_bsg(173:2:352)),log10(luminosity(1,173:2:352))','rat44');
luminosity_fit=10.^f(log10(time_bsg(173:2:352-14)'));
%plot(log10(time_bsg(173:2:352)),luminosity_fit,'-r','LineWidth',1.5),hold on,

fff= fit(log10(time_bsg(173:2:352)),log10(2*luminosity90(1,173:2:352))','rat44');
luminosity90_fit=10.^fff(log10(time_bsg(173:2:352-14)'));
%plot(log10(time_bsg(173:2:352)),luminosity90_fit,'--k','LineWidth',1.5),hold on,

ff= fit(time_axis_log(1:80),log_l(1:80),'poly4');
L_s=10.^ff(log10(time_bsg(173:2:352-14)'));
%plot(log10(time_bsg(173:2:352)),luminosity_tot_fit,'-.b','LineWidth',1.5),hold on,

f= fit(log10(time_bsg(173:2:352)),log10(2*luminosity_tot(1,173:2:352))','rat44');
luminosity_tot_fit=10.^f(log10(time_bsg(173:2:352-14)'));

load('/home/nilou/Data/processeddata/BSG/colortemp_tot.mat','t_color_bsg');
t_color_bsg=t_color_bsg(23:40);



% load('/home/nilou/Data/processeddata/BSG/colortemp_BSG_spherical.mat','t_loc_all','time_rsg','etha_all');
% 
load([path '/processeddata/BSG/colortem_1024_tot.mat'],'t_tot')
% 
ind=[173 194   201   202   203   204   205   206   208   209 210   213   214   215]
time_bsg(ind)=[];
t_tot(ind)=[];

t_sph=time_bsg(173:2:352-14)-time_bsg(173)+113.59;
T_c=zeros(83,1);
t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
% t1=13*60*(m/(10*msun))^-0.24*(r/(20*rsun))^0.94*(e/1e51)^0.29;
% t2=20*60*(m/(10*msun))^-0.77*(r/(20*rsun))^0.62*(e/1e51)^0.99;

T_c(t_sph<t_s)=50*(m/(25*msun))^-0.19*(r/(70*rsun))^0.06 *(e/1e51)^0.22* (t_sph(t_sph<t_s)/60).^-(16/45);
T_c(t_sph>t_s )=10*(m/(25*msun))^-0.11*(r/(70*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t_s)/3600).^-0.61;
% T_c(t_sph>t1 & t_sph<t2)=15*(m/(10*msun))^0.05*(r/(20*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/(15*60)).^-0.4;
% T_c(t_sph>t2)=7*(m/(10*msun))^-0.11*(r/(20*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;

t_c=T_c*11604.52;
f= fit(log10(time_bsg(173:2:352-14)),log10(t_c'/11604.52)','poly5');
T_c=(10.^f(log10(time_bsg(173:2:352-14)')))*11604.52;

fx= fit(log10(time_bsg(173:2:352-14)),log10(t_tot(173:2:352-14)/11604.52)','poly5');
t_color_bsg=(10.^fx(log10(time_bsg(173:2:352-14)')))*11604.52;
%plot(log10(time_bsg(172:352-14)),t_t,'-k','LineWidth',1.5); hold on 

%plot(log10(time_rsg(114:180)-time_rsg(114)+time_rsg(67)),t_t,'--k','LineWidth',1.5),hold on,


bandL=zeros(length(bands),83,4);
factor=zeros(1,83);
factor1=zeros(1,83);



% t_sph=time_bsg(23:40)-time_bsg(23)+113.59;
% T_c=zeros(18,1);
% t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
% t1=13*60*(m/(10*msun))^-0.24*(r/(20*rsun))^0.94*(e/1e51)^0.29;
% t2=20*60*(m/(10*msun))^-0.77*(r/(20*rsun))^0.62*(e/1e51)^0.99;
% 
% T_c(t_sph<t_s)=150*(m/(10*msun))^-1*(r/(20*rsun))^-0.1 *(e/1e51)^1.2* (t_sph(t_sph<t_s)/60).^-0.39;
% T_c(t_sph>t_s & t_sph<t1)=70*(m/(15*msun))^-0.9*(r/(20*rsun))^-0.7 *(e/1e51)^1.1* (t_sph(t_sph>t_s & t_sph<t1)/t_s).^-2.1;
% T_c(t_sph>t1 & t_sph<t2)=15*(m/(10*msun))^0.05*(r/(20*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/(15*60)).^-0.4;
% T_c(t_sph>t2)=7*(m/(10*msun))^-0.11*(r/(20*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;
% 
% L_s=zeros(length(time_bsg(23:40)),1);
% L_s(t_sph<t_s)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34*(t_sph(t_sph<t_s)./60).^(-4/3);
% L_s(t_sph>t_s)=2*10^42*(m/(15*msun))^-0.73*(r/(50*rsun))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(3600)).^-0.35;
% 
% T_c=T_c*11604.52;

for i=1:length(bands)
    if (i==length(bands))
        bandL(i,:,1)=luminosity_tot_fit;
        bandL(i,:,2)=luminosity_fit;
        bandL(i,:,3)=luminosity90_fit;
        bandL(i,:,4)=L_s;
    else 
        for j=1:83
            T=t_color_bsg(j);
            Tp=T_c(j);
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*T))-1);
            factor(j)= integral(fun,bands(i,1),bands(i,2));
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*Tp))-1);
            factor1(j)= integral(fun,bands(i,1),bands(i,2));
        end
        bandL(i,:,1)=((pi*luminosity_tot_fit')./(sigma*t_color_bsg'.^4)).*factor;
        bandL(i,:,2)=((pi*luminosity_fit')./(sigma*t_color_bsg'.^4)).*factor;
        bandL(i,:,3)=((pi*luminosity90_fit')./(sigma*t_color_bsg'.^4)).*factor;
        bandL(i,:,4)=((pi*L_s)./(sigma*T_c.^4)).*factor1';

    end

end
%set(0,'DefaultAxesColorOrder',distinguishable_colors(20));
%set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 0 0;1 0 1],'DefaultAxesLineStyleOrder','-|--|:|-.')

set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
bandM=71.197425-2.5*log10(bandL*1e-7);
close all

x0=15;
y0=15;
width=350;
height=900;
myFigure = figure('PaperPositionMode','auto','Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

ax1=subplot(3,1,1)
plot(log10(time_bsg(173:2:352-14)),reshape(bandM(1,:,:),[83 4]),'LineWidth',1.5)

%s=plotyy(log10(time_bsg(173:2:352-14)),reshape(bandM(1,:,:),[83 4]),log10(time_bsg(173:2:352-14)),reshape(bandM(1,:,:),[83 4]))
%axis([3.7 4.15 33 45])
str={'M_{Tot}','M_{obs}(\Theta==0)', 'M_{obs}(\Theta==\pi/2)', 'Spherical'};


xlabel('Log(t [sec])'); 
text(3.8,-5,'g^\prime','HorizontalAlignment','left','FontSize',11);
ax1.Position=[0.210 0.6533 0.6000 0.2717]
% s(2).YDir='reverse';
% s(2).YTickLabel=num2cell(round(log10((10.^((71.197425-(-4:-1:-10))/2.5))/1e-7),1));
set(gca,'Ydir','reverse')
set(gca,'LineWidth',1.5,'FontSize',12);
ax1_1=subplot(3,1,2)
plot(log10(time_bsg(173:2:352-14)),reshape(bandM(2,:,:),[83 4]),'LineWidth',1.5)
%ylabel('Absolute Magnitude [mag]');
ax1_1.Position=[0.210 0.6533 0.6000 0.2717]
set(gca,'LineWidth',1.5,'FontSize',11);
set(gca,'Ydir','reverse')
text(3.8,-8,'FUV','HorizontalAlignment','left','FontSize',11);
% columnlegend(1,str,'location', 'northeast');
% legend('boxoff')
%axis([3.7 4.15 34 45])
set(gca,'LineWidth',1.5,'FontSize',12);
ax1_2=subplot(3,1,3)
plot(log10(time_bsg(173:2:352-14)),reshape(bandM(3,:,:),[83 4]),'LineWidth',1.5)
%axis([3.7 4.15 32 43])
ax1_2.Position=[0.210 0.6533 0.6000 0.2717]
text(3.8,-20,'X-ray','HorizontalAlignment','left','FontSize',11);
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
z=linspace(0,1,7);
ax2.YTick=z(1:6);
ax1.YTick=-10:-5;
ax2.YTickLabel=num2cell(round(log10((10.^((71.197425-(-10:-4))/2.5))/1e-7),1));

set(gca,'LineWidth',1.5,'FontSize',12);

ax1_pos = ax1_1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',... 
    'YAxisLocation','right',...
    'Color','none');
ax2.YDir='reverse';
ax2.XTick=[]
z=linspace(0,1,8);
ax2.YTick=z(1:7);
%ax2.YLabel.String='Log(L[erg/s])'
%ax1.YTick=-13:-8;
ax2.YTickLabel=num2cell(round(log10((10.^((71.197425-(-14:-7))/2.5))/1e-7),1));

 

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
ax1_2.YTick=-24:2:-16;
ax2.YTickLabel=num2cell(round(log10((10.^((71.197425-(-22:2:-14))/2.5))/1e-7),1));

set(gca,'LineWidth',1.5,'FontSize',12);
h=legend(ax1,'Total','\Theta=0','\Theta=\pi/2','Spherical')
h.Position=[ 0.5298    0.78    0.2670    0.0950]
 name=['/home/nilou/Data/plot/BandM_BSG_1024.pdf'];    
 print('-dpdf',name) 
 export_fig(name, '-pdf')
