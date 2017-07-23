mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;

m=15*msun;
e=1e51;
r=49*rsun;

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

bands=[Z; I ; R; G; U; NUV; FUV ; X; all];



lum=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_Phi0.mat','luminosity');
luminosity_phi0=lum.luminosity(23:42);
lum=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_tot.mat','luminosity');
luminosity_tot=lum.luminosity(23:42)*2;
lum=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_Phi90.mat','luminosity');
luminosity_phi90=lum.luminosity(23:42);
load('/home/nilou/Data/processeddata/BSG/colortemp_tot.mat','t_color_bsg');
t_color_bsg=t_color_bsg(23:42);

time=load('/home/nilou/Data/timesteps.mat');
time_bsg=(time.time1*tconv_bsg);

bandL=zeros(length(bands),20,4);
factor=zeros(1,20);

t_sph=time_bsg(23:42)-time_bsg(23)+113.59;
T_c=zeros(41-23+2,1);
t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
t1=13*60*(m/(10*msun))^-0.24*(r/(20*rsun))^0.94*(e/1e51)^0.29;
t2=20*60*(m/(10*msun))^-0.77*(r/(20*rsun))^0.62*(e/1e51)^0.99;

T_c(t_sph<t_s)=150*(m/(10*msun))^-1*(r/(20*rsun))^-0.1 *(e/1e51)^1.2* (t_sph(t_sph<t_s)/60).^-0.39;
T_c(t_sph>t_s & t_sph<t1)=70*(m/(15*msun))^-0.9*(r/(20*rsun))^-0.7 *(e/1e51)^1.1* (t_sph(t_sph>t_s & t_sph<t1)/t_s).^-2.1;
T_c(t_sph>t1 & t_sph<t2)=15*(m/(10*msun))^0.05*(r/(20*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/(15*60)).^-0.4;
T_c(t_sph>t2)=7*(m/(10*msun))^-0.11*(r/(20*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;

L_s=zeros(length(time_bsg(23:41)),1);
L_s(t_sph<t_s)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34*(t_sph(t_sph<t_s)./60).^(-4/3);
L_s(t_sph>t_s)=2*10^42*(m/(15*msun))^-0.73*(r/(50*rsun))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(3600)).^-0.35;

T_c=T_c*11604.52;

for i=1:length(bands)
    if (i==length(bands))
        bandL(i,:,1)=luminosity_tot;
        bandL(i,:,2)=luminosity_phi0;
        bandL(i,:,3)=luminosity_phi90;
        bandL(i,:,4)=L_s;
    else 
        for j=1:20
            T=t_color_bsg(j);
            Tp=T_c(j);
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*T))-1);
            factor(j)= integral(fun,bands(i,1),bands(i,2));
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*Tp))-1);
            factor1(j)= integral(fun,bands(i,1),bands(i,2));
        end
        bandL(i,:,1)=((pi*luminosity_tot)./(sigma*t_color_bsg.^4)).*factor;
        bandL(i,:,2)=((pi*luminosity_phi0)./(sigma*t_color_bsg.^4)).*factor;
        bandL(i,:,3)=((pi*luminosity_phi90)./(sigma*t_color_bsg.^4)).*factor;
        bandL(i,:,4)=((pi*L_s')./(sigma*T_c'.^4)).*factor1;

    end

end
set(0,'DefaultAxesColorOrder',distinguishable_colors(20));
bandM=71.197425-2.5*log10(bandL*1e-7);
close all

x0=15;
y0=15;
width=350;
height=970;
myFigure = figure('PaperPositionMode','auto','Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

subplot(4,1,1)
plot(log10(time_bsg(23:42)'),bandM(:,:,1),'LineWidth',1.5)
%axis([3.7 4.15 33 45])
str={'z^\prime', 'i^\prime','r^\prime', 'g^\prime', 'u^\prime', ...
   'NUV', 'FUV' , 'X-ray', 'Bol'};


text(3.8,0,'Total Luminosity','HorizontalAlignment','left','FontSize',11);


set(gca,'Ydir','reverse')
set(gca,'LineWidth',1.5,'FontSize',12);
subplot(4,1,2)
plot(log10(time_bsg(23:42)'),bandM(:,:,2),'LineWidth',1.5)
ylabel('Absolute Magnitude [mag]');
set(gca,'LineWidth',1.5,'FontSize',11);
set(gca,'Ydir','reverse')
text(3.8,0,'Luminosity \Phi=0','HorizontalAlignment','left','FontSize',11);
columnlegend(3,str,'location', 'northeast');
legend('boxoff')
%axis([3.7 4.15 34 45])
subplot(4,1,3)
plot(log10(time_bsg(23:42)'),bandM(:,:,3),'LineWidth',1.5)
%axis([3.7 4.15 32 43])
text(3.8,5,'Luminosity \Phi=\pi/2','HorizontalAlignment','left','FontSize',11);
set(gca,'Ydir','reverse')
subplot(4,1,4)
plot(log10(time_bsg(23:42)'),bandM(:,:,4),'LineWidth',1.5)
%axis([3.7 4.15 34 46])
text(3.8,2,'Spherical Luminosity','HorizontalAlignment','left','FontSize',11);
set(gca,'Ydir','reverse')

samexaxis('abc','xmt','on','ytac','join','yld',1)
xlabel('Log(t [sec])'); 
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/BandL_BSG.pdf'];    
print('-dpdf',name) 
export_fig(name, '-pdf')
