mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;

m=5*msun;
e=1e51;
r=0.2*rsun;

X=0;
kappa_T=0.2*(1+X);
kappa=kappa_T;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv_ic=rconv/vconv;

all=[0 inf];
X=[4.83e16 2.41e18];
FUV=c./([1806 1340]*10^-8);
NUV=c./([3007 1693]*10^-8);
U=c./([4028 3048]*10^-8);
G=c./([5549 3783]*10^-8);
R= c./([6989 5415]*10^-8);
I= c./([8389 6689]*10^-8);
Z= c./([10833 7960]*10^-8);

bands=[ G;  FUV ; X ];



% lum=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_Phi0.mat','luminosity');
% luminosity_phi0=lum.luminosity(23:32)/10;
% lum=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_tot.mat','luminosity');
% luminosity_tot=lum.luminosity(23:32)/10;
% lum=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_Phi90.mat','luminosity');
% luminosity_phi90=lum.luminosity(23:32)/10;


load('/home/nilou/Data/processeddata/BSG/luminosity_0.mat','luminosity')
load('/home/nilou/Data/processeddata/BSG/luminosity_90.mat','luminosity90')
load('/home/nilou/Data/processeddata/BSG/luminosity_tot.mat','luminosity_tot')
load('/home/nilou/Data/processeddata/ic/spherical_ic.mat','log_l','time_axis_log','spherical_LC_NS');

luminosity(23)=4.2e42;
luminosity90(23)=2.1e42;
luminosity_tot(23)=2e42;



load('/home/nilou/Data/processeddata/ic/colortemp_ic_spherical.mat','t_loc_all','time_ic','etha_all');

fff= fit(log10(time_ic(105:140)-time_ic(105)+time_ic(66)),log10(t_loc_all(105:140)/11604.52),'rat34');



load('/home/nilou/Data/processeddata/ic/colortemp_tot_v1.mat','t_color_bsg');
t_color_ic=t_color_bsg(23:32)*5;

time=load('/home/nilou/Data/processeddata/timesteps.mat');
time_ic=(time.time1*tconv_ic);
T_c=(10.^fff(log10(time_ic(23:32)')))*11604.52;

f= fit(log10(time_ic(23:40)),log10(luminosity(23:40)/9.7)','rat44');
luminosity_fit=10.^f(log10(time_ic(23:32)'));
f= fit(log10(time_ic(23:40)),log10(2*luminosity90(23:40)/9.7)','rat44');
luminosity90_fit=10.^f(log10(time_ic(23:32)'));
f= fit(log10(time_ic(23:40)),log10(2*luminosity_tot(23:40)/9.7)','rat44');
luminosity_tot_fit=10.^f(log10(time_ic(23:32)'));

ff= fit(time_axis_log,log_l,'rat44');
L_s=10.^ff(log10(time_ic(23:32)'));



bandL=zeros(length(bands),10,4);
factor=zeros(1,10);

% t_sph=time_ic(23:32)-time_ic(23)+0.46;
% t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
% L_s=zeros(length(time_ic(23:32)),1);
% L_s(t_sph<t_s)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34*(t_sph(t_sph<t_s)./60).^(-4/3);
% L_s(t_sph>t_s)=2*10^42*(m/(15*msun))^-0.73*(r/(50*rsun))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(3600)).^-0.35;
% 
% T_c=zeros(10,1);
% t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
% t1=200*(m/(15*msun))^-0.24*(r/(5*rsun))^0.94*(e/1e51)^0.29;
% t2=400*(m/(15*msun))^-0.77*(r/(5*rsun))^0.62*(e/1e51)^0.99;
% 
% 
% T_c(t_sph<t_s)=10^3*(m/(15*msun))^-1.5*(r/(5*rsun))^-0.2 *(e/1e51)^1.4* t_sph(t_sph<t_s).^-0.4;
% T_c(t_sph>t_s & t_sph<t1)=140*(m/(15*msun))^-1.2*(r/(5*rsun))^-0.9 *(e/1e51)^1.7* (t_sph(t_sph>t_s & t_sph<t1)/t_s).^-2.2;
% T_c(t_sph>t1 & t_sph<t2)=40*(m/(15*msun))^0.05*(r/(5*rsun))^0.25 *(e/1e51)^-0.1* (t_sph(t_sph>t1 & t_sph<t2)/60).^-0.4;
% T_c(t_sph>t2)=5*(m/(15*msun))^-0.11*(r/(5*rsun))^0.38 *(e/1e51)^0.11* (t_sph(t_sph>t2)/3600).^-0.61;
% 
% 
% T_c=T_c*11604.52;

for i=1:length(bands)
    if (i==length(bands))
        bandL(i,:,1)=luminosity_tot_fit;
        bandL(i,:,2)=luminosity_fit;
        bandL(i,:,3)=luminosity90_fit;
        bandL(i,:,4)=L_s;
    else
        for j=1:10
            T=t_color_ic(j);
            Tp=T_c(j);
            
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*T))-1);
            factor(j)= integral(fun,bands(i,1),bands(i,2));
            fun = @(x) (2*h*x.^3/ c^2) .* 1./(exp(h*x/(k_B*Tp))-1);
            factor1(j)= integral(fun,bands(i,1),bands(i,2));
        end
        bandL(i,:,1)=((pi*luminosity_tot_fit')./(sigma*t_color_ic.^4)).*factor;
        bandL(i,:,2)=((pi*luminosity_fit')./(sigma*t_color_ic.^4)).*factor;
        bandL(i,:,3)=((pi*luminosity90_fit')./(sigma*t_color_ic.^4)).*factor;
        bandL(i,:,4)=((pi*L_s')./(sigma*T_c'.^4)).*factor1;

    end

end
%set(0,'DefaultAxesColorOrder',distinguishable_colors(20));
set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 0 0;1 0 1],'DefaultAxesLineStyleOrder','-|--|:|-.')

%set(0,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
bandM=71.197425-2.5*log10(bandL*1e-7);
close all
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=900;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

subplot(3,1,1)
plot(log10(time_ic(23:32)),reshape(bandM(1,:,:),[10 4]),'LineWidth',1.5)
% axis([1.07 1.35  26 42.5])
str={'M_{Tot}','M_{obs}(\Theta==0)', 'M_{obs}(\Theta==\pi/2)', 'Spherical'};

 text(1.15,-3,'g^\prime','HorizontalAlignment','left','FontSize',12);

set(gca,'Ydir','reverse')
set(gca,'LineWidth',1.5,'FontSize',12);
subplot(3,1,2)
plot(log10(time_ic(23:32)'),reshape(bandM(2,:,:),[10 4]),'LineWidth',1.5)
set(gca,'LineWidth',1.5,'FontSize',12);
% columnlegend(1,str,'location', 'northeast','FontSize',5);
% legend('boxoff')
ylabel('Absolute Magnitude [mag]');
set(gca,'Ydir','reverse')
 text(1.15,-5,'FUV','HorizontalAlignment','left','FontSize',12);
% axis([1.07 1.35  27 43])
subplot(3,1,3)
 plot(log10(time_ic(23:32)'),reshape(bandM(3,:,:),[10 4]),'LineWidth',1.5)
 text(1.15,-17,'X-ray','HorizontalAlignment','left','FontSize',12);
%axis([1.07 1.35  26 43])
set(gca,'Ydir','reverse')
% subplot(5,1,4)
%  plot(log10(time_ic(23:32)'),reshape(bandM(4,:,:),[10 4]),'LineWidth',1.5)
% %axis([1.07 1.35  27 43])
% text(1.15,-5,'FUV','HorizontalAlignment','left','FontSize',12);
% set(gca,'Ydir','reverse')
% subplot(5,1,5)
%  plot(log10(time_ic(23:32)'),reshape(bandM(5,:,:),[10 4]),'LineWidth',1.5)
% %axis([1.07 1.35  27 43])
% text(1.15,-5,'X-ray','HorizontalAlignment','left','FontSize',12);
% set(gca,'Ydir','reverse')
samexaxis('abc','xmt','on','ytac','join','yld',1)
xlabel('Log(t [sec])'); 
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/Data/plot/BandM_ic_v3.pdf'];    
%set(gcf, 'PaperSize', [2 10]); %Set the paper to have width 5 and height 5. 
print('-dpdf',name) 
 export_fig(name, '-pdf')
