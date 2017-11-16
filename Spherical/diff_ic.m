mtot=0.0048187313*2;
etot=9.1997e-4*2;
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


r=0.5*rconv;

radius=linspace(0,8,16384)*rconv;
% diffic=zeros(143,1);
% time=load('/home/nilou/bgq/output/timesteps.mat');
% 
% for t=1:143
%     t
%     name=[ '/home/nilou/bgq/processeddata/Ic/dparamIc_' int2str(t-1) '.csv'];
%     diff_frontRSG= csvread(name);
%     [x,i]=min(abs(diff_frontRSG-1));
%     diffic(t)=radius(i);
% end
close all
f = fit(time.time1,diffic,'gauss6')
 h=plot(f,time.time1,diffic),hold on,

 diffic=f(time.time1);

save('/home/nilou/bgq/processeddata/Ic/diffic.mat','diffic');