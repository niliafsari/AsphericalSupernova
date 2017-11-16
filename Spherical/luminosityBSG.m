mtot=0.0048187313*2;
etot=9.1997e-4*2;
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

r=0.5*rconv;

radius=linspace(0,8,16384)*rconv;
diffrsg=zeros(200,1);
time=load('/home/nilou/bgq/output/timesteps.mat')
time_rsg=time.time1*tconv;
diff=load('/home/nilou/bgq/processeddata/BSG/diffrsg.mat');
diff_front=diff.diffbsg;
ltot=zeros(200,1);
for t=1:200
    name=['/home/nilou/bgq/raw_data/pres16384_' int2str(t-1) '.csv'] ;
    pres= csvread(name).*pconv;
    name=['/home/nilou/bgq/raw_data/dens16384_' int2str(t-1) '.csv'] ;
    dens= csvread(name).*rhoconv;
    name=['/home/nilou/bgq/raw_data/velr16384_' int2str(t-1) '.csv'] ;
    velr= csvread(name).*vconv;
    gradp=zeros(16384,1);
    dr=8*rconv/16384;
    for i=1:16384 
        if (i==16384)
            delvel=(velr(i)-velr(i-1))/dr ;
            gradp(i)=(pres(i)-pres(i-1))/dr;
        else
            delvel=(velr(i+1)-velr(i))/dr ;
            gradp(i)=(pres(i+1)-pres(i))/dr;
        end
    end
    [u,i]=min(abs(radius-diff_front(t)));
    ltot(t)=-4*pi*diff_front(t)^2*c*gradp(i)/(kappa*dens(i));
end

m=15*msun;
e=1e51;
r=49*rsun;

t_sph=time_rsg(114:200)-time_rsg(114);
t_c=0;
L_s=zeros(length(t_sph),1);
t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
L_s(t_sph<t_c)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34;
L_s(t_sph<t_s & t_sph>t_c)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34*(t_sph(t_sph<t_s & t_sph>t_c)./60).^(-4/3);
L_s(t_sph>t_s)=2*10^42*(m/(15*msun))^-0.73*(r/(50*rsun))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(3600)).^-0.35;
%L_s(t_sph>t_s)=3.3*10^42*(m/(msun))^-0.74*(r/(1e12))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(24*3600)).^-0.34;


a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
close all

myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

f = fit(log10(time_rsg(114:200)-time_rsg(114)+time_rsg(67)),log10(ltot(114:200)),'gauss5')
% h=plot(f,log10(time_rsg(109:143)-time_rsg(109)+time_rsg(53)),log10(ltot(109:143))),hold on,

log_l=f(log10(time_rsg(114:200))-time_rsg(114)+time_rsg(67));


plot(log10(time_rsg(115:200)-time_rsg(114)+time_rsg(67)),log10(ltot(115:200)),'LineWidth',2), hold on
plot(log10(time_rsg(114:200)-time_rsg(114)+time_rsg(67)),log10(L_s),'-.','LineWidth',2)
log_l=log10(ltot(115:200));
time_axis_log=log10(time_rsg(115:200)-time_rsg(114)+time_rsg(67));
spherical_LC_NS=log10(L_s);
save('/home/nilou/Data/processeddata/BSG/spherical_BSG.mat','log_l','time_axis_log','spherical_LC_NS');

set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('Log(t [sec])'); 
ylabel('Log (L_{tot} [erg/s])');
legend('Simulation','N&S 2010','Location','northeast');
set(gca,'LineWidth',1.5,'FontSize',12);
name=['/home/nilou/bgq/plot/spherical_bsg.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
