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
a=7.566e-15;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

r=0.5*rconv;

radius=linspace(0,8,16384)*rconv;
diffrsg=zeros(143,1);
time=load('/home/nilou/bgq/output/timesteps.mat')
time_ic=time.time1*tconv;
diff=load('/home/nilou/bgq/processeddata/Ic/diffic.mat');
diff_front=diff.diffic;
ltot=zeros(143,1);
t_loc_all=zeros(267,1);
t_loc_nonth=zeros(267,1);
etha_all=zeros(267,1);
for t=1:143
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
    t_loc_all(t)=(3*pres(i)/a).^0.25;
    tdiff_all(t)=(3*kappa*dens(i)).*(pres(i).^2)./(c*gradp(i).^2);
    dens_s(t)=dens(i);
    pres_s(t)=pres(i);
    gradp_s(t)=gradp(i);
    etha_all(t)=(7e5./(tdiff_all(t))).*((dens(i)/1e-10).^-2).*(t_loc_all(t)/(100*11604.52)).^3.5;
    etha_all(t)=max(1,etha_all(t));
    t_loc_nonth(t)=t_loc_all(t)*etha_all(t)^2;
end
%save('/home/nilou/Data/processeddata/ic/colortemp_ic_spherical.mat','t_loc_all','time_ic','etha_all');

% t_sph=time_ic(116:143)-time_ic(116)+0.3;
% t_c=0.46;

% t_s=0.5*3600* (m/(15*msun))^0.41*(r/(50*rsun))^1.33*(e/1e51)^-0.58;
% L_s=zeros(length(time_ic(23:32)),1);
% %L_s(t_sph<t_c)=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34;
% L_s(t_sph<t_s )=2.5*1e44*(m/(15*msun))^-0.33*(r/(50*rsun))^2.3*(e/1e51)^0.34*(t_sph(t_sph<t_s )./60).^(-4/3);
% L_s(t_sph>t_s)=2*10^42*(m/(15*msun))^-0.73*(r/(50*rsun))*(e/1e51)^0.91.*(t_sph(t_sph>t_s)./(3600)).^-0.35;
% 
% a=get(gcf,'Position');
% x0=15;
% y0=15;
% width=350;
% height=300;
% close all
% 
% myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
% set(myFigure,'units','points','position',[x0,y0,width,height])
% % 
% %  f = fit(log10(time_ic(116:143)-time_ic(116)+time_ic(66)),log10(ltot(116:143)),'rat55')
% % 
% %  %h=plot(f,log10(time_ic(116:143)),log10(ltot(116:143))),hold on,
% % % 
% t_c=fittedmodel1(time_log_spherical);
% % % 
%  plot(log10(time_ic(116:143)-time_ic(116)+time_ic(66)),log_l,'LineWidth',2), hold on
% plot(log10(time_ic(116:143)-time_ic(116)+time_ic(66)),log10(L_s),'-.','LineWidth',2)
% 
% time_axis_log=log10(time_ic(116:143)-time_ic(116)+time_ic(66));
% spherical_LC_NS=log10(L_s);
% save('/home/nilou/Data/processeddata/ic/spherical_ic.mat','log_l','time_axis_log','spherical_LC_NS');
% 
% 
% set(gca,'LineWidth',1.5,'FontSize',12);
% xlabel('Log(t [sec])'); 
% ylabel('Log (L_{tot} [erg/s])');
% legend('Simulation','N&S 2010','Location','northeast');
% set(gca,'LineWidth',1.5,'FontSize',12);
% name=['/home/nilou/bgq/plot/spherical_ic.pdf'];
% print('-dpdf',name) 
% export_fig(name, '-pdf')
