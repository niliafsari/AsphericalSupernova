mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


ptot=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048)/rtot;

xx=zeros(2048,2048);
yy=zeros(2048,2048);

thet=zeros(2048,2048);


% for i=1:2048
%     for j=1:2048
%        xx(i,j)=radius(i)*sin(theta(j));
%        yy(i,j)=radius(i)*cos(theta(j)); 
%        thet(i,j)=theta(j);
%        rad(i,j)=radius(i);
%     end
% end
% time=load('/home/nilou/Data/timesteps.mat');
% figure
% e_cone_th=zeros(1,71-26);
% 
% cone=[0 0 0 0 0 0 0.02 0.02 0.1 0.1 0.1 0.15 0.17 0.17 0.2 0.2 0.22 0.23 0.23 0.24 0.25 0.27 0.3 0.32 ...
%     0.33 0.34 0.34 0.34 0.34 0.34 0.34 0.35 0.36 0.36 0.36 0.37 0.37 0.38 0.39  0.4 0.4 0.41 0.42 0.42 0.42];
% 
% for t=33:71
%     t
%     time_current=time.time1(t)/ttot;
%     name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t-1) '.csv'] ;
%     velr= csvread(name)/vtot;
%     name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t-1) '.csv'] ;
%     velt= csvread(name)/vtot;
%     v_y=-velt.*sin(thet) + velr.*cos(thet);
%     v_x=velt.*cos(thet) + velr.*sin(thet);
% %     t_y=(yy(v_y<0 &  rad>1 & yy<0.4)-0.4)./-v_y(v_y<0 & rad>1 & yy<0.4);
% %     r_x=1+v_x(v_y<0 & rad>1 & yy<0.4).*t_y;
% %     t_y(isinf(t_y))=median(t_y)+10^5;
% %     angleth=acos(1./r_x)*180/pi; 
%     
%     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
%     pres= csvread(name)/ptot;
%     dv=(rad.^2).*sin(thet)*2*pi*(4/2048)*(pi/(2*2048));
%     e_z=3*pres(v_y<0 & rad>1 & yy<((xx-1)*cone(t-26)/3) ).*dv(v_y<0 & rad>1 & yy<((xx-1)*cone(t-26)/3));    
%     
%     e_cone_th(t-26)=sum(e_z);
% end
% 
load('/home/nilou/Data/processeddata/Collision/all_cone.mat','e_cone_th'); 

plot(time.time1(27:71)/ttot,log10(e_cone_th/0.005));
set(gca,'LineWidth',1,'FontSize',10);

%set(h,'LineStyle','none');
%axis([0 90 -0.5 2])
xlabel(' t_{col}/t_*');
ylabel('Log dE_{th}/dt');
%caxis([log10(median(e_z)/((pi/4096)*(4/2048)))-3,log10(max(e_z)/((pi/4096)*(4/2048)))])
% caxis([-10 -7])
% colormap jet
% c=colorbar('Location','eastoutside');
% str = '$$ Log E_z $$';
% ylabel(c,str,'Interpreter','latex','FontSize',12) 
name='/home/nilou/Data/plot/collision/all_econevst_eth_v2.png';
set(gcf,'Color','w');
set(gcf, 'PaperPosition', [0 0 7 7]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 7]); %Set the paper to have width 5 and height 5.
print(gcf, '-dpng', '-r900', name)
export_fig(name, '-dpng', '-r900')
