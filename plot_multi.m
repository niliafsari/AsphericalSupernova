mtot=0.0048187313*2;
etot=0.0130949665*2;
rtot=0.5;


ptot=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048);

xx=zeros(2048,2048);
yy=zeros(2048,2048);

for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
time=load('/home/nilou/Data/processeddata/timesteps_1024.mat');
for t=200:200
%     name=['/home/nilou/Data/rawdata/pressure/pres1024_' int2str(t) '.csv'] ;
%     pres= csvread(name)/ptot;
%     name=['/home/nilou/Data/rawdata/density/dens1024_' int2str(t) '.csv'] ;
%     dens= csvread(name)/rhotot;
%     entropy=pres./(dens.^(4/3));
%     name=['/home/nilou/Data/rawdata/velocity/velr1024_' int2str(t) '.csv'] ;
%     velr= csvread(name)/vtot;
%     name=['/home/nilou/Data/rawdata/velocity/velth1024_' int2str(t) '.csv'] ;
%     velt= csvread(name)/vtot;
%     vel=sqrt(velr.^2.+velt.^2);
    close all
    a=get(gcf,'Position');
    x0=10;
    y0=10;
    width=450;
    height=550;
    myFigure = figure('PaperPositionMode','auto','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    
    s2=subplot(2,2,1)
    h(3)=surf(-xx/rtot,yy/rtot, log10(pres));
    grid off
    set(h(3),'LineStyle','none');
    %colormap jet
    c3=colorbar('Location','northoutside'); 
    c3.Position=[ax1_pos(1) ax1_pos(2)+ax1_pos(4)-0.05 c3.Position(3)/2 c3.Position(4)]
    c3.LimitsMode='manual';
    c3.Limits=[min(min(log10(pres))) max(max(log10(pres)))];
    c3.Color
    str = 'P/P_*';
    ylabel(c3,str,'FontSize',12) 
    view(2);    %**********
    
    s1=subplot(2,2,4)
    h(1)=surf(xx/rtot,-yy/rtot, log10(entropy));hold on
    grid off
    set(h(1),'LineStyle','none');
    axis( [-4 4 -4 4])
    axis equal
    axis( [-4 4 -4 4])
%     ax1 = s1; % current axes
%     ax1_pos = ax1.Position
    c1=colorbar('Location','southoutside');
    c1.Position=[ax1_pos(1)+ax1_pos(3)/2 ax1_pos(2)-c1.Position(4)+0.05 c1.Position(3)/2 c1.Position(4)]
    c1.LimitsMode='manual';
    c1.Limits=[min(min(log10(entropy),[],'omitnan')) max(max(log10(entropy),[],'omitnan'))];
    str = 'S/S*';
    ylabel(c1,str,'FontSize',12) 
    colormap jet
    %Cdata1=h(1).CData;
    view(2);
    set(gca,'LineWidth',2,'FontSize',12);
    subplot(2,2,2)
    h(2)=surf(xx/rtot,yy/rtot, log10(dens));
    grid off
    set(h(2),'LineStyle','none');
    colormap jet
    c2=colorbar('Location','northoutside'); 
    c2.Position=[ax1_pos(1)+ax1_pos(3)/2 ax1_pos(2)+ax1_pos(4)-0.05 c2.Position(3)/2 c2.Position(4)]
    c2.LimitsMode='manual';
    c2.Limits=[min(min(log10(dens))) max(max(log10(dens)))];
    c2.Color
    str = 'Log( \rho/ \rho_* )';
    ylabel(c2,str,'FontSize',12) 
    view(2);

    
    s(4)=subplot(2,2,3)
    h(4)=surf(-xx/rtot,-yy/rtot, (vel));
    grid off
    set(h(4),'LineStyle','none');
    %colormap jet
    c4=colorbar('Location','southoutside'); 
    c4.Position=[ax1_pos(1) ax1_pos(2)-c4.Position(4)+0.05 c4.Position(3)/2 c4.Position(4)]
    c4.LimitsMode='manual';
    c4.Limits=[median(median(vel,'omitnan'))-2 max(max(vel,[],'omitnan'))];
    c4.Color
    str = 'v/v_*';
    ylabel(c4,str,'FontSize',12) 
    view(2);
 
%     h(1).CDataMode='manual';
%     h(1).CData=Cdata1;
%     h(1).CDataMapping = 'direct' ;
    
    %xlabel('\theta_{tang}'); 
    ylabel('r/R_*');
%     ax1 = gca; % current axes
%     ax1_pos = ax1.Position; % position of first axes
%     ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
%     x2=linspace(0,90,10);
%     y2=linspace(-3,3,10);
%     line(x2,y2,'Parent',ax2,'Color','none')
%     set(ax2,'LineWidth',2,'FontSize',12);
%     ax2.XTick = 0:20:80 ;
%     ax2.XTickLabel =strsplit(num2str((1./cos((0:20:80)*pi/180)),2));
%     ax2.YColor = 'none';
%     xlabel(ax2,'r_{eq}/R_*');
    
    
    
    
    
    
    
    
%     set(gca,'LineWidth',2,'FontSize',12);
%     (time.time1(t+1)/ttot)
%     name=['/home/nilou/Data/plot/test_' num2str(t) '_' num2str(time.time1(t+1)/ttot) '.png'];
%     set(gcf, 'PaperPosition', [0 0 5 4.5]); %Position plot at left hand corner with width 5 and height 5.
%     set(gcf, 'PaperSize', [5 4.5]); %Set the paper to have width 5 and height 5.
%     print(gcf, '-dpng', '-r900', name)
    %print('-dpdf',name) 
    %export_fig(name, '-pdf','r10')
end

