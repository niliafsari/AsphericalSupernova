mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


pconv=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

theta=linspace(0,pi/2,256);
radius=linspace(0,2,256);

xx=zeros(256,256);
yy=zeros(256,256);

for i=1:256
    for j=1:256
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end

for t=23:23
    name=['/home/nilou/Data/processeddata/dens256_' int2str(t) '.csv'] ;
    density= csvread(name)/rhotot;
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=450;
    height=400;
    myFigure = figure('PaperPositionMode','auto','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    h=surf(xx/rtot,yy/rtot, log10(density));hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    c=colorbar('Location','eastoutside');
  
    ylabel(c,'log (\rho / \rho_*)','FontSize',12) 
    caxis auto
    view(2);
    
    %title(['Ic Model t=' int2str(t)]);
    xlabel('x/R_*'); 
    ylabel('y/R_*');
    %legend('\rho / \rho_*','Location','northeast');
    set(gca,'LineWidth',2,'FontSize',12);
    name=['density256_1.pdf'];
    set(gcf, 'PaperPosition', [0 0 5 4.5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 4.5]); %Set the paper to have width 5 and height 5.
    %print(gcf, '-dpng',name)
    %print(gcf, '-dpdf','-r10',name) 
    %export_fig(name, '-pdf','-r5')
end

