mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


ptot=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048)/0.5;

xx=zeros(2048,2048);
yy=zeros(2048,2048);

for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
time=load('/home/nilou/Data/timesteps.mat');
for t=48:48
    name=['/home/nilou/Data/gfunction_' int2str(t) '.csv'] ;
    gf= csvread(name);
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=450;
    height=400;
    myFigure = figure('PaperPositionMode','auto','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    h=surf(xx,yy, log10(gf));hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    c=colorbar('Location','eastoutside');
    axis([0,4,0,4])
    ylabel(c,'log (g_\epsilon)','FontSize',12) 
    caxis auto
    view(2);
    
    %title(['Ic Model t=' int2str(t)]);
    xlabel('x/R_*'); 
    ylabel('y/R_*');
    %legend('\rho / \rho_*','Location','northeast');
    set(gca,'LineWidth',2,'FontSize',12);
    (time.time1(t+1)/ttot)
    name=['/home/nilou/Data/gf_' num2str(time.time1(t+1)/ttot) '.pdf'];
    set(gcf, 'PaperPosition', [0 0 5 4.5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 4.5]); %Set the paper to have width 5 and height 5.
    print(gcf, '-dpdf', '-r900', name)
    %print('-dpdf',name) 
    %export_fig(name, '-pdf','r10')
end

