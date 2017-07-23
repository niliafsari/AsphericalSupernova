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


for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
       thet(i,j)=theta(j);
       rad(i,j)=radius(i);
    end
end
time=load('/home/nilou/Data/timesteps.mat');
energyth=zeros(18,2048);
Vq=zeros(2048,2048);
for t=23:40
    load(['/home/nilou/Data/processeddata/Collision/collision_' num2str(t) '.mat'],'e_z','X','Y','Vq','t_y','angleth'); 
    Vq(isnan(Vq))=0;
    energyth(t-22,:)=sum(Vq,2);
end

    t_snap=time.time1(23:40)/ttot;
    angle_tang=linspace(0,90,2048);
    [angle_tang_m,t_snap_m] = meshgrid(angle_tang,t_snap);
    
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=350;
    height=600;
    myFigure = figure('PaperPositionMode','auto','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height])  

    totalenergy=sum(e_z);
    grid off
    h=surf(angle_tang_m,t_snap_m,log10(energyth));
    set(h,'LineStyle','none');
    %caxis([log10(median(e_z)/((6/2048)*(90/2048)))-3,log10(max(e_z)/((6/2048)*(90/2048)))])
    colormap jet
    c=colorbar('Location','eastoutside');
    %axis([0 90 -3 3])
    str = '$$ log( \frac{dE_z}{E_*}) $$';
    ylabel(c,str,'Interpreter','latex','FontSize',12) 
    view(2);
    xlabel('\theta_{tang}'); 
    ylabel('Log t_{snap}/t_*');
    set(gca,'LineWidth',2,'FontSize',12);

    name=['/home/nilou/Data/plot/collision/col_energy.png'];
    set(gcf, 'PaperPosition', [0 0 3.5 6]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [3.5 6]); %Set the paper to have width 5 and height 5.
    print(gcf, '-dpng', '-r900', name)