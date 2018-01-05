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
path='/home/nilou/Data/';
for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
time=load('/home/nilou/Data/processeddata/timesteps_1024.mat');
for t=0
    if (t<=200 && t~=176)
        name=['/home/nilou/Data/rawdata/pressure/pres1024_' int2str(t) '.csv'] ;
        pres= csvread(name)/ptot;
        name=['/home/nilou/Data/rawdata/density/dens1024_' int2str(t) '.csv'] ;
        dens= csvread(name)/rhotot;
        name=['/home/nilou/Data/rawdata/velocity/velr1024_' int2str(t) '.csv'] ;
        velr= csvread(name)/vtot;
        name=['/home/nilou/Data/rawdata/velocity/velth1024_' int2str(t) '.csv'] ;
        velt= csvread(name)/vtot;   
    else
        name=[path '/rawdata/pressure/pres1024_' int2str(t) '.mat'] ;
        load(name,'pres_data');
        pres=pres_data/ptot;
        name=[path '/rawdata/density/dens1024_' int2str(t) '.mat'] ;
        load(name,'dens_data');
        dens=dens_data/rhotot;
        name=[path '/rawdata/velocity/velr1024_' int2str(t) '.mat'] ;
        load(name,'velx_data');
        velr=velx_data/vtot;
        name=[path '/rawdata/velocity/velth1024_' int2str(t) '.mat'] ;
        load(name,'vely_data');
        velt=vely_data/vtot;
    end
    entropy=pres./(dens.^(4/3));
    vel=sqrt(velr.^2+velt.^2);
    close all
    a=get(gcf,'Position');
    x0=10;
    y0=10;
    width=550;
    height=450;
    myFigure = figure('PaperPositionMode','auto','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    
    s1=subplot(2,2,1)
    h(1)=surf(-xx/rtot,yy/rtot, log10(pres));
    grid off
    set(h(1),'LineStyle','none');   
    view(2)
    c1=colorbar('Location','westoutside'); 
    axis equal
    set(gca,'LineWidth',2,'FontSize',12);
    axis([-4 0 0 4])
    s1.XTick=[]
    colormap jet
    s1.Position=[0.15 0.1+0.4 0.4 0.4];
    s1.XAxis.Color='none';
    str = 'Log( P / P_* )'; 
    ylabel(c1,str,'FontSize',12) 

    %*********
    
    s2=subplot(2,2,2)
    h(2)=surf(xx/rtot,yy/rtot, log10(dens));
    grid off
    set(h(2),'LineStyle','none');
    view(2)    
    c2=colorbar('Location','eastoutside'); 
    axis equal
    set(gca,'LineWidth',2,'FontSize',12);
    axis([0 4 0 4])
    s2.XTick=[]
    s2.YTick=[]
    colormap jet
    s2.Position=[0.15+0.325 0.1+0.4 0.4 0.4];
    s2.XAxis.Color='none';
    s2.YAxis.Color='none';
    str = 'Log( \rho/ \rho_* )'; 
    ylabel(c2,str,'FontSize',12) 
%*****
    s3=subplot(2,2,3);
    h(3)=surf(-xx/rtot,-yy/rtot, (vel));
    grid off
    set(h(3),'LineStyle','none');
    view(2)
    c3=colorbar('Location','westoutside');
    axis equal
    set(gca,'LineWidth',2,'FontSize',12);    
    axis([-4 0 -4 0])
    s3.XTick=[-4 -3 -2 -1 0]
    colormap jet
    s3.Position=[0.15 0.1 0.4 0.4];
    c3.Position=[0.109 0.1000 0.0364 0.4000]
    str = 'v/v_*';
    ylabel(c3,str,'FontSize',12) 
%**********
    s4=subplot(2,2,4)
    h(4)=surf(xx/rtot,-yy/rtot, log10(entropy));hold on
    grid off
    set(h(4),'LineStyle','none');
    axis equal
    view(2);
    c4=colorbar('Location','eastoutside'); 
    set(gca,'LineWidth',2,'FontSize',12);
    axis([0 4 -4 0])    
    colormap jet
    s4.YTick=[]
    s4.XTick=[0 1 2 3 4]
    s4.Position=[0.15+0.325 0.1 0.4 0.4];
    s4.YAxis.Color='none';
    c4.Position=[0.8570 0.1000 0.0364 0.4000]
    str = 'Log( S/ S_* )'; 
    ylabel(c4,str,'FontSize',12) 
  
   
%     (time.time1(t+1)/ttot)
      set(gcf, 'Color', 'w');
      name=['/home/nilou/Data/plot/allquant1024_' num2str(t) '_' num2str(time.time1(t+1)/ttot) '.png'];
      export_fig(name,'-m2.5')
      %      print('-dpdf',name) 
%      export_fig(name, '-pdf')
%     set(gcf, 'PaperPosition', [0 0 5 4.5]); %Position plot at left hand corner with width 5 and height 5.
%     set(gcf, 'PaperSize', [5 4.5]); %Set the paper to have width 5 and height 5.
%     print(gcf, '-dpng', '-r900', name)

end

