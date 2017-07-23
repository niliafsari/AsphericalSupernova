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
cone=[0 0 0 0 0 0 0.02 0.02 0.1 0.1 0.1 0.15 0.17 0.17 0.2 0.2 0.22 0.23 0.23 0.24 0.25 0.27 0.3 0.32 ...
    0.33 0.34 0.34 0.34 0.34 0.34 0.34 0.35 0.36 0.36 0.36 0.37 0.37 0.38 0.39  0.4 0.4 0.41 0.42 0.42 0.42];

for t=33:71
    t
    time_current=time.time1(t)/ttot;
    name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t-1) '.csv'] ;
    velr= csvread(name)/vtot;
    name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t-1) '.csv'] ;
    velt= csvread(name)/vtot;
    v_y=-velt.*sin(thet) + velr.*cos(thet);
    v_x=velt.*cos(thet) + velr.*sin(thet);
    t_y=(yy(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3))-((xx(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3))-1)*cone(t-26)/3))./-v_y(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3));
    r_x=1+v_x(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3)).*t_y;
    t_y(isinf(t_y))=median(t_y)+10^5;
    angleth=acos(1./r_x)*180/pi; 
    
    xx_c=xx(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3));
    yy_c=yy(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3));
    
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
    dens= csvread(name)/rhotot;
    dv=(rad.^2).*sin(thet)*2*pi*(4/2048)*(pi/(2*2048));
    e_z=0.5*(dens(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3)).*((v_y(v_y<0 & rad>1& yy>((xx-1)*cone(t-26)/3)).^2 )));    
    E_z=0.5*(dens(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3)).*dv(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3))).*((v_y(v_y<0 & rad>1& yy>((xx-1)*cone(t-26)/3)).^2 ));    
    angleth_next=zeros(0,1);
    if (t<71)
        dt=(time.time1(t+1)/ttot-time_current);
        rnext=sqrt((xx+v_x*dt).^2 + (yy+v_y*dt).^2);
        rnext_c=rnext(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3));

        t_colnext =t_y(t_y< dt | rnext_c>=4 )+time_current;
        r_xnext =r_x(t_y< dt | rnext_c>=4);
        e_znext =e_z(t_y< dt | rnext_c>=4);
        E_znext =E_z(t_y< dt | rnext_c>=4);        
        xx_next=xx_c(t_y< dt | rnext_c>=4);
        yy_next=yy_c(t_y< dt | rnext_c>=4);
        if (~isempty(r_xnext ))
            angleth_next =acos(1./r_xnext )*180/pi; 
        end
    else
        t_colnext =t_y+time_current;
        r_xnext =r_x;
        e_znext =e_z;
        E_znext =E_z;
        xx_next=xx_c;
        yy_next=yy_c;
        if (~isempty(r_xnext))
            angleth_next =acos(1./r_xnext )*180/pi; 
        end
    end
    size(E_znext )
    totalenergy=sum(E_z);
    totalcone=sum(E_znext);
    
%      if (~isempty(r_xnext))
%          save(['/home/nilou/Data/processeddata/Collision/collision_comp_v3_' num2str(t) '.mat'],...
%              'e_znext','t_colnext','angleth_next','r_xnext','xx_next','yy_next', 'totalenergy',...
%              'totalcone','xx_c','yy_c','e_z'); 
%      end
%     close all
%     a=get(gcf,'Position');
%     x0=15;
%     y0=15;
%     width=350;
%     height=600;
%     myFigure = figure('PaperPositionMode','auto','Position',a);    
%     set(myFigure,'units','points','position',[x0,y0,width,height]) 
%     
%     subplot(2,1,1)
%     h=scatter(angleth,log10(t_y+time_current),2,real(log10(e_z)),'filled'); hold on
%     axis([0 90 -0.5 2])
%     colormap gray
%     h=scatter(angleth_next,log10(t_colnext),2,real(log10(e_znext)),'filled'); 
%     colormap jet
%     caxis([-10 -7]) 
%     c=colorbar('Location','eastoutside');
%     axis([0 90 -0.5 2])
%     str = 'log(E_z/E_*)';
%     ylabel(c,str,'FontSize',12) 
%     view(2);
%     xlabel('\theta_{tang}'); 
%     ylabel('Log t_{col}/t_*');
%     set(gca,'LineWidth',2,'FontSize',12);
        
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=600;
    height=500;
    myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');    
    set(myFigure,'units','points','position',[x0,y0,width,height]) 

    ax1=axes;
    h3=scatter(ax1,xx(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3)),yy(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3)),3,...
        log10(e_z),'filled');
    colormap(ax1,'gray')
    axis(ax1,[0 4 0 4])
    ax2=axes;
    h4=scatter(ax2,xx_next,yy_next,3,...
    log10(e_znext),'filled');
    linkaxes([ax1,ax2])
    axis(ax2,[0 4 0 4])
    colormap(ax1,'gray')
    colormap(ax2,'jet')
    caxis(ax2,[-5 -2])
    caxis(ax1,[-5 -2])
    str = 'log(E_z/E_*)';
    set(ax1,'LineWidth',2,'FontSize',12);
    set(ax2,'LineWidth',2,'FontSize',12);
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    ax1.XTick = [0 1 2 3 4];
    ax1.YTick = [0 1 2 3 4];
    title(ax1,['E_{z,n;t}/E_*= ' num2str(totalcone,3) ', E_{z,tot;t}/E_*=' num2str(totalenergy,3) ...
        ', t/t_*=' num2str(time_current,3)]);    
    xlabel(ax1,'x/R_*'); 
    ylabel(ax1,'y/y_*');
    set([ax1,ax2],'Position',[.18 .11 .685 .815]);
    cb1 = colorbar(ax1,'Position',[.08 .11 .04 .815]);
    cb2 = colorbar(ax2,'Position',[.88 .11 .04 .815]);
    str = 'log[ (dE_{*,z}/dV_*) ]';
    ylabel(cb1,str,'FontSize',12)
    ylabel(cb2,str,'FontSize',12)
     name=['/home/nilou/Data/plot/collision/gridcol_v5_' num2str(t) '.png'];
    % set(gcf, 'PaperPosition', [0 0 3.5 6]); %Position plot at left hand corner with width 5 and height 5.
    % set(gcf, 'PaperSize', [3.5 6]); %Set the paper to have width 5 and height 5.
     print(gcf, '-dpng', '-r100', name)
     export_fig(name, '-dpng', '-r100')  
     clear e_znext t_colnext angleth_next r_xnext
end


% grid off
% axis([0 90 -0.5 2])
% %caxis([log10(median(e_z)/((pi/4096)*(4/2048)))-3,log10(max(e_z)/((pi/4096)*(4/2048)))])
% caxis([-10 -7])
% colormap jet
% c=colorbar('Location','eastoutside');
% name=['/home/nilou/Data/plot/collision/col_comp_energy.png'];
% set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
% print(gcf, '-dpng', '-r900', name)
