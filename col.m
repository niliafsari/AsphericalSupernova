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
for t=31:40
    name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t) '.csv'] ;
    velr= csvread(name)/vtot;
    name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t) '.csv'] ;
    velt= csvread(name)/vtot;
    v_y=-velt.*sin(thet) + velr.*cos(thet);
    v_x=velt.*cos(thet) + velr.*sin(thet);
    t_y=yy(v_y<0 &  rad>1)./-v_y(v_y<0 & rad>1);
    r_x=1+v_x(v_y<0 & rad>1).*t_y;
    t_y(isinf(t_y))=median(t_y)+10^5;
    angleth=acos(1./r_x)*180/pi;    
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t) '.csv'] ;
    dens= csvread(name)/rhotot;
    dv=(rad.^2).*sin(thet)*2*pi*(4/2048)*(pi/(2*2048));
    e_z=(dens(v_y<0 & rad>1).*dv(v_y<0 & rad>1)).*((v_y(v_y<0 & rad>1).^2 +v_x(v_y<0 & rad>1).^2));
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=350;
    height=600;
    myFigure = figure('PaperPositionMode','auto','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height])  

    anglethq=linspace(0,90,2048);
    t_yq=linspace(-3,3,2048);
    [X,Y] = meshgrid(anglethq,t_yq);
    Vq = griddata(real(angleth),real(log10(t_y)),e_z,X,Y, 'natural');
    for i=1:2047
        tt=min(log10(t_y( anglethq(1,i)<angleth & angleth<anglethq(1,i+1))));
        if (length(tt)>0)
            tt(tt<-3)=-3;
            tt(tt>3)=3;
            j=ceil(-(-3-tt(1))/(6/2048));
            if (j==0)
                j=1;
            end
            Vq(1:j,i)=10^-20;
        end
    end
    %Vq(X<real(min(angleth)))=10^-20;
    subplot(2,1,1)
    %h=scatter(real(angleth),real(log10(t_y)),3,log10(e_z),'filled');hold on
    totalenergy=sum(e_z);
    grid off
    h=surf(X,Y,real(log10(Vq/((6/2048)*(90/2048)))));
    set(h,'LineStyle','none');
    caxis([log10(median(e_z)/((6/2048)*(90/2048)))-3,log10(max(e_z)/((6/2048)*(90/2048)))])
    colormap jet
    c=colorbar('Location','eastoutside');
    axis([0 90 -3 3])
    str = '$$ log( \frac{dE_z}{E_*}/ d log \frac{t_{col}}{t_*} d \theta_{tang} ) $$';
    ylabel(c,str,'Interpreter','latex','FontSize',12) 
    view(2);
    xlabel('\theta_{tang}'); 
    ylabel('Log t_{col}/t_*');
    set(gca,'LineWidth',2,'FontSize',12);
    ax1 = gca; % current axes
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
    x2=linspace(0,90,10);
    y2=linspace(-3,3,10);
    line(x2,y2,'Parent',ax2,'Color','none')
    set(ax2,'LineWidth',2,'FontSize',12);
    ax2.XTick = 0:20:80 ;
    ax2.XTickLabel =strsplit(num2str((1./cos((0:20:80)*pi/180)),2));
    ax2.YColor = 'none';
    xlabel(ax2,'r_{eq}/R_*');
    subplot(2,1,2);
    h=scatter(xx(v_y<0 & rad>1),yy(v_y<0 & rad>1),3,log10(e_z/((pi/4096)*(4/2048))),'filled');hold on
    grid off
    caxis([log10(median(e_z)/((pi/4096)*(4/2048)))-3,log10(max(e_z)/((pi/4096)*(4/2048)))])
    colormap jet
    c=colorbar('Location','eastoutside');
    axis([0 4 0 4])
    str = '$$ log( \frac{dE_z}{E_*}/ d log \frac{r}{R_*} d \phi ) $$';
    ylabel(c,str,'Interpreter','latex','FontSize',12) 
    view(2);
    title(['t/t_*= ' num2str(time.time1(t+1)/ttot) ', E_{tot}=' num2str(totalenergy,3) ]);
    xlabel('x/R_*'); 
    ylabel('y/y_*');
    set(gca,'LineWidth',2,'FontSize',12);
    save(['/home/nilou/Data/processeddata/Collision/collision_' num2str(t) '.mat'],'e_z','X','Y','Vq','t_y','angleth'); 
    name=['/home/nilou/Data/plot/collision/gridcol_v2_' num2str(t) '.png'];
    set(gcf, 'PaperPosition', [0 0 3.5 6]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [3.5 6]); %Set the paper to have width 5 and height 5.
    print(gcf, '-dpng', '-r900', name)
end