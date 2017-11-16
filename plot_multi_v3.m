mtot=0.0048187313*2;
etot=0.0130949665*2;
rtot=0.5;


ptot=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

msun=1.989e33;
rsun=6.955e10;

m=5*msun;
e=1e51;
r=0.2*rsun;

kappa=0.2;
c=3e10;
h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;
Z=0.005;
a=7.566e-15;

rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

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

close all



    x0=0;
    y0=0;
    width=650;
    height=550;
    myFigure = figure ;    
    set(myFigure,'units','points','position',[299.2500   66.0000  609.7500  519.0000],'Color','w') 
path='/home/nilou/Data';
for t=234
    load([path '/processeddata/ic/diff_ic_1024_' num2str(t-1) '.mat'], 'diff_bsg'); 
    diff_bsg_t=diff_bsg(:,2:length(diff_bsg));
    diff_bsg_t=transpose(diff_bsg_t);
    %length(diff_bsg_t)
    diff_bsg=unique(diff_bsg_t,'rows');
    length(diff_bsg)
    diff_bsg=transpose(diff_bsg);
    rdiff_n=sqrt(diff_bsg(1,1:length(diff_bsg)).^2+diff_bsg(2,1:length(diff_bsg)).^2);
    xdiff_n=diff_bsg(1,1:length(diff_bsg));
    ydiff_n=diff_bsg(2,1:length(diff_bsg)); 
    phidiff_n=atan(xdiff_n./ydiff_n);
    
    size(rdiff_n)
    rdiff_n=rdiff_n(phidiff_n<=prctile(phidiff_n,98));
    xdiff_n=xdiff_n(phidiff_n<=prctile(phidiff_n,98));
    ydiff_n=ydiff_n(phidiff_n<=prctile(phidiff_n,98));
    phidiff_n=atan(xdiff_n./ydiff_n);
    
    if t>173
        phidiff_t=phidiff_n(rdiff_n>(0.5*rconv));
        xdiff_t=xdiff_n(rdiff_n>(0.5*rconv));
        ydiff_t=ydiff_n(rdiff_n>(0.5*rconv));
        rdiff_t=rdiff_n(rdiff_n>(0.5*rconv));
    else
        phidiff_t=phidiff_n;
        xdiff_t=xdiff_n;
        ydiff_t=ydiff_n;
        rdiff_t=rdiff_n;
    end

    [phidiff,I]=sort(phidiff_t);
    rdiff=rdiff_t(I);
    
    xdiff=xdiff_t(I);
    ydiff=ydiff_t(I);    
    phidiff=atan(xdiff./ydiff);
    if t>183
        x=1:length(rdiff);
        x=x';
        fit1 = fit(x,rdiff','poly6');
        fdata = feval(fit1,x);
        I = abs(fdata - rdiff') > 0.25*std(rdiff');
        outliers = excludedata(x,rdiff','indices',I);
        sum(outliers)

        rdiff(outliers)=[];
        xdiff(outliers)=[];
        ydiff(outliers)=[];
        phidiff=atan(xdiff./ydiff);


          [xmax,I]=max(xdiff);
    %     sum(phidiff>phidiff(I))

          rdiff_k=rdiff(phidiff>phidiff(I));
          xdiff_k=xdiff(phidiff>phidiff(I));
          ydiff_k=ydiff(phidiff>phidiff(I));
          phidiff_k=phidiff(phidiff>phidiff(I));
          rdiff=rdiff(phidiff<=phidiff(I));
          xdiff=xdiff(phidiff<=phidiff(I));
          ydiff=ydiff(phidiff<=phidiff(I));
          phidiff=phidiff(phidiff<=phidiff(I));

            P=[xdiff_k', ydiff_k'];
            if (length(P)>0)
                DT = delaunayTriangulation(P);
                Q = convexHull(DT);
                if (length(P)>1)
                    Q(1:2)=[];
                end
                Q=flipud(Q);
                rdiff_k=rdiff_k(Q);
                xdiff_k=xdiff_k(Q);
                ydiff_k=ydiff_k(Q);
                phidiff_k=phidiff_k(Q);
            else
                rdiff_k=[];
                xdiff_k=[];
                ydiff_k=[];
                phidiff_k=[];
            end

        rdiff=[rdiff rdiff_k];
        xdiff=[xdiff xdiff_k];
        ydiff=[ydiff ydiff_k];
        phidiff=[phidiff phidiff_k];
    end
        [xmax,I]=max(xdiff);
        [ymax,U]=max(ydiff);

%         name=[path '/rawdata/pressure/pres1024_' int2str(t-1) '.mat'] ;
%         load(name,'pres_data');
%         pres=pres_data/ptot;
%         name=[path '/rawdata/density/dens1024_' int2str(t-1) '.mat'] ;
%         load(name,'dens_data');
%         density=dens_data/rhotot;
%         name=[path '/rawdata/velocity/velr1024_' int2str(t-1) '.mat'] ;  
%         load(name,'velx_data');
%         velr= velx_data/vtot;
%         name=[path '/rawdata/velocity/velth1024_' int2str(t-1) '.mat'] ;
%         load(name,'vely_data');
%         velt= vely_data/vtot;
%     vel=sqrt(velr.^2.+velt.^2);
%     entropy=pres./density.^(4/3);
%     dens_l = griddata(xx,yy,dens,xx_l,yy_l);
%     pres_l = griddata(xx,yy,pres,xx_l,yy_l);
    
    a=get(gcf,'Position');
    close all

     myFigure = figure ;    
     set(myFigure,'units','points','position',a,'Color','white')  
    
    %**********
    

    
    s(4)=subplot(2,2,3);
    h(4)=surf(-xx/rtot,-yy/rtot, (vel));
    hold on
    alpha(0.7)
    scatter(-xdiff_n*2/rconv,-ydiff_n*2/rconv,1,'c'); 
    axis equal
    grid off
    set(h(4),'LineStyle','none');
    colormap jet
    c4=colorbar('Location','westoutside'); 
    caxis([0 25])
    str = 'v/v_*';
    ylabel(c4,str,'FontSize',12) 
    view(2);
    set(gca,'LineWidth',2,'FontSize',12);
    s(4).XAxis.TickValues=[-4 -3 -2 -1 0];
    s(4).YAxis.TickValues=[-4 -3 -2 -1 0];
%     
    s(3)=subplot(2,2,1);
    h(3)=surf(-xx/rtot,yy/rtot, log10(pres));
    
    grid off
    set(h(3),'LineStyle','none');
    hold on
    alpha(0.7)
    scatter(-xdiff_n*2/rconv,ydiff_n*2/rconv,1,'c'); 
    axis equal
    %colormap jet
    c3=colorbar('Location','westoutside'); 
    caxis([-9 0])
    str = 'Log P/P_*';
    ylabel(c3,str,'FontSize',12) 
    view(2);
    set(gca,'LineWidth',2,'FontSize',12);
    s(3).XAxis.Color='none';
    s(3).YAxis.TickValues=[0 1 2 3 4];
    s(3).XAxis.TickValues=[-4 -3 -2 -1 0];
    s(1)=subplot(2,2,4);
    %pos=s(1).Position;
    h(1)=surf(xx/rtot,-yy/rtot, log10(entropy));hold on
    grid off
    set(h(1),'LineStyle','none');
        hold on
    alpha(0.7)
    scatter(xdiff_n*2/rconv,-ydiff_n*2/rconv,1,'c'); 
    axis equal
    c1=colorbar('Location','eastoutside');
    caxis([-8 4])
    str = 'Log S/S_*';
    ylabel(c1,str,'FontSize',12) 
    colormap jet
    view(2);
    set(gca,'LineWidth',2,'FontSize',12);
    s(1).YAxis.Color='none';
    s(1).XAxis.TickValues=[0 1 2 3 4];

    
    s(2)=subplot(2,2,2);
    h(2)=surf(xx /rtot,yy /rtot, log10(density));
    grid off
    set(h(2),'LineStyle','none');
        hold on
    alpha(0.7)
    scatter(xdiff_n*2/rconv,ydiff_n*2/rconv,1,'c'); 
    axis equal
    colormap jet
    c2=colorbar('Location','eastoutside'); 
    caxis([-8.5 1])
    str = 'Log \rho/ \rho_*';
    ylabel(c2,str,'FontSize',12) 
    view(2);
    set(gca,'LineWidth',2,'FontSize',12);
    s(2).YAxis.Color='none';
    s(2).XAxis.Color='none';
    %s(2).XAxis.TickValues=[0 1 2 3 4];
    s(3).YAxis.TickValues=[0 1 2 3 4];
   s(1).Position=[0.15+0.322 0.1 0.38 0.38];
    s(2).Position=[0.15+0.322 0.1+0.38 0.38 0.38];
    s(4).Position=[0.15 0.1 0.38 0.38];
     s(3).Position=[0.15 0.1+0.38 0.38 0.38];
    
    
    c3.Position=[c4.Position(1)   c4.Position(2)+s(4).Position(3) c4.Position(3) c4.Position(4)];
    c2.Position=[c1.Position(1)   c1.Position(2)+s(4).Position(3) c1.Position(3) c1.Position(4)];
    set(gcf,'NextPlot','add');
    axes;
    %h = title(['t/t_*=' num2str(time.time1(t+1)/ttot,3)]);
    num2str(time.time1(t+1)/ttot,3)
    set(gca,'Visible','off');
%    set(h,'Visible','on'); 
%    set(h,'FontSize',10);
    %xlabel('R/R_*')
    
%     set(gca,'LineWidth',2,'FontSize',12);
%     (time.time1(t+1)/ttot)
     name=['/home/nilou/Data/plot/allquantity/1024_' num2str(t,'%02d') '.pdf'];
     % set(gcf, 'PaperPosition', [0 0 6.5 5.5]); %Position plot at left hand corner with width 5 and height 5.
     % set(gcf, 'PaperSize', [6.5 5.5]); %Set the paper to have width 5 and height 5.
     %print(gcf, '-dpng', '-r1100', name)
     %print(gcf, '-dpdf','-r100',name) 
     %export_fig(name, '-pdf','-r100')
     

%      name=['/home/nilou/Data/plot/allquantity/pdf/all_v4_' num2str(t,'%02d') '.pdf'];
%      print('-dpdf',name) 
%      export_fig(name, '-pdf')
     
end

