mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=15*msun;
e=1e51;
r=49*rsun;
kappa=0.34;
c=3e10;
rconv_bsg=r/rtot;
econv_bsg=e/etot;
mconv_bsg=m/mtot;
pconv_bsg=econv_bsg/rconv_bsg^3;
rhoconv_bsg=mconv_bsg/rconv_bsg^3;
vconv_bsg=sqrt(econv_bsg/mconv_bsg);

m=5*msun;
e=1e51;
r=0.2*rsun;

Z=0.005;
a=7.566e-15;
XX=0;
kappa_T=0.2*(1+XX);

c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;


pconv=econv/(rconv^3);
rhoconv=mconv/(rconv^3);
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048);
x=linspace(0,2,2048);
y=linspace(0,2,2048);

time=load('/home/nilou/Data/timesteps.mat');
%theta=linspace(0,pi/2,2048);
for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end

[X,Y] = meshgrid(x,y);
theta_car=asin(X./sqrt(X.^2+Y.^2));
cone=[0.02 0.02 0.1 0.1 0.1 0.15 0.17 0.17 0.2 0.2 0.22 0.23 0.23 0.24 0.25 0.27 0.3 0.32 ...
    0.33 0.34 0.34 0.34 0.34 0.34 0.34 0.35 0.36 0.36 0.36 0.37 0.37 0.38 0.39  0.4 0.4 0.41 0.42 0.42 0.42];
cone=cone*rtot;
clear xdiff_col_p xdiff_col_q xdiff_col ydiff_col_p ydiff_col_q ydiff_col
for t=41:41 
      clear xdiff_col_p xdiff_col_q xdiff_col ydiff_col_p ydiff_col_q ydiff_col
%     name1=['/home/nilou/Data/rawdata/dcodeunit/dcodeunit_' int2str(t-1) '.csv'] ;
%     d=csvread(name1);
%     indexes=find(isnan(d) | isinf(d));
%     indexes_tot=indexes((xx(indexes).^2+yy(indexes).^2)>1^2);
%     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
%     pres= csvread(name).*pconv;
%     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
%     dens= csvread(name).*rhoconv;
%     dens(indexes_tot)=0;      
%     dens_car = real(griddata(xx,yy,dens,X,Y, 'linear'));
%     pres_car = real(griddata(xx,yy,pres,X,Y, 'linear'));
%     dens_car(isnan(dens_car))=0;
%     pres_car(isnan(pres_car))=0;
%     [gradpx,gradpy]=gradient(pres_car,2*rconv/2048);
%     tdiff=(3*kappa_T*dens_car.*(pres_car.^2))./(c*gradpy.^2);
%     tdyn=(time.time1(t)-time.time1(23))*tconv;
%     
%     D=tdiff/tdyn;
%     name=['/home/nilou/Data/processeddata/colDiffFrontIc_v2_' num2str(t) '.mat'];
    save(name, 'D')
    [col_contourIc,~]=contour3(X,Y,real(tdiff/tdyn),[1 1],'LineColor','k','Linewidth',1);
    xdiff_col_q=col_contourIc(1,2:length(col_contourIc));
    ydiff_col_q=col_contourIc(2,2:length(col_contourIc));
    
     xdiff_col_p=xdiff_col_q(ydiff_col_q<cone(t-34));
     ydiff_col_p=ydiff_col_q(ydiff_col_q<cone(t-34));
    
    xdiff_col=xdiff_col_p(xdiff_col_p>0.5);
    ydiff_col=ydiff_col_p(xdiff_col_p>0.5);
%     
    [xdiff_col,I]=sort(xdiff_col);
    ydiff_col=ydiff_col(I);
    
    %X=X(X==xdiff_col);
    dv= (2/2048)^2*2*pi*sqrt(X.^2+Y.^2)*rconv^3;
     for k=1:length(xdiff_col)         
%         energy(k)= 3*pres_car(Y<ydiff_col(k) & X==xdiff_col(k)).*dv(Y<ydiff_col(k) & X==xdiff_col(k));    
%     end
    
    close;
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
    
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    h=surf(X,Y,real(log10(tdiff/tdyn)));hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    cl=colorbar;
    axis([0 2 0 2])
    %contour3(X,Y,real(tdiff/tdyn),[1 1],'LineColor','k','Linewidth',1);hold on
    %scatter(xdiff_col,ydiff_col,3,'k','filled');hold on
    yscatter=cone(t-34)*ones(1,2048);
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t/t_*=' num2str(time.time1(t),3)]); 
    axis equal   
    set(gca,'LineWidth',2,'FontSize',12);
    caxis([-2 12]) 
    ylabel(cl,'$ Log \frac{t_{diff}}{t_{dyn}} $','interpreter','latex','FontSize',24);    
    view(2); 
    name=['/home/nilou/Data/plot/collision/colDiffFrontIc_v3' num2str(t) '.png'];
    print(gcf, '-dpng', '-r100', name)
    export_fig(name, '-dpng', '-r100')
end