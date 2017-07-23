mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=15*msun;
e=1e51;
r=49*rsun;

Z=0.005;
a=7.566e-15;
XX=0.7;
kappa_T=0.2*(1+XX);

c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/(rconv^3);
rhoconv=mconv/(rconv^3);
vconv=sqrt(econv/mconv);

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

radius=zeros(1,2048);
[X,Y] = meshgrid(x,y);
theta_car=asin(X./sqrt(X.^2+Y.^2));
for t=35:71 
    name1=['/home/nilou/Data/rawdata/dcodeunit/dcodeunit_' int2str(t-1) '.csv'] ;
    d=csvread(name1);
    indexes=find(isnan(d) | isinf(d));
    indexes_tot=indexes((xx(indexes).^2+yy(indexes).^2)>0.5^2);
    name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
    pres= csvread(name).*pconv;
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
    dens= csvread(name).*rhoconv;
    dens(indexes_tot)=0;
    name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t-1) '.csv'] ;
    velr= csvread(name).*vconv;
    name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t-1) '.csv'] ;
    velt= csvread(name).*vconv;    
    
    dens_car = real(griddata(xx,yy,dens,X,Y, 'linear'));
    pres_car = real(griddata(xx,yy,pres,X,Y, 'linear'));
    velr_car = real(griddata(xx,yy,velr,X,Y, 'linear'));
    velt_car = real(griddata(xx,yy,velt,X,Y, 'linear'));
    vy_car=-velt_car.*sin(theta_car) + velr_car.*cos(theta_car);
    vx_car=velt_car.*cos(theta_car) + velr_car.*sin(theta_car);
      %vx_car=sqrt(velt_car.^2+velr_car.^2).*(Y./sqrt(X.^2+Y.^2));
    dens_car(isnan(dens_car))=0;
    pres_car(isnan(pres_car))=0;
    a=7.566e-15;
    kappa_ff=3.7e22*(1+XX)*(1-Z)*dens_car.*((3*pres_car/a).^-(7/8));
    kappa_ff(isnan(kappa_ff))=0;
    kappa_ff(isinf(kappa_ff))=0;
    kappa=kappa_T+kappa_ff;
    dtau_ff=kappa_ff.*dens_car*2/2048*rconv;
    dtau=kappa.*dens_car*2/2048*rconv;
    for i=1:2048
        tau(i,:)=sum(dtau(i:2048,:),1);
        tau_ff(i,:)=sum(dtau_ff(i:2048,:),1); 
    end
    name=[ '/home/nilou/Data/processeddata/ic/tauBSGrY_' int2str(t-1) '.csv'];
    nameff=[ '/home/nilou/Data/processeddata/ic/tauBSGYFF_' int2str(t-1) '.csv'];
    csvwrite(name,tau);
    csvwrite(nameff,tau_ff);
    
    close;
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
    
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    h=surf(X/0.5,Y/0.5, real(vy_car)/c);hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    cl=colorbar;
    axis([0 4 0 4])
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t/t_*=' num2str(time.time1(t),3)]); 
    axis equal   
    C=contour3(X/0.5,Y/0.5,(3*vy_car.*tau)/c,[1 1],'LineColor','k','Linewidth',1);
    set(gca,'LineWidth',2,'FontSize',12);
    caxis([-0.05 0.05])   
    ylabel(cl,'$ \frac{v_\perp}{c} $','interpreter','latex','FontSize',24);    
    view(2); 
    name=['/home/nilou/Data/plot/collision/colDiffFrontBSG_' num2str(t) '.png'];
    print(gcf, '-dpng', '-r100', name)
    export_fig(name, '-dpng', '-r100')
end