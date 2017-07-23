mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=5*msun;
e=1e51;
r=0.2*rsun;
%kappa=0.2;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);

%r=linspace(0,2,2048)*rconv;
r=0.5*rconv;
theta=linspace(0,pi/2,2048);

x=r*sin(theta);
y=r*cos(theta);

radius=linspace(0,2,2048)*rconv;

%[TH,R] = meshgrid(theta,radius);
%[X,Y] = pol2cart(TH,R);

% xx=linspace(0,2,2048)*rconv;
% yy=linspace(0,2,2048)*rconv;
% [xg,yg] = meshgrid(x,y);

for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end

for t=28:28
    radius=zeros(1,2048);
    name=['/home/nilou/Data/processeddata/ic/tauIcradial_' int2str(t) '.csv'] ;
    tau= csvread(name);
    nameff=['/home/nilou/Data/processeddata/ic/tauIcradialFF_' int2str(t) '.csv'] ;
    tauff= csvread(name);
    name=['/home/nilou/Data/processeddata/ic/dparam_' int2str(t) '.csv'] ;
    d= csvread(name);
    d=d;
    tau=sqrt(3.*tau.*tauff);
    close;
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
    
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    h=surf(xx,yy, log10(tau));hold on
    grid off
    set(h,'LineStyle','none');
    colormap(flipud(gray))
    colorbar
    caxis auto
    view(2);
    p3=plot(x,y,'r','Linewidth',1); hold on
    h1=contour3(xx,yy,tau,[2/3 2/3],'LineColor','b','Linewidth',1);hold on  
 %   axis equal 
    h2=contour(xx,yy,d,[1 1],'LineColor','g');
    axis([0,2*rconv,0,2*rconv])
    title(['Ic Model t=' int2str(t)]);
    xlabel('cm'); 
    ylabel('cm');
    legend('Log \tau_{eff}','Initial Progenitor','Thermalization Depth','Diffusion Front','Location','southoutside','Orientation','horizontal');
    set(gca,'LineWidth',2,'FontSize',12);
    %set(gca,'XDir','reverse')
    name=['/home/nilou/Data/plot/icdiff/tauIcContourWdiff_' int2str(t)];
    print('-dpng',name) 
    export_fig(name, '-png')
end