mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=14*msun;
e=1e51;
r=400*rsun;
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

for t=39:39
    radius=zeros(1,2048);
    name=['/home/nilou/Data/processeddata/RSG/tauRSGradial_' int2str(t) '.csv'] ;
    tau= csvread(name);
    nameff=['/home/nilou/Data/processeddata/RSG/tauRSGradialFF_' int2str(t) '.csv'] ;
    tauff= csvread(nameff);
    diff_bsg1=load(['/home/nilou/Data/processeddata/RSG/diff_rsg_smoothed55_' num2str(t+1) '.mat'], 'diff_bsg'); 

     if (t==23)
         diff_bsg_t=diff_bsg1.diff_bsg(:,3:length(diff_bsg1.diff_bsg));
     else
         diff_bsg_t=diff_bsg1.diff_bsg(:,2:length(diff_bsg1.diff_bsg));
     end
    diff_bsg_t=transpose(diff_bsg_t);
    length(diff_bsg_t)
    diff_bsg=unique(diff_bsg_t,'rows');
    length(diff_bsg)
    diff_bsg=transpose(diff_bsg);
    rdiff=sqrt(diff_bsg(1,1:length(diff_bsg)).^2+diff_bsg(2,1:length(diff_bsg)).^2);
    xdiff=diff_bsg(1,1:length(diff_bsg));
    ydiff=diff_bsg(2,1:length(diff_bsg)); 
%     phidiff=atan(xdiff./ydiff);
%     
%     rdiff=rdiff(phidiff>prctile(phidiff,0.8));
%     xdiff=xdiff(phidiff>prctile(phidiff,0.8));
%     ydiff=ydiff(phidiff>prctile(phidiff,0.8));
%     phidiff=atan(xdiff./ydiff);
%     
%      phidiff=phidiff(rdiff>(0.5*rconv));
%      xdiff=xdiff(rdiff>(0.5*rconv));
%      ydiff=ydiff(rdiff>(0.5*rconv));
%      rdiff=rdiff(rdiff>(0.5*rconv));
%     
%     [phidiff,I]=sort(phidiff);    
%     rdiff=rdiff(I);
%     xdiff=xdiff(I);
%     ydiff=ydiff(I);
%     
%     
%     xxx=1:length(rdiff);
%     xxx=xxx';
%     fit1 = fit(xxx,rdiff','poly5');
%     fdata = feval(fit1,xxx);
%     I = abs(fdata - rdiff') > 0.6*std(rdiff');
%     outliers = excludedata(xxx,rdiff','indices',I);
%     sum(outliers)
%  
%     rdiff(outliers)=[];
%     xdiff(outliers)=[];
%     ydiff(outliers)=[];
%     phidiff=atan(xdiff./ydiff);
%     
% 
%     
%     xxx=1:length(phidiff);
%     xxx=xxx';
%     fit1 = fit(xxx,phidiff','poly5');
%     fdata = feval(fit1,xxx);
%     I = abs(fdata - phidiff') > 0.6*std(phidiff');
%     outliers = excludedata(xxx,phidiff','indices',I);
%     sum(outliers)
%  
%     rdiff(outliers)=[];
%     xdiff(outliers)=[];
%     ydiff(outliers)=[];
%     phidiff=atan(xdiff./ydiff);
    
    name=['/home/nilou/Data/processeddata/RSG/dparamRSG_' int2str(t) '.csv'] ;
    d= csvread(name);
    tau_eff=sqrt(3*tau.*tauff);
    close;
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
    
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    p3=plot(x,y,'r','Linewidth',1); hold on
    h1=contour3(xx,yy,tau_eff,[2/3 2/3],'LineColor','b','Linewidth',1);hold on  
    h1=contour3(xx,yy,tau,[1 1],'LineColor','m','Linewidth',1);hold on  
 %   axis equal 
    h2=scatter(xdiff,ydiff,1,'g');hold on
    h=surf(xx,yy, log10(tau_eff));hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    alpha(.3)
    colorbar
    caxis([-4 2])
    view(2);
    view(2);
    axis([0,2*rconv,0,2*rconv])
    title(['RSG Model t=' int2str(t)]);
    xlabel('cm'); 
    ylabel('cm');
    legend('Initial Progenitor','Thermalization Depth','Photosphere','Diffusion Front','Log \tau_{eff}','Location','southoutside','Orientation','horizontal');
    set(gca,'LineWidth',2,'FontSize',12);
    %set(gca,'XDir','reverse')
    name=['/home/nilou/Data/plot/RSGdiff/tauRSGContourWdiff_' int2str(t)];
    print('-dpng',name) 
    export_fig(name, '-png')
end