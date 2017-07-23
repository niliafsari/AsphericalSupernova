mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=14*msun;
e=1e51;
r=400*rsun;
kappa=0.34;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;


theta=linspace(0,pi/2,2048);


radius=linspace(0,2,2048)*rconv;
r=linspace(0,2,2048).*rconv;
for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
gradp=zeros(2048,2048);
luminosity=zeros(1,50);
%load('luminosity.mat', 'luminosity');
time=load('/home/nilou/Data/timesteps.mat');
clear diff_bsg
for t=23:50
    t
    name=['/home/nilou/Data/processeddata/RSG/dparamRSG_' int2str(t-1) '.csv'] ;
    diff_frontBSG= csvread(name);
%     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
%     density= csvread(name)*rhoconv;
%     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
%     pres= csvread(name)*pconv;
%     name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-1) '.csv'];
%     gradp=csvread(name);
    if (t==23)
        ke=1;
    elseif (t>23 & t<28)
        ke=10;
    elseif (t>27 & t<31)
        ke=25;
    elseif (t>30 & t<33)
        ke=45;
    elseif (t>32 & t<40)
        ke=55;
    elseif (t>39 & t<51)
        ke=55;
    end
    temp=zeros(2048+2*ke,2048+ke);
    temp(1:2048,1:2048)=diff_frontBSG;
    kernel = 1*fspecial('disk', ke);
    temp1 = imfilter(log10(real(temp)),kernel,'same');
    dsmooth= temp1(1:2048,1:2048);

%     ke=40;
%     kernel = 1*fspecial('disk', ke);
%     dsmooth = imfilter(log10(real(diff_frontBSG)),kernel,'same');

      
    [diff_bsg,~]=contour3(xx,yy,dsmooth,[0 0],'LineColor','b','Linewidth',0.3);hold on  
    diff_bsg
%     rdiff=sqrt(diff_bsg{t}(1,2:length(diff_bsg{t})).^2+diff_bsg{t}(2,2:length(diff_bsg{t})).^2);
%     xdiff=diff_bsg{t}(1,2:length(diff_bsg{t}));
%     ydiff=diff_bsg{t}(2,2:length(diff_bsg{t})); 
%     phidiff=atan(xdiff./ydiff);
%     [phidiff,I]=sort(phidiff);
%     rdiff=rdiff(I);
%     xdiff=xdiff(I);
%     ydiff=ydiff(I);
%     index_r=floor(rdiff/((2/2048)*rconv));
%     index_phi=floor(atan(xdiff./ydiff)/(0.5*pi/2048));
%     index_phi(index_phi==0)=1;
%     index_r(index_r==0)=1;
%     if  sum(index_r>2048)
%         sum(index_r>2048)
%     end
%     index_r(index_r>2048)=2048;
%     flux1=zeros(1,length(diff_bsg{t})-1);
%     for k=2:length(index_r)
%         flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
% %         if (k==1)
% %             luminosity(t)= luminosity(t)+(2*pi* rdiff(k)^2 * flux1(k) *(cos(phidiff(k))-cos(phidiff(k+1))));
% %         else
%             luminosity(t)= luminosity(t)+(2*pi* rdiff(k)^2 * flux1(k) *(cos(phidiff(k-1))-cos(phidiff(k))));
%         end
%     end
   close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
%     
    
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    h=surf(xx,yy,log10(real(diff_frontBSG)));hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    cl=colorbar;
    %axis([0 4 0 4])
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']); 
    axis equal       
    
    %contour3(xx,yy,log10(real(diff_frontBSG)),[0 0],'LineColor','k','Linewidth',1);
    contour3(xx,yy,dsmooth,[0 0],'LineColor','r','Linewidth',1);
    set(gca,'LineWidth',2,'FontSize',12);
    caxis([-2 12]) 
    ylabel(cl,'$ Log \frac{t_{diff}}{t_{dyn}} $','interpreter','latex','FontSize',24);    
    view(2); 
    name=['/home/nilou/Data/plot/RSGdiff/DiffFrontRSG_smoothed_' num2str(t) '.png'];
    print(gcf, '-dpng', '-r100', name)
    export_fig(name, '-dpng', '-r100')
     save(['/home/nilou/Data/processeddata/RSG/diff_rsg_smoothed55_' num2str(t) '.mat'], 'diff_bsg') 
end

