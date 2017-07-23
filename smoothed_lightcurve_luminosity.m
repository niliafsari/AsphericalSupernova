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
%luminosity=zeros(1,50);
%load('luminosity.mat', 'luminosity');
% temp=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55.mat','luminosity');
% luminosity=temp.luminosity;
time=load('/home/nilou/Data/timesteps.mat');
luminosity=zeros(1,50);
for t=23:50
    t
    % close all
     clear xdiff
     clear ydiff
     clear Y_new
     clear X_new
    diff_bsg1=load(['/home/nilou/Data/processeddata/RSG/diff_rsg_smoothed55_' num2str(t) '.mat'], 'diff_bsg'); 

     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
     density= csvread(name)*rhoconv;
     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
     pres= csvread(name)*pconv;
     name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-2) '.csv'];
     gradp=csvread(name)*(rconv_bsg/pconv_bsg)*(pconv/rconv);
     if (t==23)
         diff_bsg_t=diff_bsg1.diff_bsg(:,3:length(diff_bsg1.diff_bsg));
     else
         diff_bsg_t=diff_bsg1.diff_bsg(:,2:length(diff_bsg1.diff_bsg));
     end
    diff_bsg_t=transpose(diff_bsg_t);
    diff_bsg=unique(diff_bsg_t,'rows');
    diff_bsg=transpose(diff_bsg);
    rdiff=sqrt(diff_bsg(1,1:length(diff_bsg)).^2+diff_bsg(2,1:length(diff_bsg)).^2);
    xdiff=diff_bsg(1,1:length(diff_bsg));
    ydiff=diff_bsg(2,1:length(diff_bsg)); 
    phidiff=atan(xdiff./ydiff);
    
    rdiff=rdiff(phidiff<=prctile(phidiff,99.2));
    xdiff=xdiff(phidiff<=prctile(phidiff,99.2));
    ydiff=ydiff(phidiff<=prctile(phidiff,99.2));
    phidiff=atan(xdiff./ydiff);
    
    phidiff=phidiff(rdiff>(0.5*rconv));
    xdiff=xdiff(rdiff>(0.5*rconv));
    ydiff=ydiff(rdiff>(0.5*rconv));
    rdiff=rdiff(rdiff>(0.5*rconv));
    
    [phidiff,I]=sort(phidiff);    
    rdiff=rdiff(I);
    xdiff=xdiff(I);
    ydiff=ydiff(I);
    
    
    x=1:length(rdiff);
    x=x';
    fit1 = fit(x,rdiff','poly5');
    fdata = feval(fit1,x);
    I = abs(fdata - rdiff') > 0.6*std(rdiff');
    outliers = excludedata(x,rdiff','indices',I);
    %sum(outliers)
 
    rdiff(outliers)=[];
    xdiff(outliers)=[];
    ydiff(outliers)=[];
    phidiff=atan(xdiff./ydiff);
    
        
    x=1:length(phidiff);
    x=x';
    fit1 = fit(x,phidiff','poly5');
    fdata = feval(fit1,x);
    I = abs(fdata - phidiff') > 0.6*std(phidiff');
    outliers = excludedata(x,phidiff','indices',I);
    %sum(outliers)
 
    rdiff(outliers)=[];
    xdiff(outliers)=[];
    ydiff(outliers)=[];
    phidiff=atan(xdiff./ydiff);
    
    [xmax,I]=max(xdiff);
    sum(phidiff>phidiff(I))
    rdiff(phidiff>phidiff(I))=[];
    xdiff(phidiff>phidiff(I))=[];
    ydiff(phidiff>phidiff(I))=[];
    phidiff=atan(xdiff./ydiff);    
        
    
    index_r=floor(rdiff/((2/2048)*rconv));
    index_phi=floor(atan(xdiff./ydiff)/(0.5*pi/2048));
    index_phi(index_phi==0)=1;
    index_r(index_r==0)=1;
    if  sum(index_r>2048)
        sum(index_r>2048)
    end
    
    index_r(index_r>2048)=2048;
    flux1=zeros(1,length(rdiff));
    factor=zeros(1,length(rdiff));
    luminosity(t)=0;
    
    for k=2:length(rdiff)
         flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
         factor(k)=4*pi* rdiff(k)^2 * (sin(phidiff(k))^2-sin(phidiff(k-1))^2);
    end
    
    luminosity(1,t)=dot(flux1,factor);
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
       
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    scatter(xdiff*2/rconv,ydiff*2/rconv,'r');
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
    axis equal
    axis([0 4 0 4])
    name=['/home/nilou/Data/plot/RSGdiff/RSG_filter_smoothed_temp_' num2str(t) '.png'];
    print(gcf, '-dpng', '-r50', name)
    export_fig(name, '-dpng', '-r50')

end
save('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55_Phi0.mat','luminosity')
