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
kappa=0.2;
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

for t=23:32
    
    t
    diff_bsg1=load(['/home/nilou/Data/processeddata/ic/diff_ic_smoothed_' num2str(t) '.mat']);
     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
     density= csvread(name)*rhoconv;
     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
     pres= csvread(name)*pconv;
     name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-1) '.csv'];
     gradp=csvread(name)*(rconv_bsg/pconv_bsg)*(pconv/rconv);
     if (t==23)
         diff_bsg_t=diff_bsg1.diff_bsg(:,3:length(diff_bsg1.diff_bsg));
     else
         diff_bsg_t=diff_bsg1.diff_bsg(:,3:length(diff_bsg1.diff_bsg));
     end
    diff_bsg_t=transpose(diff_bsg_t);
    diff_bsg=unique(diff_bsg_t,'rows');
    diff_bsg=transpose(diff_bsg);
    rdiff_n=sqrt(diff_bsg(1,1:length(diff_bsg)).^2+diff_bsg(2,1:length(diff_bsg)).^2);
    xdiff_n=diff_bsg(1,1:length(diff_bsg));
    ydiff_n=diff_bsg(2,1:length(diff_bsg)); 
    phidiff_n=atan(xdiff_n./ydiff_n);
    
    [phidiff_n,I]=sort(phidiff_n);
    rdiff_n=rdiff_n(I);
    xdiff_n=xdiff_n(I);
    ydiff_n=ydiff_n(I);    
    
    rdiff=rdiff_n(rdiff_n>(0.5*rconv));
    phidiff=phidiff_n(rdiff_n>(0.5*rconv));
    xdiff=xdiff_n(rdiff_n>(0.5*rconv));
    ydiff=ydiff_n(rdiff_n>(0.5*rconv));
    

     rdiff=rdiff(phidiff<=prctile(phidiff,99));
     xdiff=xdiff(phidiff<=prctile(phidiff,99));
     ydiff=ydiff(phidiff<=prctile(phidiff,99));
     phidiff=atan(xdiff./ydiff);
    
        
    x=1:length(rdiff);
    x=x';
    fit1 = fit(x,rdiff','poly5');
    fdata = feval(fit1,x);
    I = abs(fdata - rdiff') > 0.6*std(rdiff');
    outliers = excludedata(x,rdiff','indices',I);
    sum(outliers)
    rdiff(outliers)=[];
    xdiff(outliers)=[];
    ydiff(outliers)=[];
    phidiff=atan(xdiff./ydiff);
    
    
    [xmax,I]=max(xdiff);
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
    for k=2:length(index_r)
         densityk(k)=density(index_r(k),index_phi(k));
         gradpk(k)=gradp(index_r(k),index_phi(k));
         flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
         factor(k)=4*pi* rdiff(k)^2 * (sin(phidiff(k))^2-sin(phidiff(k-1))^2);
         %luminosity(t)= luminosity(t)+(2*pi* rdiff(k)^2 * flux1(k) *(cos(phidiff(k-1))-cos(phidiff(k))));
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
    h=surf(xx*2/rconv,yy*2/rconv,log10(gradp)); hold on
    colorbar 
    colormap jet
    alpha(.3)
    caxis([0 10])
    grid off
    set(h,'LineStyle','none');
    view(2)
    scatter(xdiff*2/rconv,ydiff*2/rconv,2,'r'); hold on
%     scatter(xdiff(index_k)*2/rconv,ydiff(index_k)*2/rconv,3,'g'); hold on
    clear index_k
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
    axis equal
    axis([0 4 0 4])
    name=['/home/nilou/Data/plot/icdiff/ic_filter_smoothed_phi0_' num2str(t) '.png'];
    print(gcf, '-dpng', '-r50', name)
    export_fig(name, '-dpng', '-r50')

end
save('/home/nilou/Data/processeddata/ic/luminosity_smoothed_ic_phi0.mat','luminosity')
