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
% temp=load('/home/nilou/Data/processeddata/BSG/luminosity_smoothedk55_temp.mat','luminosity');
%  luminosity=temp.luminosity;
 time=load('/home/nilou/Data/timesteps.mat');

for t=23:42
    
    t
    close all
     clear xdiff
     clear ydiff
     clear Y_new
     clear X_new
    diff_bsg1=load(['/home/nilou/Data/processeddata/BSG/diff_bsg_smoothed_k55_temp' num2str(t) '.mat']);
     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
     density= csvread(name)*rhoconv;
     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
     pres= csvread(name)*pconv;
     name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-2) '.csv'];
     gradp=csvread(name);
     if (t==23)
         diff_bsg_t=diff_bsg1.diff_bsg(:,3:length(diff_bsg1.diff_bsg));
     else
         diff_bsg_t=diff_bsg1.diff_bsg(:,2:length(diff_bsg1.diff_bsg));
     end
    diff_bsg_t=transpose(diff_bsg_t);
    %length(diff_bsg_t)
    diff_bsg=unique(diff_bsg_t,'rows');
    length(diff_bsg)
    diff_bsg=transpose(diff_bsg);
    rdiff_n=sqrt(diff_bsg(1,1:length(diff_bsg)).^2+diff_bsg(2,1:length(diff_bsg)).^2);
    xdiff_n=diff_bsg(1,1:length(diff_bsg));
    ydiff_n=diff_bsg(2,1:length(diff_bsg)); 
    phidiff_n=atan(xdiff_n./ydiff_n);
    

    rdiff_n=rdiff_n(phidiff_n<=prctile(phidiff_n,94));
    xdiff_n=xdiff_n(phidiff_n<=prctile(phidiff_n,94));
    ydiff_n=ydiff_n(phidiff_n<=prctile(phidiff_n,94));
    phidiff_n=atan(xdiff_n./ydiff_n);
    
    rdiff_t=rdiff_n(rdiff_n>(0.5*rconv));
    phidiff_t=phidiff_n(rdiff_n>(0.5*rconv));
    xdiff_t=xdiff_n(rdiff_n>(0.5*rconv));
    ydiff_t=ydiff_n(rdiff_n>(0.5*rconv));
    
    [phidiff,I]=sort(phidiff_t);
    rdiff=rdiff_t(I);
    xdiff=xdiff_t(I);
    ydiff=ydiff_t(I);    
  
    
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
    
%     
%      [ymax,I]=max(ydiff);
%      sum(phidiff>phidiff(I))
%      rdiff(phidiff>phidiff(I))=[];
%      xdiff(phidiff>phidiff(I))=[];
%      ydiff(phidiff>phidiff(I))=[];
%      phidiff=atan(xdiff./ydiff);  
    
    index_r=floor(rdiff/((2/2048)*rconv));
    index_phi=floor(atan(xdiff./ydiff)/(0.5*pi/2048));
    index_phi(index_phi==0)=1;
    index_r(index_r==0)=1;
    if  sum(index_r>2048)
        sum(index_r>2048)
    end
    index_r(index_r>2048)=2048;
    flux1=zeros(1,length(index_r));
    factor=zeros(1,length(index_r));
    densityk=zeros(1,length(index_r));
    gradpk=zeros(1,length(index_r));
    %luminosity(t)=0;
    for k=2:length(index_r)
         densityk(k)=density(index_r(k),index_phi(k));
         gradpk(k)=gradp(index_r(k),index_phi(k));
         flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
         rdiffm=(rdiff(k)+rdiff(k-1))/2;
         dl=sqrt((rdiffm*(phidiff(k)-phidiff(k-1)))^2+(rdiff(k)-rdiff(k-1))^2);
         rn=rdiffm*(phidiff(k)-phidiff(k-1))/dl;
         dz=ydiff(k-1)-ydiff(k);
         dR=xdiff(k)-xdiff(k-1);
         thetan=atan(dz/dR);
         if ( rn==0  || dz<0 )
             factor(k)=0;
         else    
            factor(k)=(8*rdiff(k)^2 * abs(sin(thetan))*(cos(phidiff(k-1))-cos(phidiff(k))))/rn;
            if (isnan(factor(k)))
                rn
            end
         end         
   end
    
    luminosity(1,t)=dot(flux1,factor)
    
%     index_k=find(abs(factor)>prctile(abs(factor),92));
%     factor(index_k)=[];
%     flux1(index_k)=[];
%     luminosity(t)=dot(flux1,factor);
%      clear xdiff
%      clear ydiff
%      clear Y_new
%      clear X_new
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
    %scatter(xdiff(index_k)*2/rconv,ydiff(index_k)*2/rconv,3,'g'); hold on
    clear index_k
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
    axis equal
    axis([0 4 0 4])
     name=['/home/nilou/Data/plot/BSGdiff/BSG_filter_new_Phi0' num2str(t) '.png'];
     print(gcf, '-dpng', '-r50', name)
     export_fig(name, '-dpng', '-r50')
end
save('/home/nilou/Data/processeddata/BSG/luminosity_new_Phi90.mat','luminosity')
