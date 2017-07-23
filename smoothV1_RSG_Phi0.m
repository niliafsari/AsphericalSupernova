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

% luminosity=zeros(1,50);
% luminosity90=zeros(1,50);
% luminosity_tot=zeros(1,50);

%load('luminosity.mat', 'luminosity');
% temp=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55.mat','luminosity');
% luminosity=temp.luminosity;
time=load('/home/nilou/Data/timesteps.mat');

for t=50:50
    t
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
    
    rdiff=rdiff(phidiff<=prctile(phidiff,100.0));
    xdiff=xdiff(phidiff<=prctile(phidiff,100.0));
    ydiff=ydiff(phidiff<=prctile(phidiff,100.0));
    phidiff=atan(xdiff./ydiff);
    
    phidiff=phidiff(rdiff>(0.5*rconv));
    xdiff=xdiff(rdiff>(0.5*rconv));
    ydiff=ydiff(rdiff>(0.5*rconv));
    rdiff=rdiff(rdiff>(0.5*rconv));
    
    [phidiff,I]=sort(phidiff);    
    rdiff=rdiff(I);
    xdiff=xdiff(I);
    ydiff=ydiff(I);
    
    
    x=1:length(ydiff);
    x=x';
    fit1 = fit(x,ydiff','poly6');
    fdata = feval(fit1,x);
    I = abs(fdata - ydiff') > 0.2*std(ydiff');
    outliers = excludedata(x,ydiff','indices',I);
    %sum(outliers)
 
    rdiff(outliers)=[];
    xdiff(outliers)=[];
    ydiff(outliers)=[];
    phidiff=atan(xdiff./ydiff);
    
    x=1:length(rdiff);
    x=x';
    fit1 = fit(x,rdiff','poly5');
    fdata = feval(fit1,x);
    I = abs(fdata - rdiff') > 0.2*std(rdiff');
    outliers = excludedata(x,rdiff','indices',I);
    %sum(outliers)
 
    rdiff(outliers)=[];
    xdiff(outliers)=[];
    ydiff(outliers)=[];
    phidiff=atan(xdiff./ydiff);
        
%     x=1:length(phidiff);
%     x=x';
%     fit1 = fit(x,phidiff','poly6');
%     fdata = feval(fit1,x);
%     I = abs(fdata - phidiff') > 0.4*std(phidiff');
%     outliers = excludedata(x,phidiff','indices',I);
    %sum(outliers)
 
%     rdiff(outliers)=[];
%     xdiff(outliers)=[];
%     ydiff(outliers)=[];
%     phidiff=atan(xdiff./ydiff);
    
    [xmax,I]=max(xdiff);
      
    P=[xdiff', ydiff'];

    DT = delaunayTriangulation(P);
    Q = convexHull(DT);
    Q
    if (t==26)
        Q(1:3)=[];
    elseif (t==27)
       Q(1:2)=[]; 
    elseif (t==23)
       Q(1:3)=[];    
    else
        Q(1)=[];
    end
    Q=flipud(Q);
    rdiff=rdiff(Q);
    xdiff=xdiff(Q);
    ydiff=ydiff(Q);
    phidiff=phidiff(Q);


%       rdiff_k=rdiff(phidiff>phidiff(I));
%       xdiff_k=xdiff(phidiff>phidiff(I));
%       ydiff_k=ydiff(phidiff>phidiff(I));
%       phidiff_k=phidiff(phidiff>phidiff(I));
%       rdiff=rdiff(phidiff<=phidiff(I));
%       xdiff=xdiff(phidiff<=phidiff(I));
%       ydiff=ydiff(phidiff<=phidiff(I));
%       phidiff=phidiff(phidiff<=phidiff(I));
%       
%         P=[xdiff_k', ydiff_k'];
%         if (length(P)>0)
%             DT = delaunayTriangulation(P);
%             Q = convexHull(DT);
%             Q(1)=[];
%             Q=flipud(Q);
%             rdiff_k=rdiff_k(Q);
%             xdiff_k=xdiff_k(Q);
%             ydiff_k=ydiff_k(Q);
%             phidiff_k=phidiff_k(Q);
%         else
%             rdiff_k=[];
%             xdiff_k=[];
%             ydiff_k=[];
%             phidiff_k=[];
%         end
%         
%     rdiff=[rdiff rdiff_k];
%     xdiff=[xdiff xdiff_k];
%     ydiff=[ydiff ydiff_k];
%     phidiff=[phidiff phidiff_k];
    [xmax,I]=max(xdiff); 
    [ymax,U]=max(ydiff); 
        
    
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
    factor90=zeros(1,length(index_r));
    factor_tot=zeros(1,length(index_r));
    densityk=zeros(1,length(index_r));
    gradpk=zeros(1,length(index_r));
    
    for k=2:length(index_r)
         densityk(k)=density(index_r(k),index_phi(k));
         gradpk(k)=gradp(index_r(k),index_phi(k));
         flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
         rdiffm=(rdiff(k)+rdiff(k-1))/2;
         dl=sqrt((rdiffm*(phidiff(k)-phidiff(k-1)))^2+(rdiff(k)-rdiff(k-1))^2);
         rn=rdiffm*(phidiff(k)-phidiff(k-1))/dl;
         if (rn==0)
             factor_tot(k)=0;
         else
            factor_tot(k)=(2*pi* rdiff(k)^2 * (cos(phidiff(k-1))-cos(phidiff(k))))/rn;         
         end
         dz=ydiff(k)-ydiff(k-1);
         dR=xdiff(k)-xdiff(k-1);
         thetan=atan(dz/dR);
         
        % if ((dR<0 && dz<0) || (dR<0 && dz>0) || rn==0 || (phidiff(k)>phidiff(I)))
         if ((dR<0 && dz<0) || (dR<0 && dz>0) || rn==0 || (phidiff(k)>phidiff(I)))
             factor(k)=0;
             rn1(k)=rn;
             thetan1(k)=thetan;
         else    
            thetan=atan(abs(dz)/dR);
            rn1(k)=rn;
            thetan1(k)=thetan;
            factor(k)=(8*pi* rdiff(k)^2 * cos(thetan)*(cos(phidiff(k-1))-cos(phidiff(k))))/rn;
         end
         
         
         if ((dR>0 && dz>0) || (dR<0 && dz>0) || rn==0 || (phidiff(k)<phidiff(U)))
             factor90(k)=0;
         else
            if  ((dR<0 && dz<0))
                thetan=(pi/2-atan(dz/dR))+pi/2;
            else
                thetan=atan(-dz/dR);
            end
            factor90(k)=(8*rdiff(k)^2 * sin(thetan)*(cos(phidiff(k-1))-cos(phidiff(k))))/rn;
         end
    end
    
    luminosity(1,t)=dot(flux1,factor)
    luminosity90(1,t)=dot(flux1,factor90)
    luminosity_tot(1,t)=dot(flux1,factor_tot)
    
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
       
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height]) 
    scatter(P(:,1)*2/rconv,P(:,2)*2/rconv,1,'b');  
    hold on
    plot(xdiff*2/rconv,ydiff*2/rconv,'r');
    
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
    axis equal
    axis([0 4 0 4])
    name=['/home/nilou/Data/plot/RSGdiff/RSG_k_' num2str(t) '.png'];
    print(gcf, '-dpng', '-r50', name)
    export_fig(name, '-dpng', '-r50')

end
save('/home/nilou/Data/processeddata/RSG/luminosity_0.mat','luminosity')
save('/home/nilou/Data/processeddata/RSG/luminosity_90.mat','luminosity90')
save('/home/nilou/Data/processeddata/RSG/luminosity_tot.mat','luminosity_tot')
