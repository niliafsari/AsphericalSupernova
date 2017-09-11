mtot=4.82611e-3*2;
etot=1.3402e-2*2;
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

cita=1;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end

 %luminosity=zeros(1,370);
 %luminosity90=zeros(1,370);
 %luminosity_tot=zeros(1,370);

load([path '/processeddata/BSG/luminosity_1024_0.mat'],'luminosity')
load([path '/processeddata/BSG/luminosity_1024_90.mat'],'luminosity90')
load([path '/processeddata/BSG/luminosity_1024_tot.mat'],'luminosity_tot')

%load('luminosity.mat', 'luminosity');
% temp=load('/home/nilou/Data/processeddata/RSG/luminosity_smoothed55.mat','luminosity');
% luminosity=temp.luminosity;
time=load([path '/processeddata/timesteps_1024.mat']);

for t=165:174
    t
    load([path '/processeddata/BSG/diff_bsg_1024_' num2str(t-1) '.mat'], 'diff_bsg'); 
    name=[ path '/rawdata/gradp2/gradp21024_' int2str(t-1) '.mat'];
    load(name,'gradp2'); 
    gradp=sqrt(gradp2*pconv^2/rconv^2);
    if (t<=200 && t~=176)
        name=[path '/rawdata/pressure/pres1024_' int2str(t-1) '.csv'] ;
        pres= csvread(name).*pconv;
        name=[path '/rawdata/density/dens1024_' int2str(t-1) '.csv'] ;
        density= csvread(name).*rhoconv;       
    else
        name=[path '/rawdata/pressure/pres1024_' int2str(t-1) '.mat'] ;
        load(name,'pres_data');
        pres=pres_data*pconv;
        name=[path '/rawdata/density/dens1024_' int2str(t-1) '.mat'] ;
        load(name,'dens_data');
        density=dens_data*rhoconv;
    end


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
    rdiff_n=rdiff_n(phidiff_n<=prctile(phidiff_n,99.5));
    xdiff_n=xdiff_n(phidiff_n<=prctile(phidiff_n,99.5));
    ydiff_n=ydiff_n(phidiff_n<=prctile(phidiff_n,99.5));
    phidiff_n=atan(xdiff_n./ydiff_n);
    
    if t>170
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
    if t>174
        x=1:length(rdiff);
        x=x';
        fit1 = fit(x,rdiff','poly6');
        fdata = feval(fit1,x);
        I = abs(fdata - rdiff') > 0.4*std(rdiff');
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
                Q(1)=[];
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
    thetan2=zeros(1,length(index_r));
    factor_tot=zeros(1,length(index_r));
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
         if (rn==0)
             factor_tot(k)=0;
         else
            factor_tot(k)=(2*pi* rdiff(k)^2 * (cos(phidiff(k-1))-cos(phidiff(k))))/rn;         
         end
         dz=ydiff(k)-ydiff(k-1);
         dR=xdiff(k)-xdiff(k-1);
         thetan=atan(dz/dR);
         
        % if ((dR<0 && dz<0) || (dR<0 && dz>0) || rn==0 || (phidiff(k)>phidiff(I)))
         if ((dR<0 && dz<0) || (dR<0 && dz>0) || rn==0 || (phidiff(k)>phidiff(I)) )
             factor(k)=0;
             rn1(k)=rn;
             thetan1(k)=thetan;
         else    
            thetan=atan(abs(dz)/dR);
            rn1(k)=rn;
            thetan1(k)=thetan;
            factor(k)=(8*pi* rdiff(k)^2 * cos(thetan)*(cos(phidiff(k-1))-cos(phidiff(k))))/rn;
         end
         
         
         if ((dR>=0 && dz>=0) || (dR<=0 && dz>=0) || rn==0 || (phidiff(k)<phidiff(U)))
             factor90(k)=0;
         else
            if  ((dR<0 && dz<=0))
                thetan=(pi/2-atan(dz/dR))+pi/2;
            else
                thetan=atan(-dz/dR);
            end
            thetan2(k)=thetan;
            if (thetan2(k) ==0 && phidiff(k)>0.9)
                dz
                dR
            end
            factor90(k)=(8*rdiff(k)^2 * sin(thetan)*(cos(phidiff(k-1))-cos(phidiff(k))))/rn;
         end
    end
    
    luminosity(1,t)=dot(flux1,factor)
    luminosity90(1,t)=dot(flux1,factor90)
    luminosity_tot(1,t)=dot(flux1,factor_tot)
    
    %close all
%     a=get(gcf,'Position');
%     x0=15;
%     y0=15;
%     width=650;
%     height=650;
%        
%     myFigure = figure('PaperPositionMode','auto','Position',a);
%       
%     set(myFigure,'units','points','position',[x0,y0,width,height]) 
%     h=surf(xx*2/rconv,yy*2/rconv,log10(density));hold on
%     grid off
%     set(h,'LineStyle','none');
%     colormap jet
%     cl=colorbar;
%     scatter(P(:,1)*2/rconv,P(:,2)*2/rconv,1,'b');  
%     hold on
    plot(xdiff_t*2/rconv,ydiff_t*2/rconv,'r');
    view(2)
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
    axis equal
    view(2)
   axis([0 4 0 4])
%     name=[path '/plot/RSGdiff/BSG_1024_' num2str(t) '.png'];
%     print(gcf, '-dpng', name)
%     export_fig(name, '-dpng')

end
save([path '/processeddata/BSG/luminosity_1024_0.mat'],'luminosity')
save([path '/processeddata/BSG/luminosity_1024_90.mat'],'luminosity90')
save([path '/processeddata/BSG/luminosity_1024_tot.mat'],'luminosity_tot')
