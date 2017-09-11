mtot=4.82611e-3;
etot=1.3402e-2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=14*msun;
e=1e51;
r=400*rsun;
kappa=0.34;

c=3e10;
h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;
Z=0.005;
a=7.566e-15;
X=0.7;

rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/(rconv^3);
rhoconv=mconv/(rconv^3);
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

r=linspace(0,2,2048).*rconv;
theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048)*rconv;
for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end

cita=1;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end
time=load([path '/processeddata/timesteps_1024.mat']);
load([path '/f.mat'],'f')
time_rsg=time.time1*tconv;
C=3.4e36;
prefac=C*2.7*k_B*rhoconv^2/(3^(7/8)*a^(1/8)*pconv^(7/8))
for t=172:420
    t
    diff_bsg1=load([path '/processeddata/RSG/diff_rsg_1024_' num2str(t-1) '.mat'], 'diff_bsg'); 
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
    name=[path '/rawdata/gfunction/gfun1024_' int2str(t-1) '.mat'] ;
    load(name,'gfun_data');
    gfun=gfun_data*prefac;
    
    t_loc_all=(3*pres/a).^0.25;
    tdiff_all=(3*kappa*density).*(pres.^2)./(c*gradp.^2);
    %etha_all=((density.^-2).*(t_loc_all.^3.5))./(3.7e22*(1+X)*(1-Z)*c*min(time_rsg(t),tdiff_all));
    etha_all=(7e5./min(time_rsg(t),tdiff_all)).*((density/1e-10).^-2).*(t_loc_all/(100*11604.52)).^3.5;

    
    
    diff_bsg_t=diff_bsg1.diff_bsg(:,2:length(diff_bsg1.diff_bsg));

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
    
    if t>172
        phidiff=phidiff(rdiff>(0.5*rconv));
        xdiff=xdiff(rdiff>(0.5*rconv));
        ydiff=ydiff(rdiff>(0.5*rconv));
        rdiff=rdiff(rdiff>(0.5*rconv));
    end     
    
    [phidiff,I]=sort(phidiff);    
    rdiff=rdiff(I);
    xdiff=xdiff(I);
    ydiff=ydiff(I);
    
    
if (t>178)    
    x=1:length(ydiff);
    x=x';
    fit1 = fit(x,ydiff','poly6');
    fdata = feval(fit1,x);
    I = abs(fdata - ydiff') > 0.15*std(ydiff');
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
    I = abs(fdata - rdiff') > 0.15*std(rdiff');
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
    Q(1)=[]
%     if (t==26)
%         Q(1:3)=[];
%     elseif (t==27)
%        Q(1:2)=[]; 
%     elseif (t==23)
%        Q(1:3)=[];    
%     else
%         Q(1)=[];
%     end
    Q=flipud(Q);
    rdiff=rdiff(Q);
    xdiff=xdiff(Q);
    ydiff=ydiff(Q);
    phidiff=phidiff(Q);
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
    gfunction=zeros(1,length(index_r));
    flux1=zeros(1,length(index_r));
    factor=zeros(1,length(index_r));
    factor90=zeros(1,length(index_r));
    factor_tot=zeros(1,length(index_r));
    densityk=zeros(1,length(index_r));
    gradpk=zeros(1,length(index_r));
    dL=zeros(1,length(index_r));
    tempBB=zeros(1,length(index_r));
    etha=zeros(1,length(index_r));
    t_c=zeros(1,length(index_r));
    for k=2:length(index_r)
         gfunction(k)=gfun(index_r(k),index_phi(k));
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
         dL(k)=factor_tot(k)*flux1(k);
         tempBB(k)=t_loc_all(index_r(k),index_phi(k));
         etha(k)=etha_all(index_r(k),index_phi(k));
         t_c(k)=((10.^f(gfunction(k)))+1)*tempBB(k)/(1+(1/(etha(k)^2)))^(1/17);
    end
    
    t_tot(1,t)= dot(dL,t_c)/dot(flux1,factor_tot);
    
      %  close all
%     a=get(gcf,'Position');
%     x0=15;
%     y0=15;
%     width=650;
%     height=650;
%        
%     myFigure = figure('PaperPositionMode','auto','Position',a);
%       
%     set(myFigure,'units','points','position',[x0,y0,width,height]) 
%     h=surf(xx*2/rconv,yy*2/rconv,log10(etha_all));hold on
%     grid off
%     set(h,'LineStyle','none');
%     colormap jet
%     cl=colorbar;
%     scatter(P(:,1)*2/rconv,P(:,2)*2/rconv,1,'b');  
%     hold on
%     %plot(xdiff*2/rconv,ydiff*2/rconv,'r');
%     view(2)
%     xlabel('x/R_*');
%     ylabel('y/R_*');
%     title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
%     axis equal
%     axis([0 4 0 4])
%     name=[path '/plot/RSGdiff/RSG_1024_' num2str(t) '.png'];
%     print(gcf, '-dpng', name)
%     export_fig(name, '-dpng')
end
%save([path '/processeddata/RSG/colortem_1024_tot.mat'],'t_tot')
