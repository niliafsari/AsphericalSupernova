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
h=6.62e-27;
k_B=1.38e-16;
sigma=5.6e-5;
Z=0.005;
a=7.566e-15;
X=0.7;

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
r=0.5*rconv;
xxx=r*sin(theta);
y=r*cos(theta);

gradp=zeros(2048,2048);
cita=0;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end
time=load([path '/processeddata/timesteps_1024.mat']);
load([path '/f.mat'],'f')
time_rsg=time.time1*tconv;
C=3.4e36;
prefac=C*2.7*k_B*rhoconv^2/(3^(7/8)*a^(1/8)*pconv^(7/8));
for t=211:211
    t
    load([path '/processeddata/BSG/diff_bsg_1024_' num2str(t-1) '.mat'], 'diff_bsg'); 
    name=[ path '/rawdata/gradp2/gradp21024_' int2str(t-1) '.mat'];
    load(name,'gradp2'); 
    gradp=sqrt(gradp2*pconv^2/rconv^2);
    
    name=[path '/rawdata/gfunction/gfun1024_' int2str(t-1) '.mat'] ;
    load(name,'gfun_data');
    gfun=gfun_data*prefac;
    
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


    t_loc_all=(3*pres/a).^0.25;
    tdiff_all=(3*kappa*density).*(pres.^2)./(c*gradp.^2);
    %etha_all=((density.^-2).*(t_loc_all.^3.5))./(3.7e22*(1+X)*(1-Z)*c*min(time_rsg(t),tdiff_all));
    etha_all=(7e5./min(time_rsg(t),tdiff_all)).*((density/1e-10).^-2).*(t_loc_all/(100*11604.52)).^3.5;

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
    ke=15;
    temp=zeros(2048+2*ke,2048+ke);
    temp(1:2048,1:2048)=etha_all;
    kernel = 1*fspecial('disk', ke);
    temp1 = imfilter(log10(real(temp)),kernel,'same');
    dsmooth= temp1(1:2048,1:2048); 
    
        close;
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=350;
    height=300;
    myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
    set(myFigure,'units','points','position',[x0,y0,width,height])
        
    p3=plot(xxx,y,'r','Linewidth',1.5); hold on
    %h1=contour3(xx,yy,tau_eff,[2/3 2/3],'LineColor','b','Linewidth',1);hold on
    %h1=contour3(xx,yy,tau,[1 1],'LineColor','m','Linewidth',1.5);hold on
 %   axis equal
    scatter(xdiff_n,ydiff_n,0.8,'g');hold on
    
    %h2=contour(xx,yy,d,[1 1],'LineColor','g');
    %h1=contour3(xx,yy,etha_all,[1 1],'LineColor','k','Linewidth',1);hold on 
    contour3(xx,yy,dsmooth,[1 1],'LineColor','k','Linewidth',1.5);  
    h=surf(xx,yy, log10(density));hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    alpha(.2)
    c=colorbar
    title(c,['Log (\rho [g/cm^{-3}])']);
    caxis([-14 -1])
    view(2);
    axis([0,2*rconv,0,2*rconv])
    title(['BSG t=' num2str(time.time1(t)*tconv,3) '[s]']);
    xlabel('cm'); 
    ylabel('cm');
    legend('R/R_*=1','Diffusion Front','\eta =1','Location','best');
    %set(gca,'LineWidth',2,'FontSize',12);
    %set(gca,'XDir','reverse')
    %name=['/home/nilou/Data/plot/BSGdiff/colortempContourWdiff_v2_' int2str(t)];
    %print('-dpng',name) 
    %export_fig(name, '-png')
    set(gca,'LineWidth',2,'FontSize',12);
    name=['/home/nilou/Data/rth_BSG_1024.pdf'];
    print('-dpdf',name) 
    export_fig(name, '-pdf')
%        close all
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
%     contour(xx*2/rconv,yy*2/rconv,log10(etha_all),[0 0],'-r');hold on
%     %scatter(P(:,1)*2/rconv,P(:,2)*2/rconv,1,'b');  
%     hold on
%     scatter(xdiff*2/rconv,ydiff*2/rconv,2,'m');
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
%save([path '/processeddata/BSG/colortem_1024_tot.mat'],'t_tot')
