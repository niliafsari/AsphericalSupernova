mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=15*msun;
e=1e51;
r=49*rsun;

Z=0.005;
a=7.566e-15;
X=0.7;
kappa_T=0.2*(1+X);
kappa=kappa_T;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

r=0.5*rconv;
theta=linspace(0,pi/2,2048);

xxx=r*sin(theta);
y=r*cos(theta);

radius=linspace(0,2,2048)*rconv;

for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
time=load('/home/nilou/Data/timesteps_1024.mat');
%load('/home/nilou/Data/processeddata/BSG/eta_com_bsg.mat','eta_com');
time_bsg=(time.time1*tconv);

cita=0;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end

C=3.4e36;
prefac=C*2.7*k_B*rhoconv^2/(3^(7/8)*a^(1/8)*pconv^(7/8));
for t=230:230
    a=7.566e-15;
    t
    radius=zeros(1,2048);

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
    etha_all=(7e5./min(time_ic(t),tdiff_all)).*((density/1e-10).^-2).*(t_loc_all/(100*11604.52)).^3.5;
%    etha_all=((density.^-2).*(t_loc_all.^3.5))./(3.7e22*(1+X)*(1-Z)*c*min(time_bsg(t),tdiff_all));
%    etha_all=eta_com(:,:,t-22);
%    dsmooth=eta_com(:,:,t-22);
    ke=30;
    if (t>48)
        ke=75;
    end
    temp=zeros(2048+2*ke,2048+ke);
    temp(1:2048,1:2048)=etha_all;
    kernel = 1*fspecial('disk', ke);
    temp1 = imfilter(log10(real(temp)),kernel,'same');
    dsmooth= temp1(1:2048,1:2048);   
    [therm_bsg,~]=contour3(xx,yy,dsmooth,[1 1],'LineColor','b','Linewidth',0.3);  
     rtherm=sqrt(therm_bsg(1,1:length(therm_bsg)).^2+therm_bsg(2,1:length(therm_bsg)).^2);
     xtherm=therm_bsg(1,1:length(therm_bsg));
     ytherm=therm_bsg(2,1:length(therm_bsg)); 
     phitherm=atan(xtherm./ytherm);
    dsmooth=10.^dsmooth;
    [phitherm,I]=sort(phitherm);
    rtherm=rtherm(I);
    xtherm=xtherm(I);
    ytherm=ytherm(I);  
    
    
    phitherm=phitherm(rtherm>(0.5*rconv));
    xtherm=xtherm(rtherm>(0.5*rconv));
    ytherm=ytherm(rtherm>(0.5*rconv));
    rtherm=rtherm(rtherm>(0.5*rconv)); 
    
    if t>24
        x=1:length(rtherm);
        x=x';
        fit1 = fit(x,rtherm','poly5');
        fdata = feval(fit1,x);
        I = abs(fdata-rtherm') > 0.4*std(rtherm');
        outliers = excludedata(x,rtherm','indices',I);
        sum(outliers)

        rtherm(outliers)=[];
        xtherm(outliers)=[];
        ytherm(outliers)=[];
        phitherm=atan(xtherm./ytherm);
    end
    
     index_therm_r=floor(rtherm/((2/2048)*rconv));
     index_therm_phi=floor(atan(xtherm./ytherm)/(0.5*pi/2048));
     index_therm_phi(index_therm_phi==0)=1;
     index_therm_r(index_therm_r==0)=1;
    
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
    
    rdiff_n=rdiff_n(phidiff_n<=prctile(phidiff_n,99));
    xdiff_n=xdiff_n(phidiff_n<=prctile(phidiff_n,99));
    ydiff_n=ydiff_n(phidiff_n<=prctile(phidiff_n,99));
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
    factor_phi0=zeros(1,length(index_r));
    t_loc=zeros(1,length(index_r));
    etha=zeros(1,length(index_r));
    densityk=zeros(1,length(index_r));
    gradpk=zeros(1,length(index_r));
    luminosity_tot(1,t)=0;
    t_color_tot(1,t)=0;
    luminosity_phi0(1,t)=0;
    t_color_phi0(1,t)=0;
    

    for k=2:length(index_r)
         densityk(k)=density(index_r(k),index_phi(k));
         gradpk(k)=gradp(index_r(k),index_phi(k));
         flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
         factor(k)=2*pi* rdiff(k)^2 * (cos(phidiff(k-1))-cos(phidiff(k)));         
         factor_phi0(k)=4*pi* rdiff(k)^2 * (sin(phidiff(k))^2-sin(phidiff(k-1))^2);

         tmp = abs(phitherm-phidiff(k));
         [~,k_therm] = min(tmp); 
         if (dsmooth(index_r(k),index_phi(k))<=1)
            t_loc(k)=(3*pres(index_therm_r(k_therm),index_therm_phi(k_therm))/a)^0.25;       
            etha(k)=1;
         else
            if (isnan(dsmooth(index_r(k),index_phi(k))))
                continue
            end
            %t_loc(k)=(3*pres(index_therm_r(k_therm),index_therm_phi(k_therm))/a)^0.25;
            t_loc(k)=(3*pres(index_r(k),index_phi(k))/a)^0.25;
            etha(k)=dsmooth(index_r(k),index_phi(k));
         end         
    end
    luminosity_tot(1,t)=dot(flux1,factor);
    t_color_tot(1,t)= sum((flux1.*factor.*t_loc.*etha.^2)/luminosity_tot(1,t))
    luminosity_phi0(1,t)=dot(flux1,factor_phi0);
    t_color_phi0(1,t)= sum((flux1.*factor_phi0.*t_loc.*etha.^2)/luminosity_phi0(1,t));
    
    tau_eff=sqrt(3.*tau.*tauff);
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
    h1=contour3(xx,yy,tau,[1 1],'LineColor','m','Linewidth',1.5);hold on
 %   axis equal
    scatter(xdiff,ydiff,0.8,'g');hold on
    
    %h2=contour(xx,yy,d,[1 1],'LineColor','g');
    %h1=contour3(xx,yy,etha_all,[1 1],'LineColor','k','Linewidth',1);hold on 
    contour3(xx,yy,dsmooth,[1 1],'LineColor','k','Linewidth',1.5);  
    h=surf(xx,yy, log10(tau));hold on
    grid off
    set(h,'LineStyle','none');
    colormap jet
    alpha(.2)
    c=colorbar
    title(c,['Log \tau']);
    caxis auto
    view(2);
    axis([0,2*rconv,0,2*rconv])
    title(['BSG t=' num2str(time.time1(t)*tconv,3) '[s]']);
    xlabel('cm'); 
    ylabel('cm');
    legend('R/R_*=1','\tau =1 ','Diffusion Front','\eta =1','Location','best');
    %set(gca,'LineWidth',2,'FontSize',12);
    %set(gca,'XDir','reverse')
    %name=['/home/nilou/Data/plot/BSGdiff/colortempContourWdiff_v2_' int2str(t)];
    %print('-dpng',name) 
    %export_fig(name, '-png')
    set(gca,'LineWidth',2,'FontSize',12);
    name=['/home/nilou/Data/rth_BSG.pdf'];
    print('-dpdf',name) 
    export_fig(name, '-pdf')
end
% save('/home/nilou/Data/processeddata/BSG/colortemp_tot_v3.mat','t_color_tot')
% save('/home/nilou/Data/processeddata/BSG/colortemp_phi0_v3.mat','t_color_phi0')