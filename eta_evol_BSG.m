mtot=0.0048187313;
etot=0.0130949665;
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
% load('/home/nilou/Data/processeddata/BSG/eta.mat','eta');
% load('/home/nilou/Data/processeddata/BSG/m_coord.mat','m_coord');
eta_com=zeros(2048,2048,42-22);
time=load('/home/nilou/Data/timesteps.mat');
time_bsg=(time.time1*tconv);
for t=23:42
    a=7.566e-15;
    t   
    if (t==23)
         eta_com(:,:,t-22)=eta(:,:,t-22);
         continue       
    else
        first=reshape(m_coord(:,:,t-1-22),[1 2048*2048]);
        second=reshape(m_coord(:,:,t-22),[1 2048*2048]);
        row=zeros(1,2048*2048);
        col=zeros(1,2048*2048);
        for i=1:2048
            id=knnsearch(first((i-1)*2048+1:i*2048)',second((i-1)*2048+1:i*2048)');
            col((i-1)*2048+1:i*2048)=i;
            row((i-1)*2048+1:i*2048)=id;
        end
        idx=sub2ind([2048 2048], row,col);
        first_eta=eta(:,:,t-22);
        old_eta=reshape(eta_com(:,:,t-1-22),[1 2048^2]);
    end
    old_eta=reshape(old_eta(idx),[2048 2048]);
    sum(sum(first_eta>old_eta))
    eta_com(:,:,t-22)=min(first_eta, old_eta);
    
    
%     
%     name1=['/home/nilou/Data/rawdata/dcodeunit/dcodeunit_' int2str(t-1) '.csv'] ;
%     d=csvread(name1);
%     indexes=find(isnan(d) | isinf(d));
%     indexes_tot=indexes(((xx(indexes).^2+yy(indexes).^2)/rconv)>0.5^2);
%     name=['/home/nilou/Data/processeddata/BSG/tauBSGradial_' int2str(t-1) '.csv'] ;
%     tau= csvread(name);
%     nameff=['/home/nilou/Data/processeddata/BSG/tauBSGradialFF_' int2str(t-1) '.csv'] ;
%     tauff= csvread(nameff);
%     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
%     density= csvread(name)*rhoconv;
%     density(indexes_tot)=0;
%     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
%     pres= csvread(name)*pconv;
%     name=['/home/nilou/Data/processeddata/BSG/dparamBSG_v2_' int2str(t-1) '.csv'] ;
%     d= csvread(name);
%     name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-2) '.csv'];
%     gradp=csvread(name);
%     diff_bsg1=load(['/home/nilou/Data/processeddata/BSG/diff_bsg_smoothed_k55_temp' num2str(t) '.mat']);
%      
%     t_loc_all=(3*pres/a).^0.25;
%     tdiff_all=(3*kappa*density).*(pres.^2)./(c*gradp.^2);
%     etha_all=((density.^-2).*(t_loc_all.^3.5))./(3.7e22*(1+X)*(1-Z)*c*min(time_bsg(t),tdiff_all));
%     
%     dm=4*pi*(xx.^2+yy.^2).*density;
%     for i=1:2048
%         m_coord(i,:,t-22)=trapz(dm(i:2048,:),1); 
%     end
%     
%     
%     ke=30;
%     if (t>48)
%         ke=75;
%     end
%     temp=zeros(2048+2*ke,2048+ke);
%     temp(1:2048,1:2048)=etha_all;
%     kernel = 1*fspecial('disk', ke);
%     temp1 = imfilter(log10(real(temp)),kernel,'same');
%     dsmooth= temp1(1:2048,1:2048);   
%     [therm_bsg,~]=contour3(xx,yy,dsmooth,[0 0],'LineColor','b','Linewidth',0.3);  
%      rtherm=sqrt(therm_bsg(1,1:length(therm_bsg)).^2+therm_bsg(2,1:length(therm_bsg)).^2);
%      xtherm=therm_bsg(1,1:length(therm_bsg));
%      ytherm=therm_bsg(2,1:length(therm_bsg)); 
%      phitherm=atan(xtherm./ytherm);
%     dsmooth=10.^dsmooth;
%     eta(:,:,t-22)=dsmooth;
%     
% %     
% %     [phitherm,I]=sort(phitherm);
% %     rtherm=rtherm(I);
% %     xtherm=xtherm(I);
% %     ytherm=ytherm(I);  
% %     
% %     
% %     phitherm=phitherm(rtherm>(0.5*rconv));
% %     xtherm=xtherm(rtherm>(0.5*rconv));
% %     ytherm=ytherm(rtherm>(0.5*rconv));
% %     rtherm=rtherm(rtherm>(0.5*rconv)); 
% %     
% %     if t>23
% %         x=1:length(rtherm);
% %         x=x';
% %         fit1 = fit(x,rtherm','poly5');
% %         fdata = feval(fit1,x);
% %         I = abs(fdata-rtherm') > 0.4*std(rtherm');
% %         outliers = excludedata(x,rtherm','indices',I);
% %         sum(outliers)
% % 
% %         rtherm(outliers)=[];
% %         xtherm(outliers)=[];
% %         ytherm(outliers)=[];
% %         phitherm=atan(xtherm./ytherm);
% %     end
% %     
% %      index_therm_r=floor(rtherm/((2/2048)*rconv));
% %      index_therm_phi=floor(atan(xtherm./ytherm)/(0.5*pi/2048));
% %      index_therm_phi(index_therm_phi==0)=1;
% %      index_therm_r(index_therm_r==0)=1;
% %     
% %     if (t==23)
% %          diff_bsg_t=diff_bsg1.diff_bsg(:,3:length(diff_bsg1.diff_bsg));
% %      else
% %          diff_bsg_t=diff_bsg1.diff_bsg(:,2:length(diff_bsg1.diff_bsg));
% %      end
% %     diff_bsg_t=transpose(diff_bsg_t);
% %     %length(diff_bsg_t)
% %     diff_bsg=unique(diff_bsg_t,'rows');
% %     length(diff_bsg)
% %     diff_bsg=transpose(diff_bsg);
% %     rdiff_n=sqrt(diff_bsg(1,1:length(diff_bsg)).^2+diff_bsg(2,1:length(diff_bsg)).^2);
% %     xdiff_n=diff_bsg(1,1:length(diff_bsg));
% %     ydiff_n=diff_bsg(2,1:length(diff_bsg)); 
% %     phidiff_n=atan(xdiff_n./ydiff_n);
% %     
% %     rdiff_n=rdiff_n(phidiff_n<=prctile(phidiff_n,99));
% %     xdiff_n=xdiff_n(phidiff_n<=prctile(phidiff_n,99));
% %     ydiff_n=ydiff_n(phidiff_n<=prctile(phidiff_n,99));
% %     phidiff_n=atan(xdiff_n./ydiff_n);
% %     
% %     rdiff_t=rdiff_n(rdiff_n>(0.5*rconv));
% %     phidiff_t=phidiff_n(rdiff_n>(0.5*rconv));
% %     xdiff_t=xdiff_n(rdiff_n>(0.5*rconv));
% %     ydiff_t=ydiff_n(rdiff_n>(0.5*rconv));
% %     
% %     
% %     
% %     [phidiff,I]=sort(phidiff_t);
% %     rdiff=rdiff_t(I);
% %     xdiff=xdiff_t(I);
% %     ydiff=ydiff_t(I);    
% % 
% %     
% %     x=1:length(rdiff);
% %     x=x';
% %     fit1 = fit(x,rdiff','poly5');
% %     fdata = feval(fit1,x);
% %     I = abs(fdata - rdiff') > 0.6*std(rdiff');
% %     outliers = excludedata(x,rdiff','indices',I);
% %     sum(outliers)
% %  
% %     rdiff(outliers)=[];
% %     xdiff(outliers)=[];
% %     ydiff(outliers)=[];
% %     phidiff=atan(xdiff./ydiff);
% %     
% % 
% %     index_r=floor(rdiff/((2/2048)*rconv));
% %     index_phi=floor(atan(xdiff./ydiff)/(0.5*pi/2048));
% %     index_phi(index_phi==0)=1;
% %     index_r(index_r==0)=1;
% %     if  sum(index_r>2048)
% %         sum(index_r>2048)
% %     end
% %     index_r(index_r>2048)=2048;
% %     flux1=zeros(1,length(index_r));
% %     factor=zeros(1,length(index_r));
% %     factor_phi0=zeros(1,length(index_r));
% %     t_loc=zeros(1,length(index_r));
% %     etha=zeros(1,length(index_r));
% %     densityk=zeros(1,length(index_r));
% %     gradpk=zeros(1,length(index_r));
% %     luminosity_tot(1,t)=0;
% %     t_color_tot(1,t)=0;
% %     luminosity_phi0(1,t)=0;
% %     t_color_phi0(1,t)=0;
% %     
% % 
% %     for k=2:length(index_r)
% %          densityk(k)=density(index_r(k),index_phi(k));
% %          gradpk(k)=gradp(index_r(k),index_phi(k));
% %          flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
% %          factor(k)=2*pi* rdiff(k)^2 * (cos(phidiff(k-1))-cos(phidiff(k)));         
% %          factor_phi0(k)=4*pi* rdiff(k)^2 * (sin(phidiff(k))^2-sin(phidiff(k-1))^2);
% % 
% %          tmp = abs(phitherm-phidiff(k));
% %          [~,k_therm] = min(tmp); 
% %          if (dsmooth(index_r(k),index_phi(k))<=1)
% %             t_loc(k)=(3*pres(index_therm_r(k_therm),index_therm_phi(k_therm))/a)^0.25;       
% %             etha(k)=1;
% %          else
% %             %t_loc(k)=(3*pres(index_therm_r(k_therm),index_therm_phi(k_therm))/a)^0.25;
% %             t_loc(k)=(3*pres(index_r(k),index_phi(k))/a)^0.25;
% %             etha(k)=dsmooth(index_r(k),index_phi(k));
% %          end         
% %     end
% %     luminosity_tot(1,t)=dot(flux1,factor);
% %     t_color_tot(1,t)= sum((flux1.*factor.*t_loc.*etha.^2)/luminosity_tot(1,t))
% %     luminosity_phi0(1,t)=dot(flux1,factor_phi0);
% %     t_color_phi0(1,t)= sum((flux1.*factor_phi0.*t_loc.*etha.^2)/luminosity_phi0(1,t));
% %     
% %     tau_eff=sqrt(3.*tau.*tauff);
% %     close;
% %     a=get(gcf,'Position');
% %     x0=15;
% %     y0=15;
% %     width=650;
% %     height=650;
% %     
% %     myFigure = figure('PaperPositionMode','auto','Position',a);
% %       
% %     set(myFigure,'units','points','position',[x0,y0,width,height])   
% %     p3=plot(xxx,y,'r','Linewidth',1); hold on
% %     h1=contour3(xx,yy,tau_eff,[2/3 2/3],'LineColor','b','Linewidth',1);hold on
% %     h1=contour3(xx,yy,tau,[1 1],'LineColor','m','Linewidth',1);hold on
% %  %   axis equal
% %     scatter(xdiff,ydiff,0.5,'g');hold on
% %     
% %     %h2=contour(xx,yy,d,[1 1],'LineColor','g');
% %     %h1=contour3(xx,yy,etha_all,[1 1],'LineColor','k','Linewidth',1);hold on 
% %     contour3(xx,yy,dsmooth,[1 1],'LineColor','k','Linewidth',1);  
% %     h=surf(xx,yy, log10(tau_eff));hold on
% %     grid off
% %     set(h,'LineStyle','none');
% %     colormap jet
% %     alpha(.2)
% %     colorbar
% %     caxis auto
% %     view(2);
% %     axis([0,2*rconv,0,2*rconv])
% %     title(['BSG model at t=' num2str(time.time1(t)*tconv,3) ' (sec)']);
% %     xlabel('cm'); 
% %     ylabel('cm');
% %     legend('R/R_*=1','Thermalization Depth','\tau =1 ','Diffusion Front','\eta =1','Log \tau_{eff}','Location','best');
% %     set(gca,'LineWidth',2,'FontSize',12);
% %     %set(gca,'XDir','reverse')
% %     name=['/home/nilou/Data/plot/BSGdiff/colortempContourWdiff_' int2str(t)];
% %     print('-dpng',name) 
% %     export_fig(name, '-png')
end
%save('/home/nilou/Data/processeddata/BSG/colortemp_tot_v2.mat','t_color_tot')
%save('/home/nilou/Data/processeddata/BSG/colortemp_phi0_v2.mat','t_color_phi0')