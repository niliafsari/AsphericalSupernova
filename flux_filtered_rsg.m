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
kappa=0.2;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);


r=0.5*rconv;
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
time=load('/home/nilou/Data/timesteps.mat');
load('diff_rsg.mat','diff_rsg');
for t=1:50
    t
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
    density= csvread(name)*rhoconv;
    name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
    pres= csvread(name)*pconv;
%     for i=1:2048
%         for j=1:2048
%             if (i==1)
%                 pr=r(i); 
%                 pindex_r=i;
%             else
%                 pr=r(i-1);
%                 pindex_r=i-1;
%             end
%             if (j==1)
%                 pt=theta(j);
%                 pindex_t=j;
%             else
%                 pt=theta(j-1);
%                 pindex_t=j-1;
%             end
%             if (i==2048)
%                 nr=r(i);
%                 nindex_r=i;
%             else
%                 nr=r(i+1);
%                 nindex_r=i+1;
%             end
%             if (j==2048)
%                 nt=theta(j);
%                 nindex_t=j;
%             else
%                 nt=theta(j+1);
%                 nindex_t=j+1;
%             end
%             dr=nr-pr;
%             cr=pr+ dr/2;
%             dt=nt-pt;
%             ct=pt+ dt/2;
%             gradp(i,j)=sqrt(((pres(nindex_r,j)-pres(pindex_r,j))/dr)^2+ ((pres(i,nindex_t)-pres(i,pindex_t))/(cr*dt))^2);
%         end
%     end
    name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-1) '.csv'];
    gradp=csvread(name)*(rconv_bsg/pconv_bsg)*(pconv/rconv);
    rdiff=sqrt(diff_rsg{t}(1,2:length(diff_rsg{t})).^2+diff_rsg{t}(2,2:length(diff_rsg{t})).^2);
    xdiff=diff_rsg{t}(1,2:length(diff_rsg{t}));
    ydiff=diff_rsg{t}(2,2:length(diff_rsg{t})); 
    phidiff=atan(xdiff./ydiff);
    [phidiff,I]=sort(phidiff);
    rdiff=rdiff(I);
    xdiff=xdiff(I);
    ydiff=ydiff(I);
    index_r=floor(rdiff/((2/2048)*rconv));
    index_phi=floor(atan(xdiff./ydiff)/(0.5*pi/2048));
    index_phi(isnan(index_phi))=2048;
    index_phi(index_phi==0)=1;
    index_r(index_r==0)=1;
    if  sum(index_r>2048)
        sum(index_r>2048)
    end
    index_r(index_r>2048)=2048;
    flux1=zeros(1,length(diff_rsg{t})-1);
    for k=1:length(index_r)
        flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
        if (k==1)
            luminosity(t)= luminosity(t)+(2*pi* rdiff(k)^2 * flux1(k) *(cos(phidiff(k))-cos(phidiff(k+1))));
        else
            luminosity(t)= luminosity(t)+(2*pi* rdiff(k)^2 * flux1(k) *(cos(phidiff(k-1))-cos(phidiff(k))));
        end
    end
     close all
%     a=get(gcf,'Position');
%     x0=15;
%     y0=15;
%     width=450;
%     height=450;
%     myFigure = figure('PaperPositionMode','auto','Position',a);    
%     set(myFigure,'units','points','position',[x0,y0,width,height])   
%      
% 
%     grid off  
%     %colormap jet
%     %c=colorbar('Location','eastoutside');
%     axis([0,4,0,4])
%     %ylabel(c,'log (g_\epsilon)','FontSize',12) 
%     %caxis auto
%     view(2);
%     legend_entries{1}='RSG';
%     legend_entries{2}='BSG';
%     legend_entries{3}='Ic';
%     
%     %title(['Ic Model t=' int2str(t)]);
%     xlabel('x/R_*'); 
%     ylabel('y/R_*');
%     %legend([hcont3 hcont2 hcont1],legend_entries,'Orientation','horizontal','Location','northoutside');
%     set(gca,'LineWidth',1,'FontSize',12);
%     (time.time1(t+1)/ttot)
%     name=['/home/nilou/Data/diff_all_' num2str(time.time1(t+1)/ttot) '.pdf'];
%     set(gcf, 'PaperPosition', [0 0 5 4.5]); %Position plot at left hand corner with width 5 and height 5.
%     set(gcf, 'PaperSize', [5 4.5]); %Set the paper to have width 5 and height 5.
%     print(gcf, '-dpdf', '-r900', name)
%     %print('-dpdf',name) 
%     %export_fig(name, '-pdf','r10')
end
save('/home/nilou/Data/luminosity_rsg_final.mat', 'luminosity')
