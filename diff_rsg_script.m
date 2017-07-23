clear all
mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
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
clear xdiff
clear ydiff
clear rdiff
clear phidiff
clear diff_rsg
diff_rsg=cell(50);
for t=1:50
    t
    name=['/home/nilou/Data/processeddata/RSG/dparamRSG_' int2str(t-1) '.csv'] ;
    diff_frontRSG= csvread(name);
    [C,~]=contour3(xx,yy,diff_frontRSG,[1 1],'LineColor','b','Linewidth',0.3);hold on     
    rdiff=sqrt(C(1,2:length(C)-1).^2+C(2,2:length(C)-1).^2);
    xdiff=C(1,2:length(C)-1);
    ydiff=C(2,2:length(C)-1); 
    phidiff=atan(xdiff./ydiff);
    [phidiff,I]=sort(phidiff);
    rdiff=rdiff(I);
    xdiff=xdiff(I);
    ydiff=ydiff(I);
    index_r=floor(rdiff/((2/2048)*rconv));
    index_phi=floor(atan(xdiff./ydiff)/(0.5*pi/2048));
    index_phi(index_phi==0)=1;
    index_r(index_r==0)=1;
    close all
    
    y0=0.5*rconv;
    x0=0;
    size(xdiff)
     rdiff_new=sqrt((xdiff-x0).^2+(ydiff-y0).^2);
     phidiff_new=atan((xdiff-x0)./(ydiff-y0));
     oi=0;
     num=1000;
     rem=num;
    for i=1:ceil(length(rdiff)/num) 
       if (i==ceil(length(rdiff)/num))
           rem=mod(length(rdiff),num)
       end
       X=phidiff_new((i-1)*num+1:(i-1)*num+rem);
       Y=rdiff_new((i-1)*num+1:(i-1)*num+rem);     
        X=X(Y<=prctile(Y,60));
        Y=Y(Y<=prctile(Y,60));
         X=X(Y>= prctile(Y,15));
         Y=Y(Y>= prctile(Y,15));
%         f=fit(X',Y','fourier1');
%         Y=f(X);
%        Y(1:length(Y))=mean(Y);
        X_new(oi+1:oi+length(X))=X;
        Y_new(oi+1:oi+length(Y))=Y;
        oi=length(X_new);
    end
    X_new(X_new<0)=-((X_new(X_new<0)+(pi/2)));
    xdiff=Y_new.*sin(X_new);
    ydiff=Y_new.*cos(X_new);     
    xdiff(X_new<0)=Y_new(X_new<0).*cos(X_new(X_new<0));
    ydiff(X_new<0)=Y_new(X_new<0).*sin(X_new(X_new<0));
    

    size(xdiff)
     ydiff=ydiff+y0;
    ydiff=ydiff(ydiff>=0);
    xdiff=xdiff(ydiff>=0);
    ydiff=ydiff(xdiff>=0);
    xdiff=xdiff(xdiff>=0);
    ydiff=ydiff(ydiff<=2*rconv);
    xdiff=xdiff(ydiff<=2*rconv);
    ydiff=ydiff(xdiff<=2*rconv);
    xdiff=xdiff(xdiff<=2*rconv);
    
    size(xdiff)
%     a=get(gcf,'Position');
%     x0=15;
%     y0=15;
%     width=450;
%     height=450;
%     myFigure = figure('PaperPositionMode','auto','Position',a);    
%     set(myFigure,'units','points','position',[x0,y0,width,height]) ;
%     scatter(xdiff,ydiff,1,'.');
%     axis([0,2*rconv,0,2*rconv])
%     grid off  
     
    diff_rsg{t}(1,:)=xdiff;
    diff_rsg{t}(2,:)=ydiff;
     
     clear xdiff
     clear ydiff
     clear Y_new
     clear X_new

% %     %colormap jet
% %     %c=colorbar('Location','eastoutside');
     
%     %ylabel(c,'log (g_\epsilon)','FontSize',12) 
%     %caxis auto
%     view(2);
%     legend_entries{1}='RSG';
%     legend_entries{2}='BSG';
%     legend_entries{3}='Ic';
%     
%     %title(['Ic Model t=' int2str(t)]);
     xlabel('x'); 
     ylabel('y');
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
save('diff_rsg.mat','diff_rsg');