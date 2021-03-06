mtot= 4.8261095e-3*2;
etot=1.31610869e-2*2;
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
    

luminosity=zeros(1,501-201);
%load('luminosity.mat', 'luminosity');
time=load([path '/processeddata/timestepsni_1024.mat']);
clear diff_bsg
for t=289:303
    t
    name=[ path '/rawdata/dcodeunit/dcodeunit_ni1024_' int2str(t-1) '.mat'];
    load(name,'d');
    diff_frontBSG= d*(kappa/c)*rhoconv*vconv*rconv;

    if (t<302)
        ke=5;
    else
        ke=35;
    end
    temp=zeros(2048+2*ke,2048+ke);
    temp(1:2048,1:2048)=diff_frontBSG;
    kernel = 1*fspecial('disk', ke);
    temp1 = imfilter(log10(real(temp)),kernel,'same');
    dsmooth= temp1(1:2048,1:2048);

%    close all
%  
%     h=surf(xx,yy,log10(real(diff_frontBSG)));hold on
%     grid off
%     set(h,'LineStyle','none');
%     colormap jet
%     cl=colorbar;
%     %axis([0 4 0 4])
%     xlabel('x/R_*');
%     ylabel('y/R_*');
%     title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']); 
%     axis equal       
    
    [diff_bsg,~]=contour3(xx,yy,dsmooth,[0 0],'LineColor','r','Linewidth',1);
    set(gca,'LineWidth',2,'FontSize',12);
    caxis([-2 7]) 
    %ylabel(cl,'$ Log \frac{t_{diff}}{t_{dyn}} $','interpreter','latex','FontSize',24);    
    view(2); 
    name=[path '/plot/RSGdiff/DiffFrontRSGni_1024_' num2str(t-1) '.png'];
    print('-dpng', name)
    export_fig(name, '-dpng')
    save([path '/processeddata/RSG/diff_rsgni_1024_' num2str(t-1) '.mat'], 'diff_bsg') 
end