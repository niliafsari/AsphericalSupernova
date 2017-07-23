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
time=load('/home/nilou/Data/timesteps.mat');
for t=35:41
    t
    diff_bsg=load(['/home/nilou/Data/processeddata/BSG/diff_bsg_smoothed_' num2str(t) '.mat']);
     name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
     density= csvread(name)*rhoconv;
     name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
     pres= csvread(name)*pconv;
     name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-1) '.csv'];
     gradp=csvread(name);

    rdiff=sqrt(diff_bsg(1,2:length(diff_bsg)).^2+diff_bsg(2,2:length(diff_bsg)).^2);
    xdiff=diff_bsg(1,2:length(diff_bsg));
    ydiff=diff_bsg(2,2:length(diff_bsg)); 
    phidiff=atan(xdiff./ydiff);
    [phidiff,I]=sort(phidiff);
    rdiff=rdiff(I);
    xdiff=xdiff(I);
    ydiff=ydiff(I);
    scatter(xdiff,ydiff);
    index_r=floor(rdiff/((2/2048)*rconv));
    index_phi=floor(atan(xdiff./ydiff)/(0.5*pi/2048));
    index_phi(index_phi==0)=1;
    index_r(index_r==0)=1;
    if  sum(index_r>2048)
        sum(index_r>2048)
    end
    index_r(index_r>2048)=2048;
    flux1=zeros(1,length(diff_bsg)-1);
    for k=2:length(index_r)
         flux1(k)= (c/(kappa* density(index_r(k),index_phi(k))))*gradp(index_r(k),index_phi(k));
         luminosity(t)= luminosity(t)+(2*pi* rdiff(k)^2 * flux1(k) *(cos(phidiff(k-1))-cos(phidiff(k))));
    end
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=650;
    height=650;
       
    myFigure = figure('PaperPositionMode','auto','Position',a);
      
    set(myFigure,'units','points','position',[x0,y0,width,height])   
    scatter(xdiff,ydiff);
    xlabel('x/R_*');
    ylabel('y/R_*');
    title(['t=' num2str(time.time1(t)*tconv,3) ' (sec)']); 
    axis equal       
    set(gca,'LineWidth',2,'FontSize',12);

end
save('/home/nilou/Data/processeddata/BSG/luminosity_smoothed.mat','luminosity')
