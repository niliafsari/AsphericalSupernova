mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
ttot=rtot/sqrt(etot/mtot);
msun=1.989e33;
rsun=6.955e10;
m=14*msun;
e=1e51;
r=400*rsun;
c=3e10;
rconv_rsg=r/rtot;
econv=e/etot;
mconv=m/mtot;
v=sqrt(e/m);
rho=m/(r^3);
Z=0.005;
a=7.566e-15;
X=0.7;
kappa_T=0.2*(1+X);

m=15*msun;
e=1e51;
r=49*rsun;
c=3e10;
rconv_bsg=r/rtot;

m=5*msun;
e=1e51;
r=0.2*rsun;
kappa=0.2;
rconv_ic=r/rtot;


theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048)/0.5;

xx=zeros(2048,2048);
yy=zeros(2048,2048);

for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
time=load('/home/nilou/Data/timesteps.mat');
for t=23:50
    temp=load(['/home/nilou/Data/processeddata/ic/diff_ic_smoothed_' num2str(t) '.mat'], 'diff_bsg'); 
    diff_ic=temp.diff_bsg;
    temp=load(['/home/nilou/Data/processeddata/BSG/diff_bsg_smoothed_k55' num2str(t) '.mat'], 'diff_bsg');
    diff_bsg=temp.diff_bsg;
    temp=load(['/home/nilou/Data/processeddata/RSG/diff_rsg_smoothed55_' num2str(t) '.mat'], 'diff_bsg');
    diff_rsg=temp.diff_bsg;
    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=450;
    height=450;
    myFigure = figure('PaperPositionMode','auto','Color','w','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height])   
     
    scatter(diff_ic(1,:)/(rconv_ic*0.5),diff_ic(2,:)/(rconv_ic*0.5),1,'k');hold on 
    scatter(diff_bsg(1,:)/(rconv_bsg*0.5),diff_bsg(2,:)/(rconv_bsg*0.5),1,'b');hold on 
    scatter(diff_rsg(1,:)/(rconv_rsg*0.5),diff_rsg(2,:)/(rconv_rsg*0.5),1,'r');hold on 

    axis equal
    axis([0,4,0,4])
    view(2);
    legend_entries{1}='Ic';
    legend_entries{2}='BSG';
    legend_entries{3}='RSG';
    set(gca,'LineWidth',2,'FontSize',12);
    xlabel('x/R_*'); 
    ylabel('y/R_*');
    title(['t/t_*=' num2str(time.time1(t)/ttot,3) ' ']); 
    
    name=['/home/nilou/Data/plot/allquantity/diff_all_v6_' num2str(t) '.png'];
    set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [10 10]); %Set the paper to have width 5 and height 5.
    print(gcf, '-dpng', '-r300', name)
    %print('-dpdf',name) 
    export_fig(name, '-dpng','-r300')
end

