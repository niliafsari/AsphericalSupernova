mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
ttot=rtot/sqrt(etot/mtot);
rhotot=mtot/rtot^3;
msun=1.989e33;
rsun=6.955e10;

poly_orig=csvread('poly.csv');
poly_changed=csvread('poly_gamma133.csv');
close all;
a=get(gcf,'Position');
x0=15;
y0=15;
width=380;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

plot(poly_changed(:,1)/rtot,poly_changed(:,2)/rhotot,'-b' ,'LineWidth',2); hold on
plot(poly_orig(:,1)/rtot,poly_orig(:,2)/rhotot,'-.r' ,'LineWidth',2)
axis([0 1 0 80])
set(gca,'LineWidth',2,'FontSize',13);

xlabel('r/R_*'); 
ylabel('\rho/\rho_*');
legend('Sedov at t/t_*=0.005','Sedov n=3','northeast');
name='/home/nilou/Data/poly.pdf';
print('-dpdf',name) 
export_fig(name, '-pdf')