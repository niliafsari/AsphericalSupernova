mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
ttot=rtot/sqrt(etot/mtot);
gridsize=2048;
vtot=sqrt(etot/mtot);
 pattern = csvread('/home/nilou/Data/patternspeed_raw.csv');
 shock = csvread('/home/nilou/Data/shockspeed_raw.csv');
 phi_array = csvread('/home/nilou/Data/phi_array.csv');
close all

figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

x=phi_array(0.17*gridsize/4:1*gridsize/4);
y=log10(pattern(gridsize/4-6,0.17*gridsize/4:1*gridsize/4));
y(isinf(y))=1;
 f = fit(phi_array(0.17*gridsize/4:1*gridsize/4),y','poly8')


fdata = feval(f,x);
I = abs(fdata - y') > 0.3*std(y');
outliers = excludedata(x,y','indices',I);
sum(outliers)
y(outliers)=[];
x(outliers)=[];

h1=plot(f,x',y','.b'),hold on,
set(h1,'MarkerSize',6,'color','b','LineWidth',1.5);
 
%plot(phi_array(0.17*gridsize/4:1*gridsize/4),log10(pattern(gridsize/4-10,0.17*gridsize/4:1*gridsize/4)),'b'), hold on,
plot(phi_array(0.17*gridsize/4:1*gridsize/4),log10(shock(gridsize/4-10,0.17*gridsize/4:1*gridsize/4)),'--r','LineWidth',1.5);


set(gca,'LineWidth',2,'FontSize',12);

xlabel('\theta [rad]'); 
ylabel('log(v/v_*)');
%title('Time of Maximum Compression Rate');
legend('v_\phi', 'fitted','v_s','Location','northeast');
set(gca,'LineWidth',1.5,'FontSize',10);
name=['/home/nilou/Data/pattern.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')