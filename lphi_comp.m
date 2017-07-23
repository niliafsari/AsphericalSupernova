mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


pconv=etot/rtot^3;
rhoconv=mtot/rtot^3;
vconv=sqrt(etot/mtot);
data=load('/home/nilou/Data/lphi_data_v3.mat');
data=data.cursor_info;
dat=zeros(83,3);
for i=1:length(data)
    dat(i,:)=data(i).Position;
end

lphi=(0.5-(sqrt(dat(:,1).^2+dat(:,2).^2)))/0.5;

phi=atan(dat(:,1)./dat(:,2));
[phi,I]=sort(phi);
dat=dat(I,:);
lphi=lphi(I,:);
save('lphi_v3.mat','lphi');

dens=load('/home/nilou/Data/rawdata/density/dens2048_0.csv');
for i=1:length(dat)
    i
    %phi=round(atan(dat(i,1)./dat(i,2))*2048/(pi/2));
    %r=round(sqrt(dat(i,1).^2+dat(i,2).^2)*2048/2);
    phi=atan(dat(i,1)./dat(i,2));
    r=round((0.5-0.5*f(phi))*2048/2);
    phi=round(atan(dat(i,1)./dat(i,2))*2048/(pi/2));
    dens_comp(i,1)=dens(r,phi);
end

close
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])
%set(h,'LineStyle','none');
%caxis([-3 10])
%p2=plot(X,Y,'m','Linewidth',1), hold on,
%scatter(phi,lphi,'+r');
close all
phi=atan(dat(:,1)./dat(:,2));
f=fit(phi,lphi,'poly3')
f1=fit(phi,dens_comp/rhoconv,'poly7')

[hAx,hLine1,hLine2]=plotyy([phi,phi],[lphi,f(phi)],[phi,phi],[log10(dens_comp/rhoconv),log10(f1(phi))])
hLine2(1).LineStyle='none';
hLine1(2).LineStyle='--';
hLine1(1).LineStyle='none';
hLine1(1).Marker='.';
hLine2(1).Marker='*';
hLine2(1).MarkerSize=4;
hLine1(1).MarkerEdgeColor='k';
hLine2(1).MarkerEdgeColor='k';
hLine1(2).Color='r';
hLine2(2).Color='b';
hLine1(2).LineWidth=1.5;
hLine2(2).LineWidth=1.5;

hAx(2).YColor=[0 0 1];
ylabel(hAx(1),'l_\phi / R_*') % left y-axis
ylabel(hAx(2),'Log (\rho_\phi/\rho_*)') % right y-axis

set(gca,'LineWidth',1.5,'FontSize',10);

%plot(f1,phi,dens_comp/rhoconv)
%axis equal 
%axis([0,2,0,2])
xlabel('\theta [rad]'); 
%ylabel('l_\phi / R_*');
%title('Time of Maximum Compression Rate');
legend('l_\phi','fitted l_\phi','\rho_\phi','fitted \rho_\phi','Location','northwest');
set(gca,'LineWidth',1.5,'FontSize',10);
%set(gca,'XDir','reverse')
name=['/home/nilou/Data/lphi_rho_comp'];
print('-dpdf',name) 
export_fig(name, '-pdf')