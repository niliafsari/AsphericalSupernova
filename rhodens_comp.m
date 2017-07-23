mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


pconv=etot/rtot^3;
rhoconv=mtot/rtot^3;
vconv=sqrt(etot/mtot);
data=load('/home/nilou/Data/lphi_data_v3.mat');
data=data.cursor_info;
dat=zeros(83,3);
dens_comp=zeros(83,1);
I1=zeros(83,1);
for i=1:length(data)
    dat(i,:)=data(i).Position;
end
phi=atan(dat(:,1)./dat(:,2));
[phi,I]=sort(phi);
dat=dat(I,:);
dens=load('/home/nilou/Data/rawdata/density/dens2048_0.csv');
for i=1:length(dat)
    i
    %phi=round(atan(dat(i,1)./dat(i,2))*2048/(pi/2));
    %r=round(sqrt(dat(i,1).^2+dat(i,2).^2)*2048/2);
    phi=atan(dat(i,1)./dat(i,2))
    r=round((0.5-0.5*f(phi))*2048/2)
    phi=round(atan(dat(i,1)./dat(i,2))*2048/(pi/2));
    dens_comp(i,1)=dens(r,phi)
end
close all
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
phi=atan(dat(:,1)./dat(:,2));
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');

set(myFigure,'units','points','position',[x0,y0,width,height])  
%set(myFigure,'units','points','position',[x0,y0,width,height])
%set(h,'LineStyle','none');
%caxis([-3 10])
%p2=plot(X,Y,'m','Linewidth',1), hold on,

%plot(phi,dens_comp/rhoconv,'r','Linewidth',1),hold on,

f1=fit(phi,dens_comp/rhoconv,'poly4')
plot(f1,phi,dens_comp/rhoconv)

rho=dens_comp/rhoconv
save('rhophi_v3.mat','rho');
csvwrite('rhophi_v3.csv',rho)
%plot(phi,3*I1*1e-4,'b','Linewidth',1);
%axis equal 
%axis([0,2,0,2])
xlabel('\phi (rad)'); 
ylabel('\rho_\phi / \rho_*');
%title('Time of Maximum Compression Rate');
legend('Location','northwest');
set(gca,'LineWidth',2,'FontSize',10);
%set(gca,'XDir','reverse')
name=['/home/nilou/Data/rhophi_comp_v3'];
print('-dpdf',name) 
export_fig(name, '-pdf')