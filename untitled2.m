mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


pconv=etot/rtot^3;
rhoconv=mtot/rtot^3;
vconv=sqrt(etot/mtot);
data=load('/home/nilou/Data/lphi_data.mat');
data=data.cursor_info;
dat=zeros(89,3);
dens_comp=zeros(89,1);
I1=zeros(89,1);
for i=1:length(data)
    dat(i,:)=data(i).Position;
end
[~,I]=sort(dat(:,3));
dat=dat(I,:);
time=load('/home/nilou/Data/timesteps.mat');
for i=1:length(dat)
    i
    [~,I1(i,1)]=min(abs(time.time1-dat(i,3)));
    name=['/home/nilou/Data/dens2048_' int2str(I1(i,1)-1) '.csv'];
    dens=csvread(name);
    phi=round(atan(dat(i,1)./dat(i,2))*2048/(pi/2));
    r=round(sqrt(dat(i,1).^2+dat(i,2).^2)*2048/2);
    dens_comp(i,1)=dens(r,phi)
end
figure
%myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
%set(myFigure,'units','points','position',[x0,y0,width,height])
%set(h,'LineStyle','none');
%caxis([-3 10])
%p2=plot(X,Y,'m','Linewidth',1), hold on,
phi=atan(dat(:,1)./dat(:,2));
plot(phi,dens_comp/rhoconv,'r','Linewidth',1),hold on,

%plot(phi,3*I1*1e-4,'b','Linewidth',1);
%axis equal 
%axis([0,2,0,2])
xlabel('\phi (rad)'); 
ylabel('\rho_\phi / \rho_*');
%title('Time of Maximum Compression Rate');
%legend('Diffiusion Front','Progenitor','Location','southoutside','Orientation','horizontal');
set(gca,'LineWidth',2,'FontSize',12);
%set(gca,'XDir','reverse')
name=['/home/nilou/Data/rhophi_comp'];
print('-dpng',name) 
export_fig(name, '-png')