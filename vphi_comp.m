speed=load('/home/nilou/trapezium/patternspeed10_raw.csv');

%speed=log10(speed);
mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


pconv=etot/rtot^3;
rhoconv=mtot/rtot^3;
vconv=sqrt(etot/mtot);


speed=speed/vconv;
theta=linspace(0,3*pi/8,2048*3/4)';
speed=speed(theta>0.158);
theta=theta(theta>0.158);
speed=speed(theta<0.9);
theta=theta(theta<0.9);
speed=log10(speed)
f=fit(theta,speed,'gauss7')



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
lphi=lphi(I);
vphi=f(phi);
save('vphi_v3.mat','vphi');
close all
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

plot(f,theta,speed)
%axis equal 
%axis([0,2,0,2])
xlabel('\phi (rad)'); 
ylabel('log (v_\phi / v_*)');
%title('Time of Maximum Compression Rate');
legend('v_\phi','fitted','Location','northeast');
set(gca,'LineWidth',1.5,'FontSize',10);
%set(gca,'XDir','reverse')
name=['/home/nilou/Data/vphi_comp_v3'];
print('-dpdf',name) 
export_fig(name, '-pdf')