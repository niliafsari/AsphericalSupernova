data=load('/home/nilou/Data/cursor_info1.mat');
data=data.cursor_info1;
dat=zeros(64,3);
for i=1:length(data)-1
    dat(i,:)=data(i).Position;
end
[~,I]=sort(dat(:,3));
dat=dat(I,:);
lphi=(0.5-(sqrt(dat(:,1).^2+dat(:,2).^2)))/0.5
save('lphi_v2.mat','lphi');
phi=atan(dat(:,1)./dat(:,2));
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
%p2=plot(X,Y,'m','Linewidth',31), hold on,
%scatter(phi,lphi,'+r');
f=fit(phi,lphi,'poly2')
plot(f,phi,lphi)
%axis equal 
%axis([0,2,0,2])
xlabel('\theta (rad)'); 
ylabel('l_\phi / R_*');
%title('Time of Maximum Compression Rate');
legend('l_\phi','fitted','Location','northwest');
set(gca,'LineWidth',1.5,'FontSize',10);
%set(gca,'XDir','reverse')
name=['/home/nilou/Data/lphi_comp_v2'];
print('-dpdf',name) 
export_fig(name, '-pdf')