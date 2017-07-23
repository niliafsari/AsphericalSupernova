
theta=linspace(0,pi/2,2048);


radius=linspace(0,2,2048);

%[TH,R] = meshgrid(theta,radius);
%[X,Y] = pol2cart(TH,R);

xx=zeros(2048,2048);
yy=zeros(2048,2048);
%theta=linspace(0,pi/2,2048);
for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end

for t=23:23
    radius=zeros(1,2048);
    name=['dens2048_' int2str(t) '.csv'] ;
    d= csvread(name);      
%     for i=1:2048
%         radius(i)=find(abs(d(:,i)-1)==min(abs(d(:,i)-1)),1)*(2/2048)*rconv;
%     end
%     d(isnan(d)| isinf(d))=1e-3;
%     X=radius.*sin(theta);
%     Y=radius.*cos(theta);
    %d(ii,jj)=1;
    close;
    a=get(gcf,'Position');
    x0=10;
    y0=10;
    width=550;
    height=400;
    
    myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
    %figure;
    set(myFigure,'units','points','position',[x0,y0,width,height])
%    h=surf(xx,yy, d);hold on
%    set(h,'LineStyle','none');
    h=surf(xx,yy, log10(d));hold on
    set(h,'LineStyle','none');
    colormap jet
    colorbar
    %caxis([-3 10])
    view(2);
    %p2=plot(X,Y,'m','Linewidth',1), hold on,
    %p3=plot(x,y,'r','Linewidth',1);
    axis equal 
    axis([0,2,0,2])
    xlabel('x'); 
    ylabel('y');
    %legend('Diffiusion Front','Progenitor','Location','southoutside','Orientation','horizontal');
    set(gca,'LineWidth',2,'FontSize',12);
    %set(gca,'XDir','reverse')
    name=['/home/nilou/Data/densmtime_' int2str(t)];
    print('-dpng',name, '-r1200') 
    export_fig(name, '-png','-r1200')
end