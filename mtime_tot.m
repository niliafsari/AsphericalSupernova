theta=linspace(0,pi/2,2048);


radius=linspace(0,2,2112);

%[TH,R] = meshgrid(theta,radius);
%[X,Y] = pol2cart(TH,R);

xx=zeros(2048,2048);
yy=zeros(2048,2048);
mtime_tot=zeros(2048,2048);
mtime_index=zeros(2048,2048);
%theta=linspace(0,pi/2,2048);
for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end

for t=0:50
    radius=zeros(1,2048);
    name=['mtime2048_' int2str(t) '.csv'] ;
    name1=['dcodeunit_' int2str(t) '.csv'] ;
    mtime= csvread(name);   
    d=csvread(name1);
    if (t==0)
        mtime_tot=mtime;
        mtime_index=find(isnan(d) | isinf(d));  
    else
        mtime_tot(mtime_index)=mtime(mtime_index);
        mtime_index=find(d(isnan(d) | isinf(d)));          
    end

end
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
    h=surf(xx,yy, mtime_tot);hold on
    set(h,'LineStyle','none');
    colormap jet
    %colorbar
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
    name=['/home/nilou/Data/mtimetot'];
    print('-dpng',name) 
    export_fig(name, '-png')