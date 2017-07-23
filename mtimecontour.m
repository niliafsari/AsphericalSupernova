theta=linspace(0,pi/2,2048);


radius=linspace(0,2,2048);
r=0.5;
%[TH,R] = meshgrid(theta,radius);
%[X,Y] = pol2cart(TH,R);
x=r*sin(theta);
y=r*cos(theta);

xx=zeros(2048,2048);
yy=zeros(2048,2048);
mtime_tot=zeros(2048,2048);
mtime_index=zeros(2048,2048);
% %theta=linspace(0,pi/2,2048);
 for i=1:2048
     for j=1:2048
        xx(i,j)=radius(i)*sin(theta(j));
        yy(i,j)=radius(i)*cos(theta(j)); 
     end
 end
for t=1:45
    t    
    radius=zeros(1,2048);
    name=['/home/nilou/Data/rawdata/mtime/mtime2048_' int2str(t) '.csv'] ;
    name1=['/home/nilou/Data/rawdata/dcodeunit/dcodeunit_' int2str(t) '.csv'] ;
    mtime= csvread(name);   
    d=csvread(name1);
    if (t==1)
        mtime_index=find(isnan(d) | isinf(d));  
        mtime_tot(mtime_index)=mtime(mtime_index);
    else
        mtime_tot(mtime_index)=mtime(mtime_index);
        if (mod(t,10)==0 ) 
            clear mtime_index
            mtime_index=find(isnan(d) | isinf(d)); 
        end
    end

end
csvwrite('/home/nilou/Data/rawdata/mtime/mtime_v2.csv',mtime_tot);
%mtime=csvread('/home/nilou/Data/rawdata/mtime/mtime_tot4.csv');
%close;
a=get(gcf,'Position');
x0=10;
y0=10;
width=550;
height=400;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])
v=linspace(0.08,0.14,100);
% mtime1=zeros(2048,2048);
% mtime1(isnan(d))=mtime(isnan(d));
h=contour(xx,yy, mtime_tot,v);hold on
%set(h,'LineStyle','none');
colormap jet
colorbar
%caxis([-3 10])
view(2);
%p2=plot(X,Y,'m','Linewidth',1), hold on,
p3=plot(x,y,'r','Linewidth',1);
axis equal 
axis([0,2,0,2])
xlabel('x'); 
ylabel('y');
title('Time of Maximum Compression Rate');
%legend('Diffiusion Front','Progenitor','Location','southoutside','Orientation','horizontal');
set(gca,'LineWidth',2,'FontSize',12);
%set(gca,'XDir','reverse')
%name=['/home/nilou/Data/mtimetot0'];
%print('-dpng',name) 
%export_fig(name, '-png')