mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=5*msun;
e=1e51;
r=0.2*rsun;
kappa=0.2;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);

%r=linspace(0,2,2048)*rconv;
r=0.5*rconv;
%theta=linspace(0,pi/2,2048);

x=r*sin(theta);
y=r*cos(theta);


%[TH,R] = meshgrid(theta,r);
%[X,Y] = pol2cart(TH,R);

theta=linspace(0,pi/2,2048);
for t=0:11
    radius=zeros(1,2048);
    name=['dparam_' int2str(t) '.csv'] ;
    d= csvread(name);
    for i=1:2048
        
        radius(i)=find(abs(d(:,i)-1)==min(abs(d(:,i)-1)),1)*(2/2048)*rconv;
    end
    mean(min(abs(d-1)))
    xx=radius.*sin(theta);
    yy=radius.*cos(theta);
%     a=get(gcf,'Position');
%     close;
%     myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
%     plot(x,y,'b','Linewidth',1.5),hold on,
%     plot(xx,yy,'r','Linewidth',1.5)
%     axis([0,2*rconv,0,2*rconv])
%     title(['Ic Model t=' int2str(t)]);
%     xlabel('cm'); 
%     ylabel('cm');
%     legend('Progenitor','Diffiusion Front',1);
%     set(gca,'LineWidth',2,'FontSize',17);
%     %set(gca,'XDir','reverse')
%     name=['diffIc_' int2str(t)];
%     print('-dpng',name) 
%     export_fig(name, '-png')
end



%min(d)
%d(d>1e10) = 1e10;
%figure
%surf(X,Y,log10(d))