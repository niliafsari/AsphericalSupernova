syms x
f=@(x) 1./(x.^(3/2).*(1-x.^-4))

close all
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=400;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])

y = logspace(0,4,1000);


% for i=1:length(y)
%     q(i) = integral(f,inf,y(i))
% end

h=loglog(y,q)
h.LineWidth=2
xlabel('x'); 
ylabel('G(x)-G(\infty)');
set(gca,'LineWidth',2,'FontSize',12);
name=['/home/nilou/Data/Gfun.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
