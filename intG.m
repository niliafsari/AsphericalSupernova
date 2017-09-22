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

y = logspace(-5,4,50);

q=zeros(1,50);
for i=1:length(y)
    q(i) = integral(f,1+y(i),1.02)
end

%h=semilogx(y,q)

% 
f= fit(q',log10(y)','poly7');
yy=f(-100:1000);

plot(f,q',log10(y)')


%axis([10^-5 10^4 -5 5])

ylabel('x-1'); 
xlabel('Log [G(10^4)-G(x)]');
set(gca,'LineWidth',2,'FontSize',12);
name=['/home/nilou/Data/Gfun_v1.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')
