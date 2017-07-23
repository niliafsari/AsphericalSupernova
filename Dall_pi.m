
DRSG=load('DRSG_v3.mat','D');
DBSG=load('DBSG_v3.mat','D');
DIc=load('DIc_v3.mat','D');

DRSG_pi=load('DRSG_pi_v3.mat','D');
DBSG_pi=load('DBSG_pi_v3.mat','D');
DIc_pi=load('DIc_pi_v3.mat','D');


close all
figure
a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])
semilogy(phi,DRSG.D,'-.r','MarkerSize',3,'Linewidth',1),hold on,
semilogy(phi,DBSG.D,':','Linewidth',1.5),hold on,
semilogy(phi,DIc.D,'-<k','MarkerSize',3,'Linewidth',1),hold on,
semilogy(phi,DRSG_pi.D,'-<r','MarkerSize',3,'Linewidth',1),hold on,
semilogy(phi,DBSG_pi.D,'-b','MarkerSize',3,'Linewidth',1),hold on,
semilogy(phi,DIc_pi.D,'-xk','MarkerSize',3,'Linewidth',1),hold on,
semilogy(phi,ones(1,83),'--m','MarkerSize',3,'Linewidth',1.5),hold on,

xlabel('\theta [rad]'); 
ylabel('$\mathcal{D}$','Interpreter','latex');
axis([0.2,1.2,10^-7.5,10^7.5])
%title('Time of Maximum Compression Rate');
legend({'RSG $\alpha_f=\pi/2$','BSG  $\alpha_f=\pi/2$','Ic  $\alpha_f=\pi/2$','RSG $\alpha_f=\pi$','BSG  $\alpha_f=\pi$','Ic  $\alpha_f=\pi$', '$\mathcal{D}_{crit}=1$'},'Interpreter','Latex','Location','southeast');
legend('boxoff') 
set(gca,'LineWidth',1.5,'FontSize',10);
name=['/home/nilou/Data/Dall_v3'];
print('-dpdf',name) 
export_fig(name, '-pdf')








