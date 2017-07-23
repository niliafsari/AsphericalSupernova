mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


ptot=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048)/rtot;

xx=zeros(2048,2048);
yy=zeros(2048,2048);

thet=zeros(2048,2048);


for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
       thet(i,j)=theta(j);
       rad(i,j)=radius(i);
    end
end
time=load('/home/nilou/Data/timesteps.mat');
figure
% counter=0;
% for t=27:71
%     load(['/home/nilou/Data/processeddata/Collision/collision_comp' num2str(t) '.mat'],'e_znext','t_colnext','angleth_next','r_xnext'); 
%     t_accum(counter+1:length(t_colnext)+counter)=t_colnext;
%     r_accum(counter+1:length(t_colnext)+counter)=r_xnext;
%     e_accum(counter+1:length(t_colnext)+counter)=e_znext;
%     th_accum(counter+1:length(t_colnext)+counter)=angleth_next;
%     counter=length(t_colnext)+counter;
%     clear e_znext t_colnext angleth_next r_xnext
% end

load('/home/nilou/Data/processeddata/Collision/all_col_v2.mat','t_accum','r_accum','e_accum','th_accum'); 

anglethq=linspace(0,90,2048);
t_yq=linspace(-0.5,2,2048);

[X,Y] = meshgrid(anglethq,t_yq);
Vq = griddata(real(th_accum),real(log10(t_accum)),e_accum,X,Y, 'natural');

h=surf(X,Y,real(log10(Vq)));

set(h,'LineStyle','none');
axis([0 90 -0.5 2])
xlabel('\theta_{col}');
ylabel('Log t_{col}/t_*');
%caxis([log10(median(e_z)/((pi/4096)*(4/2048)))-3,log10(max(e_z)/((pi/4096)*(4/2048)))])
caxis([-10 -7])
colormap jet
c=colorbar('Location','eastoutside');
str = '$$ Log E_z/E_* $$';
ylabel(c,str,'Interpreter','latex','FontSize',12) 
name='/home/nilou/Data/plot/collision/all_surf.png';
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
print(gcf, '-dpng', '-r150', name)
export_fig(name, '-dpng', '-r150')
