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
time=load('/home/nilou/Data/processeddata/timesteps.mat');
figure
counter=0;

% for t=33:71
%     t
%     load(['/home/nilou/Data/processeddata/Collision/collision_comp_v3_' num2str(t) '.mat'],'e_znext','t_colnext','angleth_next','r_xnext'); 
%     t_accum(counter+1:length(t_colnext)+counter)=t_colnext;
%     r_accum(counter+1:length(t_colnext)+counter)=r_xnext;
%     e_accum(counter+1:length(t_colnext)+counter)=e_znext;
%     th_accum(counter+1:length(t_colnext)+counter)=angleth_next;
%     counter=length(t_colnext)+counter;
%     clear e_znext t_colnext angleth_next r_xnext
% end
% 
% 
% 

load('/home/nilou/Data/processeddata/Collision/all_col_v3.mat','e_zq','t_accum','r_accum','e_accum','th_accum'); 

 t_yq=linspace(-0.5,2,2048);
 angle_thq=linspace(0,90,2048);

[angle_tang_m,t_snap_m] = meshgrid(angle_thq,t_yq);
load('/home/nilou/Data/processeddata/e_zqq.mat','e_zqq')
e_zqq=e_zqq';

[row,col,v]=find(e_zqq ~=0);

mu=[row(:), col(:)];

%z   = zeros(size(angle_tang_m));
load('/home/nilou/Data/processeddata/z.mat','z');
z= z* ((90/2048)*(2.5/2048));
z(z<10^-13)=10^-13;
% for i = 109668:length(mu)
%      i
%      a= e_zqq(mu(i,1),mu(i,2))/(2*pi*0.0439*0.0012);
%      z= z + a*exp(-(((angle_tang_m-angle_tang_m(mu(i,1),mu(i,2))).^2)/(2*0.0439^2)+((t_snap_m-t_snap_m(mu(i,1),mu(i,2))).^2)/( 2* 0.0012^2)));
% end

% for i=1693:2047
%     for j=1:2047
%         e_zqq(i,j)=sum(e_accum( t_yq(1,j)<log10(t_accum) & log10(t_accum)<t_yq(1,j+1) & angle_thq(1,i)<th_accum &  th_accum< angle_thq(1,i+1)));
%     end
% end
% 
load('/home/nilou/Data/processeddata/e_zqq.mat','e_zqq')
close all
angle_tang_m=90-angle_tang_m;
%plot(t_yq(2:2048),log10(cumsum(e_zq/(2.5/2048))));
%h=scatter(th_accum,log10(t_accum),2,log10(e_accum),'filled');
h=surf(angle_tang_m,t_snap_m,log10(z));
set(h,'LineStyle','none');
ylabel('Log t_{*,col}');
xlabel('\theta_{c}');
%ylabel('Log [(dE_z t_*)/ (E_* dlogt_{col}) ]');

set(gca,'LineWidth',2,'FontSize',10);

axis([0 90 -0.5 2])
% xlabel('Log t_{col}/t_*');
% ylabel('Log E_z/E_*');
%caxis([log10(median(e_z)/((pi/4096)*(4/2048)))-3,log10(max(e_z)/((pi/4096)*(4/2048)))])
%caxis([-9 -3])
colormap jet
c=colorbar('Location','eastoutside');
%str = '$$ Log E_z $$';
%ylabel(c, 'Log E_z/E_*');
ylabel(c,'Log [dE_{*,z}/  dlogt_{*,col} d\theta_{c} ]');
name='/home/nilou/Data/plot/collision/Col_E.png';
set(gca,'LineWidth',2,'FontSize',16);
set(gcf,'Color','w');
set(gcf, 'PaperPosition', [0 0 12 9]);
set(gcf, 'PaperSize', [14 10]);
%print('-dpdf',name, '-r300') 
% export_fig(name, '-pdf')

