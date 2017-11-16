m=4.825996760367126832e-3*2;
e=9.199700591776255711e-4*2;
r=0.5;
rho_c=m/r^3;
v_c=r/sqrt(e/m);
t=1;
name=['/home/nilou/bgq/raw_data/dens16384_' int2str(t-1) '.csv'] ;
dens1= csvread(name)/rho_c;
name=['/home/nilou/bgq/raw_data/velr16384_' int2str(t-1) '.csv'] ;
vel1= csvread(name)/v_c;

t=97;
name=['/home/nilou/bgq/raw_data/dens16384_' int2str(t-1) '.csv'] ;
dens= csvread(name)/rho_c;

name=['/home/nilou/bgq/raw_data/velr16384_' int2str(t-1) '.csv'] ;
vel= csvread(name)/v_c;

name=['/home/nilou/bgq/raw_data/pres16384_' int2str(t-1) '.csv'] ;
pres= csvread(name);

theta=linspace(0,pi/2,2048);


radius=linspace(0,8,16384)/0.5;
rr=linspace(0,1,2048)';
dens_save=dens(1:2048);
vel_save=vel(1:2048);
pres_save=pres(1:2048);
size(rr)
size(dens_save)
data=[rr dens_save pres_save vel_save];
csvwrite('/home/nilou/poly_0.75.csv',data);

for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end

a=get(gcf,'Position');
x0=15;
y0=15;
width=350;
height=300;
close all
% 
myFigure = figure('PaperPositionMode','auto','Position',a,'Color','w');
set(myFigure,'units','points','position',[x0,y0,width,height])
y = logspace(-16,2,50);
rad=ones(1,50)*0.75;
plot(radius,log10((radius'.^3).*dens1),'r');hold on
plot(radius,log10((radius'.^3).*dens),'b');hold on
plot(radius,log10(vel1),'m');hold on
plot(radius,log10(vel),'g');hold on

xlabel('r/R_*')
ylabel('log r^3 \rho/\rho_* ')
legend('Original','Spherical at r/R_*=0.75','log v/v_* at t=0','log v/v_* at r=0.75' )
plot(rad,log10(y),'k');hold on
axis([0 2 -16 2])
