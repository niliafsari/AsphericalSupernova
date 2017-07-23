mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
ttot=rtot/sqrt(etot/mtot);
gridsize=2048;
vtot=sqrt(etot/mtot);
pattern = csvread('/home/nilou/Data/patternspeed_raw.csv');
shock = csvread('/home/nilou/Data/shockspeed_raw.csv');
phi_array = csvread('/home/nilou/Data/phi_array.csv');

plot(pattern(gridsize/4-10,1:3*gridsize/4),'b'), hold on,
plot(shock(gridsize/4-10,1:3*gridsize/4),'r');


set(gca,'LineWidth',2,'FontSize',12);

xlabel('\phi'); 
ylabel('v/v_*');
%title('Time of Maximum Compression Rate');
legend('v_\phi','v_{sh}','Location','northeast');
set(gca,'LineWidth',1.5,'FontSize',10);
name=['/home/nilou/Data/pattern.pdf'];
print('-dpdf',name) 
export_fig(name, '-pdf')