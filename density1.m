
cell=2048;
densR1=zeros(1,76);
densR2=zeros(1,76);
densR3=zeros(1,76);
densR4=zeros(1,76);
densR1_b=zeros(1,76);
densR2_b=zeros(1,76);
densR3_b=zeros(1,76);
densR4_b=zeros(1,76);

for i=0:75
    name=['dens2048_' int2str(i) '.csv'] ;
    dens= csvread(name);
    densR1(i+1)=dens(cell/4,cell);
    densR2(i+1)=dens(2*cell/4,cell);
    densR3(i+1)=dens(3*cell/4,cell);
    densR4(i+1)=dens(4*cell/4,cell);
    densR1_b(i+1)=dens(cell/4,floor(cell/2));
    densR2_b(i+1)=dens(2*cell/4,floor(cell/2));
    densR3_b(i+1)=dens(3*cell/4,floor(cell/2));
    densR4_b(i+1)=dens(4*cell/4,floor(cell/2));
end
t=linspace(0,0.35,76);
figure;
subplot(2,1,1);
semilogy(t,densR1,'r','LineWidth',2);
hold on
semilogy(t,densR2,'b','LineWidth',2);
hold on
semilogy(t,densR3,'g','LineWidth',2);
hold on
semilogy(t,densR4,'c','LineWidth',2);
legend('R=0.5', 'R=1', 'R=1.5', 'R=2');
xlabel('time')
ylabel('density')
title('Density at \theta=90')
subplot(2,1,2);
plot(t,log10(densR1./densR1_b),'r','LineWidth',2);
hold on
plot(t,log10(densR2./densR2_b),'b','LineWidth',2);
hold on
plot(t,log10(densR3./densR3_b),'g','LineWidth',2);
hold on
plot(t,log10(densR4./densR4_b),'c','LineWidth',2);
legend('R=0.5', 'R=1', 'R=1.5', 'R=2')
xlabel('time')
ylabel('density')
title('$\rho_\theta$=90 / $\rho_\theta$=45')
