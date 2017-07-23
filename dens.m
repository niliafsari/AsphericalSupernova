cell=2048;
densR1=zeros(75);
densR2=zeros(75);
densR3=zeros(75);
densR4=zeros(75);
densR1_b=zeros(75);
densR2_b=zeros(75);
densR3_b=zeros(75);
densR4_b=zeros(75);

for i=0:75
    name=['dens2048_' int2str(i) '.csv'] ;
    dens= csvread(name);
    densR1(i+1)=dens(cell/4,cell);
    densR2(i+1)=dens(2*cell/4,cell);
    densR3(i+1)=dens(3*cell/4,cell);
    densR4(i+1)=dens(4*cell/4,cell);
    densR1_b(i+1)=dens(cell/4,cell/3);
    densR2_b(i+1)=dens(2*cell/4,2*cell/3);
    densR3_b(i+1)=dens(3*cell/4,2*cell/3);
    densR4_b(i+1)=dens(4*cell/4,2*cell/3);
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
subplot(2,2,1);
plot(t,densR1./densR1_b,'r','LineWidth',2);
hold on
plot(t,densR2./densR2_b,'b','LineWidth',2);
hold on
plot(t,densR3./densR3_b,'g','LineWidth',2);
hold on
plot(t,densR4./densR4_b,'c','LineWidth',2);
legend('R=0.5', 'R=1', 'R=1.5', 'R=2')
xlabel('time')
ylabel('density')
title('$\rho_\theta$=90 / $\rho_\theta$=60')
