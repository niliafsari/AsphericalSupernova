

t=linspace(0,0.35,76);
figure
subplot(2,1,1);
semilogy(t,densR1,'r','LineWidth',2)
hold on
semilogy(t,densR2,'b','LineWidth',2)
hold on
semilogy(t,densR3,'g','LineWidth',2)
hold on
semilogy(t,densR4,'c','LineWidth',2)
legend('R=0.5', 'R=1', 'R=1.5', 'R=2')
xlabel('time')
ylabel('Entropy')
title('Entroy at \theta=90')
ax = gca; % handle of current axes
ax.FontSize = 18;
subplot(2,1,2);
semilogy(t,densR1./densR1_b,'r','LineWidth',2)
hold on
semilogy(t,densR2./densR2_b,'b','LineWidth',2)
hold on
semilogy(t,densR3./densR3_b,'g','LineWidth',2)
hold on
semilogy(t,densR4./densR4_b,'c','LineWidth',2)
legend('R=0.5', 'R=1', 'R=1.5', 'R=2')
xlabel('time')
ylabel('Entropy')
title('S_{\theta=90} / S_{\theta=45}')
ax = gca; % handle of current axes
ax.FontSize = 18;