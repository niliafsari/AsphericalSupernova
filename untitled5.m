%% Create two axes
s1=subplot(2,1,1);
ax1 = s1;
[x,y,z] = peaks;
surf(ax1,x,y,z), hold on;
view(2)
ax2 = s1;
scatter(ax2,randn(1,120),randn(1,120),50,randn(1,120),'filled')
%% Link them together
linkaxes([ax1,ax2])
%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax1,'hot')
colormap(ax2,'cool')
%% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);