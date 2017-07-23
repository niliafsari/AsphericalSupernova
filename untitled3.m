    close all
    a=get(gcf,'Position');
    x0=15;
    y0=15;
    width=580;
    height=500;
    myFigure = figure('PaperPositionMode','auto','Position',a);    
    set(myFigure,'units','points','position',[x0,y0,width,height]) 

    ax1=axes;
    h3=scatter(ax1,xx(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3)),yy(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3)),3,...
        log10(e_z),'filled');
    colormap(ax1,'gray')
    axis(ax1,[0 4 0 4])
    ax2=axes;
    h4=scatter(ax2,xx_next,yy_next,3,...
    log10(e_znext),'filled');
    linkaxes([ax1,ax2])
    axis(ax2,[0 4 0 4])
    colormap(ax1,'gray')
    colormap(ax2,'jet')
    caxis(ax2,[-10 -7])
    caxis(ax1,[-10 -7])
    str = 'log(E_z/E_*)';
    set(ax1,'LineWidth',2,'FontSize',12);
    set(ax2,'LineWidth',2,'FontSize',12);
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    ax1.XTick = [0 1 2 3 4];
    ax1.YTick = [0 1 2 3 4];
    title(ax1,['E_{cone}/E_*= ' num2str(sum(e_znext),3) ', E_{tot}/E_*=' num2str(totalenergy,3) ]);
    xlabel(ax1,'x/R_*'); 
    ylabel(ax1,'y/y_*');
    set([ax1,ax2],'Position',[.19 .11 .685 .815]);
    cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
    cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);