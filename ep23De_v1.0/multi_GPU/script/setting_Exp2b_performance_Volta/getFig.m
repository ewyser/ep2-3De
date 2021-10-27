function [] = getFig(x,data,char,clim,np,lz,t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fig=figure(1);
set(fig,'Units','pixels','Position',[687 494 2*459 327]);
clf
subplot(2,2,1)
    scatter3(x(np,1),x(np,2),x(np,3),4.0,data,'filled');
    caxis([0 clim])
    axis equal;
    zlabel('$z$ [m]');
    zticks([0 5 10]);
    zlim([0 lz]);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    cb=colormap(jet);
        cb=colorbar('FontSize',10,'TickLabelInterpreter','latex','Fontsize',10,'Location','south');
        cb.Position         =[0.5 0.2 0.1 0.05];
        cb.Label.String     =char;
        cb.Label.Units      ='normalized';
        cb.Label.Position   =[0.5 2.2];
        cb.Label.FontSize   =10;
        cb.Label.Interpreter='Latex';
    box on;
    view(0,0);
    title('a) vertical profile','FontSize',12);
subplot(2,2,3);
    scatter3(x(np,1),x(np,2),x(np,3),4.0,data,'filled');
    caxis([0 clim])
    axis equal;
    xlabel('$x$ [m]');
    ylabel('$y$ [m]');
    yticks([0 5 10]);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    box on;
    view(0,90);
    title('b) aerial view','FontSize',12);
subplot(2,2,[2,4]);
    a = ones(size(data)); a(data<1) = 1;
    s = scatter3(x(np,1),x(np,2),x(np,3),4.0,data,'filled');
        s.MarkerFaceAlpha = 'flat';
        s.AlphaData       = a;
    axis equal;
    xlabel('$x$ [m]');
    ylabel('$y$ [m]');
    zlabel('$z$ [m]');
    yticks([0 5 10]);
    zticks([0 5 10]);
    zlim([0 lz]);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    box on;
        ax = gca;
        ax.BoxStyle = 'full';
    view(47,13);
    caxis([0 clim])
    title(['c) perspective at $t=',num2str(t,'%.2f'),'$ [s]'],'FontSize',12);
set(gcf,'color','white');
end

