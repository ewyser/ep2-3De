clear all,close all,clc,clf
set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;
%opengl hardware

nmp  = 64000
name = ['Exp1a_D_',num2str(nmp),'np.mat']
load(name)


%% DISPLAY
c1=meD.L(3)/2-meD.h(1)/2;
c2=meD.L(3)/2+meD.h(1)/2;
[np] = find(x(:,2)>0);

fig1=figure(1);
set(fig1,'Units','pixels','Position',[41 703 918 253]);
subplot(2,2,[1,3])
    scatter3(x(np,1)*1000,x(np,2)*1000,x(np,3)*1000,4.0,log10(1.0+epII(np)),'filled');
    
    axis equal;
    xlabel('$x$ [mm]');
    zlabel('$z$ [mm]');
    zticks([0 0.05 0.1]*1000);
    xlim([0 0.5]*1000);
    zlim([0 lz]*1000);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    cb=colormap(jet);
        cb=colorbar('FontSize',10,'TickLabelInterpreter','latex','Fontsize',10,'Location','south');
        cb.Position         =[0.13 0.7 0.1 0.075];
        cb.Label.String     ='$\log_{10}(1+\epsilon^p_{\mathrm{eqv}})$ [-]';
        cb.Label.Units      ='normalized';
        cb.Label.Position   =[1.7 0.0];
        cb.Label.FontSize   =10;
        cb.Label.Interpreter='Latex';
        cb.Ticks            = [0 0.5 1]
    caxis([0 1]);
    box on;
    view(0,0);
subplot(2,2,[2,4]);
    scatter3(x(np,1)*1000,x(np,2)*1000,x(np,3)*1000,4.0,log10(1.0+epII(np)),'filled');
    caxis([0 1]);
    axis equal;
    xlabel('$x$ [mm]');
    ylabel('$y$ [mm]');
    zlabel('$z$ [mm]');
    xlim([0 0.5]*1000);
    zlim([0 lz]*1000);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    box on;
        ax = gca;
        ax.BoxStyle = 'full';
    view(47,13);
name = ['fig3DElastoPlasticCollapseEquivStr' num2str(mpD.n)]
print(gcf,name,'-depsc');
print(gcf,name,'-dpng');


%% DISPLAY
c1=meD.L(3)/2-meD.h(1)/2;
c2=meD.L(3)/2+meD.h(1)/2;
[np] = find(x(:,2)>0);

fig2=figure(2);
set(fig2,'Units','pixels','Position',[400 461 509 253]);
    scatter3(x(np,1)*1000,x(np,2)*1000,x(np,3)*1000,4.0,log10(1.0+epII(np)),'filled');
    caxis([0 1]);
    axis equal;
    xlabel('$x$ [mm]');
    ylabel('$y$ [mm]');
    zlabel('$z$ [mm]');
    xlim([0 0.5]*1000);
    zlim([0 lz]*1000);
        cb=colormap(jet);
        cb=colorbar('FontSize',10,'TickLabelInterpreter','latex','Fontsize',10,'Location','south');
        cb.Position         =[0.6 0.8 0.15 0.075];
        cb.Label.String     ='$\log_{10}(1+\epsilon^p_{\mathrm{eqv}})$ [-]';
        cb.Label.Units      ='normalized';
        cb.Label.Position   =[1.8 0.0];
        cb.Label.FontSize   =10;
        cb.Label.Interpreter='Latex';
        cb.Ticks            = [0 0.5 1]
    caxis([0 1]);
    set(gca,'FontSize',12,'TickLabelInterpreter','latex');
    box on;
    grid on,grid minor
        ax = gca;
        ax.BoxStyle = 'full';
    view(47,13);
name = ['fig3DElastoPlasticCollapseEquivStr' num2str(mpD.n)]
print(gcf,name,'-depsc');
print(gcf,name,'-dpng');

%%
fs=load('Data_Buietal_2008_experimental_failure_surface.txt');
fs= [fs(:,1) zeros(size(fs,1),1) fs(:,2)].*1000;
S =load('Data_Buietal_2008_experimental_surface.txt');
S= [S(:,1) zeros(size(S,1),1) S(:,2)].*1000;


disp = du;
minu = 5e-4;
disp(disp>minu)=1;
disp(disp<minu)=0;
minu = minu/1e-3;
fig3=figure(3)
set(fig3,'Units','pixels','Position',[85 604 560 248]);
hold on
ax3=scatter3(x(np,1).*1000,x(np,2).*1000,x(np,3).*1000,4.0,(disp(np)),'filled');
map = [0 1 0;1 0 0]; colormap(map);
%pp.cbchar='Region';
pp.cbpos =[0.3200 0.7500 0.4009 0.0592];
pos            = pp.cbpos;
pp.cblpos=[pos(1) pos(2)+2];

cb=colorbar('FontSize',10,'TickLabelInterpreter','latex','Fontsize',10,'Location','SouthOutside');
cb.Position         =pp.cbpos;
%cb.Label.String     =pp.cbchar;
cb.Label.FontSize   =10;
cb.Label.Interpreter='Latex';
%cb.Label.Position   =pp.cblpos;
cb.Ticks      = [0.25 0.75];
cb.TickLabels = {['$u_p \leq',num2str(minu),'$ mm'],['$u_p>',num2str(minu),'$ mm']};
        cb.Title.String     =['ep3De v1.0, $n_{mp}=',num2str(mpD.n),'$']
        cb.Title.Interpreter='Latex';

view(0,0);
ax2=plot3(fs(:,1),fs(:,2),fs(:,3),'b:','LineWidth',2)
ax1=plot3(S(:,1),S(:,2),S(:,3),'b--','LineWidth',2)
box on;
tit = {'Experiment: final geometry','Experiment: failure surface'};
h1=legend([ax1 ax2 ax3],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.4407 0.5033 0.4514 0.1730],'NumColumns',1);
hold off;
axis equal;
xlabel('$x$ [mm]');
zlabel('$z$ [mm]');
xlim([0 0.5].*1000);
    zticks([0 0.05 0.1].*1000);
    zlim([0 lz].*1000);
grid on;
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
name = ['fig3DElastoPlasticCollapse_nmp_' num2str(mpD.n)];
print(gcf,name,'-depsc');
print(gcf,name,'-dpng');
