clear all,close all,clc,clf
set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;
%opengl hardware
nmp  = 12800;
name = ['Exp1a_D_',num2str(nmp),'np.mat'];
load(name);

%%
fig1=figure(1);
set(fig1,'Units','pixels','Position',[681 760 560 185]);
    d = sortrows([x(:,1) x(:,2) log10(1+epII)],3,'descend');
    d = ([x(:,1)*1000 x(:,2)*1000 log10(1+epII)]);
    hold on;
    scatter(d(:,1),d(:,2),5,d(:,3),'filled');
    X = [50 160 160 50 50];
    Y = [30 30 80 80 30];
    ax1=plot(X,Y,'r--','LineWidth',2)
    colormap('jet');
    hold off;
    axis equal;
	axis tight;
    xlabel('$x$ [mm]');
    ylabel('$z$ [mm]');
    xlim([0 0.5*1000]);
    box on;
    grid on;
    cb=colormap(jet);
        cb=colorbar('FontSize',12,'TickLabelInterpreter','latex','Fontsize',12,'Location','south');
        cb.Position         =[0.5 0.5 0.15 0.1];
        cb.Label.String     ='$\log_{10}(1+\epsilon^p_{\mathrm{eqv}})$ [-]';
        cb.Label.Units      ='normalized';
        cb.Label.Position   =[1.8 0.0];
        cb.Label.FontSize   =12;
        cb.Label.Interpreter='Latex';
        caxis([0 1]);
    yticks([0 0.05*1000 0.1*1000]);
    ylim([0 lz*1000]);
    drawnow
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
name = ['figElastoPlasticCollapseEquivStr2D' num2str(mpD.n)];
print(gcf,name,'-depsc');
print(gcf,name,'-dpng');


























fig2=figure(2);
np = find(x(:,1)>min(X/1000)&x(:,1)<max(X/1000)&x(:,2)>min(Y/1000)&x(:,2)<max(Y/1000))
set(fig2,'Units','pixels','Position',[681 760 560 219]);
    d = sortrows([x(:,1) x(:,2) log10(1+epII)],3,'descend');
    d = ([x(np,1)*1000 x(np,2)*1000 log10(1+epII(np))]);
    hold on;
    scatter(d(:,1),d(:,2),5,d(:,3),'filled');
    ax1=plot(X,Y,'r--','LineWidth',2)
    colormap('jet');
    hold off;
    axis equal;
	axis tight;
    xlabel('$x$ [mm]');
    ylabel('$z$ [mm]');
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);
    box on;
    grid on; grid minor;
    cb=colormap(jet);
        cb=colorbar('color','white','FontSize',12,'TickLabelInterpreter','latex','Fontsize',12,'Location','south');
        cb.Position         =[0.2500 0.3767 0.1500 0.1000];
        cb.Label.String     ='$\log_{10}(1+\epsilon^p_{\mathrm{eqv}})$ [-]';
        cb.Label.Units      ='normalized';
        cb.Label.Position   =[1.8 0.9];
        cb.Label.FontSize   =12;
        cb.Label.Interpreter='Latex';
        caxis([0 1]);
    drawnow
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
name = ['figElastoPlasticCollapseEquivStr2DZoomedIn' num2str(mpD.n)];
print(gcf,name,'-depsc');
print(gcf,name,'-dpng');



























%% 
fs=load('Data_Buietal_2008_experimental_failure_surface.txt');
fs= [fs(:,1) fs(:,2)];
S =load('Data_Buietal_2008_experimental_surface.txt');
S= [S(:,1) S(:,2)];

np = 1:mpD.n;
disp = sqrt(u(:,1).^2+u(:,2).^2);;
minu = 5e-4;
disp(disp>minu)=1;
disp(disp<minu)=0;
minu = minu/1e-3;
fig3=figure(3);
set(fig3,'Units','pixels','Position',[85 604 560 248]);
hold on
ax3=scatter(x(np,1)*1000,x(np,2)*1000,4.0,(disp(np)),'filled');
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
cb.TickLabels = {['$u_p \leq ',num2str(minu),'$ mm'],['$u_p>',num2str(minu),'$ mm']};
        cb.Title.String     =['ep2De v1.0, $n_{mp}=',num2str(mpD.n),'$'];
        cb.Title.Interpreter='Latex';
ax2=plot(fs(:,1)*1000,fs(:,2)*1000,'b:','LineWidth',2);
ax1=plot(S(:,1)*1000,S(:,2)*1000,'b--','LineWidth',2);
box on;
grid on;
tit = {'Experiment: final geometry','Experiment: failure surface'};
h1=legend([ax1 ax2 ax3],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.4407 0.5033 0.4514 0.1730],'NumColumns',1);
hold off;
axis equal;
xlabel('$x$ [mm]');
ylabel('$z$ [mm]');
xlim([0 0.5*1000]);
    yticks([0 0.05*1000 0.1*1000]);
    ylim([0 lz*1000]);
set(gca,'FontSize',15,'TickLabelInterpreter','latex');

name = ['figElastoPlasticCollapse2D_nmp_' num2str(mpD.n)]
print(gcf,name,'-depsc');
print(gcf,name,'-dpng');
























