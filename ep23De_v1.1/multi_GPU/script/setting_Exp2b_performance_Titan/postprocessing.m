clear all

GPU = load('Performance_oneGPU.mat');
np             = GPU.Perf(1,:)
perf_mpi(1,:) = [1,GPU.Perf(2,:)]

it = 1
for nprocess = [2 4 8]
    it = it+1
nGPUs = load(['Performance_',num2str(nprocess),'GPUs.mat']);

if(nprocess<8)
perf_mpi(it,:) = [nprocess,nGPUs.Perf(2,:)]
else
perf_mpi(it,:) = [nprocess,nGPUs.Perf(2,1:8)]    
end

end

fig1=figure(1);
set(fig1,'Units','pixels','Position',[113 592 559 325]);
hold on
ax1 = plot(np(1:end-1),perf_mpi(1,2:end-1),'k-*')
ax2 = plot(np(1:end-1),perf_mpi(2,2:end-1),'-o')
ax3 = plot(np(1:end-1),perf_mpi(3,2:end-1),'-s')
ax4 = plot(np(1:end-1),perf_mpi(4,2:end-1),'-x')
hold off
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('$n_{mp}$ [-]')
ylabel('Wall-clock time [s]')
yticks([10 60 120])
yticklabels({'10','60','120'})
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
tit = {['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(1,1)),'$'],...
       ['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(2,1)),'$'],...
       ['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(3,1)),'$'],...
       ['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(4,1)),'$']};
h1=legend([ax1 ax2 ax3 ax4],tit);
title(h1,'Tesla V100 NVlink 32 GB');
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.1487 0.6929 0.4134 0.2034],'NumColumns',2);
box on






e = 1./[2;4;8].*(perf_mpi(1,2:end-1)./perf_mpi(2:end,2:end-1))
figure
plot(np(1:end-1),e)




fig2=figure(2);
set(fig2,'Units','pixels','Position',[113 592 559 325]);
hold on
ax1 = plot(np(1:end-1),perf_mpi(1,2:end-1)./perf_mpi(1,2:end-1),'k-*')
ax2 = plot(np(1:end-1),perf_mpi(1,2:end-1)./perf_mpi(2,2:end-1),'-o')
ax3 = plot(np(1:end-1),perf_mpi(1,2:end-1)./perf_mpi(3,2:end-1),'-s')
ax4 = plot(np(1:end-1),perf_mpi(1,2:end-1)./perf_mpi(4,2:end-1),'-x')
hold off
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('$n_{mp}$ [-]')
ylabel('Speed-up')
yticks([1 2 4 8])
yticklabels({'1x','2x','4x','8x'})
ylim([0 10])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
tit = {['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(1,1)),'$'],...
       ['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(2,1)),'$'],...
       ['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(3,1)),'$'],...
       ['$n_{\mathrm{GPU}} = ',num2str(perf_mpi(4,1)),'$']};
h1=legend([ax1 ax2 ax3 ax4],tit);
title(h1,'Tesla V100 NVlink 32 GB');
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.1487 0.6929 0.4134 0.2034],'NumColumns',2);
box on
% writematrix([x p0],'myData.txt','Delimiter',';'); 
% type myData.txt;






