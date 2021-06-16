%% MATLAB SCRIPT TO RUN COMPILATION FOR GPU AND CPU EXECUTABLES (.EXE)
% delete what needs to be deleted
delete('*.dat','*.out','*.exe','*.avi','*.lib','*.exp');
% set precision arithmetic
typeD = 'double';
BCK1  = 32;
BCK2  = 256;
% compile gpu code
if(ismac || isunix)
    call = ['nvcc ../../src/DP_gpu_main.cu -o gpu.out --gpu-architecture=sm_70 -w -O3 '];
elseif(ispc)
    call = ['nvcc ../../src/DP_gpu_main.cu -o gpu.exe --gpu-architecture=sm_70 -w -O3 '];
else
    disp('Outer space OS');
end
if(typeD=='single')
    arg  = ['-DDAT=float  -DSIM=',int2str(sim),' -DFPS=1 -DD=0.05 -DBCK1=',num2str(BCK1),' -DBCK2=',num2str(BCK2),''];
    system([call arg]); % trouble ? type nvcc --help in terminal
elseif(typeD=='double')
    arg  = ['-DDAT=double -DSIM=',int2str(sim),' -DFPS=1 -DD=0.05 -DBCK1=',num2str(BCK1),' -DBCK2=',num2str(BCK2),''];
    system([call arg]); % trouble ? type nvcc --help in terminal
end
