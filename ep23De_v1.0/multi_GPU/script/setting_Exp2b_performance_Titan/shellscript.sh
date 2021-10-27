#!/bin/bash
#------------------------------------------------------------------
# initialization
module load cuda/11.4
module load openmpi/215_g65_c10

rm *.out

mpirun=$(which mpirun)
if [ $# -lt 1 ]; then
    echo $0: usage: shellscript.sh nprocs
    exit 1
fi
nprocs=$7
#------------------------------------------------------------------
# compile GPU source code
nvcc ../../src/gpu_main.cu -o gpu.out -arch=sm_52 --compiler-options -O3 -w --compiler-bindir mpic++ -DSIM=$1 -DDAT=$2 -DFPS=$3 -DD=$4 -DBCK1=$5 -DBCK2=$6 -DDIMS_X=$7 -DBFS=$8
#------------------------------------------------------------------
# run GPU code 
run_cmd="-np $nprocs -rf ../../rkf/gpu_rankfile_128 ./gpu.out"
echo $mpirun $run_cmd
$mpirun $run_cmd
#------------------------------------------------------------------
## nvidia profiler
#$mpirun $run_cmd $nvprof --profile-from-start off --export-profile mpm_3D.%q{OMPI_COMM_WORLD_RANK}.prof -f ./gpu.out