#!/bin/bash
# shellscript.sh bash has to be located in ./script/setting_Exp2b/
#------------------------------------------------------------------
# load what should be loaded
module load cuda/11.0
module load openmpi/gcc54-200c
#------------------------------------------------------------------
# remove what should be removed
rm *.dat
rm *.out
#------------------------------------------------------------------
#compile GPU code
nvcc ../../src/gpu_main.cu -o gpu.out -arch=sm_70 -O3 -DSIM=$1 -DDAT=$2 -DFPS=$3 -DD=$4 -DBCK1=$5 -DBCK2=$6

