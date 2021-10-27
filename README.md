
# ep2-3De v1.1 #

  

## Software

  

****ep2-3De v1.1**** is a software suite implementation of the material point method (e.g., GIMP) to model elasto-plastic problems on a GPU, using a Drucker-Prager rheology, considering either 2D or 3D configuration. It is written in CUDA C, visualisation and input data generations are written in MATLAB. The latter generates data needed, compile source codes and runs the executables. Post-processing routines can be used for output visualization. The new version v1.1 comes with a multi-GPU implementation of the initial release v1.0.

As a disclaimer, the interested user should know that we did not test the MPI implementation on other GPU clusters. Hence, we do not guarantee the multi-GPU implementation to be stable or even working on other computing systems. 

  

## Distributed software, directory content

- single_GPU/
    - ep2De_v1.0
        - scripts
            - setting_Exp_1b
        - src

    - ep3De_v1.0

        - scripts

            - setting_Exp1a

            - setting_Exp2a

            - setting_Exp2b

        - src

        - bash_ep3De_v1.0
- multi_GPU/
    - rkf
    - script
        - setting_Exp2b_performance_Titan
        - setting_Exp2b_performance_Volta
    - src

  

## Usage example (single GPU)

To execute the GPU code, you need an Nvidia GPU. Older architectures (i.e., Kepler or Maxwell) mights results in an odd behaviour of gpu codes. We recommand the usage of the most recent GPU architectures. The GPU architecture ````--gpu-architecture=sm_XX```` can be configured in the MATLAB routine ````compileExpXX````, which is located in subfolder ````./scripts/setting_ExpXX````.

  

ep2-3De v.1.0 was compiled and executed on the following architecture:

- Volta architecture, i.e., V100 with ````sm_70````,

- Turing architecture, i.e., GTX1650 and RTX2080ti with ````sm_75````,

- Ampere architecture, i.e., A100 and RTX3090 with ````sm_80````.

  

As an example, to run Exp2b, you can:

- Select the 3D solver folder (ep3De_v1.0) and navigate in ````./ep3De_v1.0/scripts/setting_Exp2b/````

- Open the MATLAB routine ````compileExp2b.m```` and:

- a) select the GPU architecture installed on the computer, i.e., ````--gpu_architecture=sm_xx````,

- b) define the number of threads per block ````BCK1```` and ````BCK2````, i.e., ````BCK1<=16```` and ````BCK2<=1024````

- c) set the precision arithmetic ````typeD = 'float'```` or ````typeD = 'double'````,

- d) (opt.) define ````FPS```` (frame per seconds) and ````D```` (local damping).

- Open the MATLAB routine Exp2b.m and run it.

  

## Miscellaneous

  

- For GPUs with compute capability lower than 6.0 (i.e., Kepler or Maxwell architecture), ````atomicAdd```` only supports floating point arithmetic. Shall ep2-3De_v1.0 compiled on a GPU with a lower compute capability, the arithmetic precision should be set to ````float````.

  

## Reference

Please cite us if you use our routines (will be updated upon final publication):

```tex

@article{Wyser2021,

title={An explicit GPU-based material point method solver for

elastoplastic problems (ep2-3De v1.0)},

author={Wyser, Emmanuel and Alkhimenkov, Yury and Jaboyedoff, Michel and Podladchikov, Yury Y.},

year={2021},

}

```

  

### License

Copyright (C) 2021  Emmanuel Wyser, Yury Alkhimenkov, Michel Jaboyedoff, Yury Y. Podladchikov.

  

It is free software: you can redistribute it and/or modify

it under the terms of the GNU General Public License as published by

the Free Software Foundation, either version 3 of the License, or

(at your option) any later version.

  

ep2-3De v1.0 is distributed in the hope that it will be useful,

but WITHOUT ANY WARRANTY; without even the implied warranty of

MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

GNU General Public License for more details.

  

You should have received a copy of the GNU General Public License

along with ep2-3De v1.0.  If not, see <http://www.gnu.org/licenses/>.

  

### Contact: manuwyser@gmail.com
