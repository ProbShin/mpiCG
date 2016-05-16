#########################################################################
# File Name: run.sh
# Author: Shin Cheng
# Descriptions: 
# Created Time: Fri May 13 15:43:52 2016
#########################################################################
#!/bin/bash

#mpirun -np 8 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./gm ./tinyA_dd.mtx ./tinyb_dd.mtx ./mysolution.mtx
mpirun -np 8 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./a.out ./tinyA_dd.mtx ./tinyb_dd.mtx ./mysolution.mtx
