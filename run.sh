#########################################################################
# File Name: run.sh
# Author: Shin Cheng
# Descriptions: 
# Created Time: Fri May 13 15:43:52 2016
#########################################################################
#!/bin/bash


#mpirun -np 8 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./tinyA.mtx ./tinyb.mtx ./mysolution.mtx 

#mpirun -np 8 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./bcsstk18.mtx ./b11948.mtx ./mysolution.mtx 


#mpirun -np 4  -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./LFAT5/LFAT5.mtx ./b_LFAT.mtx ./mysolution.mtx 

#mpirun -np 8  -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./A_32.mtx ./b_32.mtx ./mysolution.mtx 


#mpirun -np 8  -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./A_10000.mtx ./b_10000.mtx ./mysolution.mtx 


mpirun -np 8 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ../../flan/flan_orig/flan.mtx ./b_flan.mtx ./mysolution.mtx 


