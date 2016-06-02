#########################################################################
# File Name: run.sh
# Author: Shin Cheng
# Descriptions: 
# Created Time: Fri May 13 15:43:52 2016
#########################################################################
#!/bin/bash


#mpirun -np 8 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./tinyA.mtx ./tinyb.mtx ./mysolution.mtx 

#mpirun -np 8 -f ./myhost -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../cfd2/cfd2_order.mtx ./b_cfd.mtx ./cfd_sol.mtx 

mpirun -np 48 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../Emilia/Emilia_923_order.mtx ./b_Emilia_order.mtx ./Emilia_sol.mtx 

#mpirun -np 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../cfd2/cfd2_order.mtx ./b_cfd.mtx ./cfd_sol.mtx 
#mpirun -np 8 -f ./myhost -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../stockf/stockf.mtx ./b_stockf.mtx ./stockf_sol.mtx 

#mpirun -np 4 -f ./myhost -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ./misc/bcsstk18.mtx ./misc/b_bcsstk18.mtx ./sol.mtx

#mpirun -np 4  -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./LFAT5/LFAT5.mtx ./b_LFAT.mtx ./mysolution.mtx 

#mpirun -np 8  -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./A_32.mtx ./b_32.mtx ./mysolution.mtx 


#mpirun -np 8  -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./A_10000.mtx ./b_10000.mtx ./mysolution.mtx 


#mpirun -np 8 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ../../flan/flan_orig/flan.mtx ./b_flan.mtx ./mysolution.mtx 


#mpirun -np 2 -f ./myhost -perhost 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./LFAT5.mtx ./b_LFAT5.mtx ./mysolution.mtx 
