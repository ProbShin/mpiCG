#########################################################################
# File Name: run.sh
# Author: Shin Cheng
# Descriptions: 
# Created Time: Fri May 13 15:43:52 2016
#########################################################################
#!/bin/bash


#mpirun -np 8 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./tinyA.mtx ./tinyb.mtx ./mysolution.mtx 

#mpirun -np 20  -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../cfd2/cfd2.mtx ./b_cfd2.mtx ./cfd_sol.mtx 

#mpirun -np 1  -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../cfd2/cfd2_order.mtx ./b_cfd2_order.mtx ./cfd_sol.mtx 

#mpirun -np 8 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../Emilia/Emilia_923.mtx ./b_Emilia.mtx ./Emilia_sol.mtx 

#mpirun -np 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../Emilia/Emilia_923.mtx ./b_Emilia.mtx ./Emilia_sol.mtx 

#mpirun -np 15 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../Geo/Geo_1438.mtx ./b_Geo.mtx ./Geo_sol.mtx 

#mpirun -np 1 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../cfd2/cfd2_order.mtx ./b_cfd.mtx ./cfd_sol.mtx 
#mpirun -np 8 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../stockf/stockf.mtx ./b_stockf.mtx ./stockf_sol.mtx 

mpirun -np 1  -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 1 ./mycg ../../bcsstk18/bcsstk18_aug16.mtx ../../CG_advance/b_bcsstk18_order.mtx ./sol.mtx

#mpirun -np 4  -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./LFAT5/LFAT5.mtx ./b_LFAT.mtx ./mysolution.mtx 

#mpirun -np 8  -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ./A_32.mtx ./b_32.mtx ./mysolution.mtx 




#mpirun -np 8 -genv I_MPI_DEVICE=ssm  -genv OMP_NUM_THREADS 8 ./mycg ../../flan/flan_orig/flan.mtx ./b_flan.mtx ./mysolution.mtx 


