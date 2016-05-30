
#include<iostream>
#include<cstdio>     /* stderr */
#include<cstdlib>    /* exit */
#include<fstream>    /* fstream */

#include <unistd.h>  /* sleep */
#include <cmath>     /* fabs */
#include "mpi.h"
#include "mmio.h"    /* matrix market */
#include "mkl.h"
#include "omp.h"

#include <vector>
#include <map>

using namespace std;


#define ROOT_ID   0
#define FILE_A_READER_ID  ROOT_ID
#define FILE_B_READER_ID  ROOT_ID

#define MAX_ITER 100000 


#define ACRED     "\x1b[31m"
#define ACGREEN   "\x1b[32m"
#define ACYELLOW  "\x1b[33m"
#define ACBLUE    "\x1b[34m"
#define ACMAGENTA "\x1b[35m"
#define ACCYAN    "\x1b[36m"
#define ACRESET   "\x1b[0m"

#define cxTMP 2

bool check_head_sparse(FILE* fp, int &m, int &n, int&nnz );
bool check_head_dense (FILE* fp, int &m, int &n          );
void check_head_show  (FILE* fp);


//============Solve Ax=b using parallel CG========================//
void myCG(int n, double* vrr, int* colind, int* rbegin, double *b_loc, double*x_loc,  int rank, int nproc ){
  double t1,t2,t3,t4,t5;
  
  int np = n/nproc;
  int i,j,k,l;
  int counter=0;
  double t_total, t_time;
  t_total = 0;
  double *p_loc =  new double[np];
  double *p_glb =  new double[n];
  double *&x = x_loc;
  double *&r = b_loc; 
  double *ap= new double[np];
  double r2norm_loc = 0, r2norm_glb = 0;
  double pap_loc    = 0, pap_glb    = 0;
  double alpha = 0;
  double beta  = 0;
  double tmp   = 0;
  int offset = rank*np;
  double b2norm = 0; 


  for(int i=0; i<np; i++){
    ap[i]=0;
    //r[i]=b_loc[i];
    p_loc[i]=r[i];
    r2norm_loc+=r[i]*r[i];
    x[i]=0;
  }


  MPI_Allreduce(&r2norm_loc, &r2norm_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgather(p_loc, np, MPI_DOUBLE, p_glb, np, MPI_DOUBLE, MPI_COMM_WORLD);  

  b2norm=r2norm_glb;
  int iter=0;
  while((r2norm_glb>b2norm*1e-12) && (iter<MAX_ITER)){
    pap_loc = 0;
    //mat vec
    //if(rank==0 && iter==cxTMP){
    //  printf("before MATVEC:\n");
    //  for(int i=0;i<n; i++)
    //    printf("p[%d]=%.10f\n",i,p_glb[i]);
    //}
    
    t_time = omp_get_wtime();
    for(int i=0; i<np; i++){
      ap[i]=0;
      for(j=rbegin[i+offset]; j<rbegin[i+1+offset]; j++){
        ap[i] += vrr[j]*p_glb[colind[j]];   //ap[i] += vrr[j]*p[colind[j]];
        counter++;
      }
      pap_loc+=p_loc[i]*ap[i];   //pap_loc+=p[offset+i]*ap[i];
    }
    
    t_total+= -(t_time-=omp_get_wtime());
    //if(rank==0 && iter==cxTMP){
    //  printf("after MATVEC:\nrank%d iter%dAP:\n",rank,iter);
    //  for(int i=0;i<np; i++)
    //    printf("AP[%d]_iter%d=%f\n",i+offset, iter ,ap[i]);
    // 
    //}

    MPI_Allreduce(&pap_loc, &pap_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //if(rank==0 && iter==cxTMP){
    //  printf("\npAp_iter%d=%f\n",iter,pap_glb);
    //}
  
    //if(rank==0 && iter==cxTMP){
    //  printf("\nres_old_iter%d=%f\n", iter, r2norm_glb);
    //}

    alpha = r2norm_glb/pap_glb;
    //if(rank==0 && iter==cxTMP){
    //  printf("\nalpha_iter%d=%.30f\n",iter,alpha);
    //}
    
    beta  = r2norm_glb; 
    r2norm_loc = 0;
    for(int i=0; i<np; i++) {
      x[i]=x[i]+alpha*p_loc[i]; 
      r[i]=r[i]-alpha*ap[i];
      r2norm_loc += r[i]*r[i];
    }
    
    //if(rank==0 && iter==cxTMP) {
    //  printf("rank%dxnew_iter%d:\n",rank,iter);
    //  for(int i=0;i<np; i++)
    //    printf("xnew[%d]=%lg\n",i, x[i]);
    //}
    
    
    MPI_Allreduce(&r2norm_loc, &r2norm_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    beta = r2norm_glb/beta;
    
    //if(rank==0 && iter==cxTMP){
    //  printf("\nbeta_iter%d=%f\n",iter,beta);
    //}
    


    //if(rank==0 && iter==cxTMP){
    //  printf("\nBEFOR p=r+beta*p:\n");
    //  printf("p:\n");
    //  for(int i=0;i<np; i++){
    //    printf("%lg\n", p_loc[i]);
    //  }
    //  printf("r:\n");
    //  for(int i=0;i<np; i++){
    //    printf("%lg\n", r[i]);
    //  }
    //}


    for(int i=0;i<np; i++){
      p_loc[i]=r[i]+beta*p_loc[i];
    }
    
    //if(rank==0 && iter==cxTMP){
    //  printf("\nAFTER p=r+beta*p:\n");
    //  printf("pnew:\n");
    //  for(int i=0;i<np; i++)
    //    printf("%lg\n", p_loc[i]);
    //}



    MPI_Allgather(p_loc, np, MPI_DOUBLE, p_glb, np, MPI_DOUBLE, MPI_COMM_WORLD);  
    
    //if(rank==0 && iter==cxTMP){
    //  printf("\nrank%d pnew_iter%d:\n",rank, iter);
    //  for(int i=0;i<n; i++)
    //    printf("pnew[%d]_iter%d=%.10f\n",i, iter, p_glb[i] );
    //}
      
    iter ++;
  }//END OF WHILE
  
  //MPI_Gather(x, np, MPI_DOUBLE, xout, np, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

  //printf("rank%d t1=%f t2=%f t3=%f t4=%f t5=%f\n",rank, -t1,-t2,-t3,-t4,-t5);
  //printf("Process %d  counter=%d, iter=%d(max%d), fnorm**2=%.01f\n", rank, counter, iter, MAX_ITER, rfnorm_glb);
  printf("Process %d  counter=%d, iter=%d(max%d) r2norm**2=%f\n", rank, counter, iter, MAX_ITER, r2norm_glb);
  printf("total matvec time %f\n", t_total);
  //delete []r;
  //delete []p;
  delete []p_glb;
  delete []ap;
  return;
}// END OF myCG




// =========================================
// solve Ax=b with CG 
// A input
// b input
// b will be overwrite by x
// =========================================
int main(int argc, char const *argv[]){
  int nproc, rank;
  double rtime;
  double j_start, j_end, j_total = 0.0;
  MPI_Init(&argc, (char***) &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  

  if( argc!=4 ){
    if(rank==ROOT_ID)
      fprintf(stderr,  ACRED"usage: ""%s "ACGREEN"<inA.mtx> <inb.mtx>"ACRED" <outx.mtx>"ACRESET"\n", argv[0]);
    MPI_Finalize();
    exit(1);

  }

  int    *colind=NULL; //new int [annz];   for csr
  int    *rbegin=NULL; //new int [m+1];    for csr
  double *avrr  =NULL; //new double[annz]; for csr
  
  double *brr   =NULL; //new double[m];    b
  double *b_loc =NULL; //new double[m/p]   local b
  
  double *&xrr  =brr ; //*&xrr = brr ;     x
  double *x_loc =NULL; //*&xrr = brr ;     x
  int m, annz, mrhs;   //mtx size; actual nnz; number of right hand sides
  int mp;              //m/p -> mm

  double io_time=0;
  FILE *f1 = NULL;
  FILE *f2 = NULL;


  if(rank==ROOT_ID)
    io_time = omp_get_wtime();

  if(rank==FILE_A_READER_ID)
  {
    int tn, tnnz;
    f1 = fopen(argv[1], "r");
    if(!f1||!check_head_sparse(f1, m, tn, tnnz)){
      printf("cannot check head!");
      printf("%d %d %d\n", m,tn, tnnz);
    }

    if( m!=tn || tnnz==0 ) {
      fprintf(stderr, "read "ACRED"%s!"ACRESET" failed\n", argv[1]);
      exit(2);
    }
    
    //read from symm to adjacency matrix
    vector< vector<int> >      nodelist(m, vector<int>());
    vector< vector<double> >   valulist(m, vector<double>());

    int r,c;
    double v;
    annz=0;
    for(int i=0; i<tnnz; i++){
      fscanf(f1, "%d %d %lf\n", &r, &c, &v);
      r--; c--;
      if(r==c){
        nodelist[c].push_back(r);
        valulist[c].push_back(v);
        annz+=1;
      }
      else{ //r >= b
        nodelist[r].push_back(c);
        nodelist[c].push_back(r);

        valulist[r].push_back(v);
        valulist[c].push_back(v);

        annz+=2;
      }
    }

    //now construct the graph
    avrr   = new double[annz];
    colind = new int   [annz];
    rbegin = new int   [m+1 ];

    int cnt=0;
    rbegin[0]=0;
    for(int i=0; i<m; i++){
      rbegin[i+1] = rbegin[i] + nodelist[i].size();
      for(int j=0; j<nodelist[i].size(); j++){
        colind[cnt] = nodelist[i][j]; 
        avrr[cnt] = valulist[i][j];
        cnt++;
      }
    }

  } //end of file a read

  if(rank==FILE_B_READER_ID)
  {
    f2 = fopen(argv[2], "r");
    int tn;
    if(!f1       ||
        !check_head_dense(f2, mrhs, tn) ||
        tn != 1) {
      fprintf(stderr, "read "ACRED"%s!"ACRESET" failed\n", argv[2]);
      exit(2);
    }
    brr  = new double[mrhs]; 

    for(int i=0; i<mrhs; i++) {
      fscanf(f2, "%lf\n", brr+i);
    }
  }


  MPI_Bcast(&m,    1, MPI_INT, FILE_A_READER_ID, MPI_COMM_WORLD);
  MPI_Bcast(&mrhs, 1, MPI_INT, FILE_B_READER_ID, MPI_COMM_WORLD);
  MPI_Bcast(&annz, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);

  mp = m/nproc;
  if(m!=mrhs) exit(3);

  if(avrr  ==NULL ) avrr  = new double[annz];
  if(colind==NULL ) colind= new int [annz];
  if(rbegin==NULL ) rbegin= new int [m+1 ];
  
  if(b_loc ==NULL ) b_loc = new double[mp ];
  if(x_loc ==NULL ) x_loc = new double[mp ];

  //broadcast matrix A
  MPI_Bcast  (avrr,   annz, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
  MPI_Bcast  (colind, annz, MPI_INT,    ROOT_ID, MPI_COMM_WORLD);
  MPI_Bcast  (rbegin, m+1 , MPI_INT,    ROOT_ID, MPI_COMM_WORLD);

  //scatter right hand side b
  MPI_Scatter(brr, mp, MPI_DOUBLE, b_loc, mp, MPI_DOUBLE, FILE_B_READER_ID, MPI_COMM_WORLD);
  
  fprintf(stdout,"rank%d arrived\n",rank);
  if(rank==ROOT_ID)
    fprintf(stdout,"file_io(preparation) takes %lg(s)\n", -(io_time-=omp_get_wtime()));
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==ROOT_ID)
    fprintf(stdout,"start count, lineID%0.0f\n", rtime=omp_get_wtime() );

  myCG(m, avrr, colind, rbegin, b_loc, x_loc, rank, nproc);

  if(rank==ROOT_ID)
    fprintf(stdout,"time "ACRED"%f"ACRESET"\n",-(rtime-=omp_get_wtime()));

  //if(rank==ROOT_ID){
  //  printf("the result are:\n");
  //  for(int i=0; i<m; i++)
  //    fprintf(stdout, "%f ", xrr[i]);
  //  printf("\n"); 
  //}

  if(rank==0){ 
    printf("x of rank%d\n",rank);
    for(int i=0;i<20;i++)
      printf("%f\n",x_loc[i]);
    printf("\n");
  }

  MPI_Finalize();

  //delete[] colind;
  //delete[] rbegin;
  //delete[] avrr;
  //delete[] brr;
  
}




bool check_head_sparse(FILE* fp, int &m, int &n, int&nnz ){
  MM_typecode matcode;
  //check_head_show(fp);
  if( mm_read_banner(fp, &matcode)!=0 || 
      mm_is_matrix(matcode)!=1        ||
      mm_is_sparse(matcode)!=1        ||
      mm_is_coordinate(matcode)!=1    ||
      mm_is_real(matcode)!=1          ||
      mm_is_symmetric(matcode)!=1     ||
      //mm_is_general(matcode)!=1       ||
      mm_read_mtx_crd_size(fp, &m, &n, &nnz)!=0 ) 
    return false;
  return true;
}


bool check_head_dense (FILE* fp, int &m, int &n          ){
  MM_typecode matcode;
  if( mm_read_banner(fp, &matcode)!=0 || 
      mm_is_matrix(matcode)!=1        ||
      mm_is_dense(matcode)!=1         ||
      mm_is_array(matcode)!=1         ||
      mm_is_real(matcode)!=1          ||
      mm_is_general(matcode)!=1       ||
      mm_read_mtx_array_size(fp, &m, &n)!=0 ) 
    return false;
  return true;
}



void check_head_show(FILE* fp){
  MM_typecode matcode;
  if(mm_read_banner(fp, &matcode) != 0) {
    printf("Could not process Matrix Market banner.\n");
    return;
  }
  if( mm_is_matrix(matcode) ) printf("is matrix\n");
  if( mm_is_sparse(matcode) ) printf("is sparse\n");
  if( mm_is_coordinate(matcode) ) printf("is coordinate\n");
  if( mm_is_dense(matcode) ) printf("is dense\n");
  if( mm_is_array(matcode) ) printf("is array\n");
  if( mm_is_complex(matcode) ) printf("is complex\n");
  if( mm_is_real(matcode) ) printf("is real\n");
  if( mm_is_pattern(matcode) ) printf("is pattern\n");
  if( mm_is_integer(matcode) ) printf("is integer\n");
  if( mm_is_symmetric(matcode) ) printf("is symm\n");
  if( mm_is_general(matcode) ) printf("is general\n");
  if( mm_is_skew(matcode) ) printf("is skew\n");
  if( mm_is_hermitian(matcode) ) printf("is hermitian\n");
}













