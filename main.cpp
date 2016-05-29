
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
#define FILEA_READER_ID  ROOT_ID
#define FILEB_READER_ID  ROOT_ID

#define MAX_ITER 100000 


#define ACRED     "\x1b[31m"
#define ACGREEN   "\x1b[32m"
#define ACYELLOW  "\x1b[33m"
#define ACBLUE    "\x1b[34m"
#define ACMAGENTA "\x1b[35m"
#define ACCYAN    "\x1b[36m"
#define ACRESET   "\x1b[0m"



bool check_head_sparse(FILE* fp, int &m, int &n, int&nnz );
bool check_head_dense (FILE* fp, int &m, int &n          );
void check_head_show  (FILE* fp);


//============Solve Ax=b using parallel CG========================//
void myCG(int n, double* vrr, int* colind, int* rbegin, double *b, double*xout,  int rank, int nproc ){
  int np = n/nproc;
  int i,j,k,l;
  int counter=0;
  double j_total, j_start, j_end;
  j_total = 0;
  double *p_loc =  new double[np];
  double *p = new double[n ];
  
  
  double *x = new double[np]; 
  double *r = new double[np];

  double *ap= new double[np];
  double rnorm = 0;
  double pap   = 0;
  double alpha = 0;
  double beta  = 0;
  double tmp   = 0;
  double maxnorm=0;
  int offset = rank*np;

  for(int i=0; i<np; i++){
    x[i]=0; ap[i]=0;
    r[i]=b[i];
    p_loc[i]=r[i];
    rnorm+=r[i]*r[i];
    if(fabs(r[i])>maxnorm) {
      maxnorm = fabs(r[i]);
    }
  }

  MPI_Allgather(p_loc, np, MPI_DOUBLE, p, np, MPI_DOUBLE, MPI_COMM_WORLD);

  double xc=rnorm; 
  rnorm=0;
  MPI_Allreduce(&xc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  xc = maxnorm;
  MPI_Allreduce(&xc, &maxnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  int iter=0;
  while((maxnorm>1e-6) && (iter<MAX_ITER)){
    pap = 0;
    //mat vec
    j_start = omp_get_wtime();
    for(i=0; i<np; i++){
      ap[i]=0;
      for(j=rbegin[i+offset]; j<rbegin[i+1+offset]; j++){
        ap[i] += vrr[j]*p[colind[j]];
        counter++;
      }
      pap+=p[offset+i]*ap[i];
    }
    xc = pap;
    
    j_end = omp_get_wtime();
    j_total+= j_end- j_start;
    
    MPI_Allreduce(&xc, &pap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    alpha = rnorm/pap;
    beta  = rnorm; 
    rnorm = 0;
    maxnorm = 0;
    for(i=0; i<np; i++) {
      x[i]=x[i]+alpha*p[offset+i];
      r[i]=r[i]-alpha*ap[i];
      rnorm += r[i]*r[i];
      if( (fabs(r[i])>maxnorm)){
        maxnorm = fabs(r[i]);
      }
    }
    xc=rnorm;
    MPI_Allreduce(&xc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    xc=maxnorm;
    MPI_Allreduce(&xc, &maxnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    beta = rnorm/beta;
    for(i=0;i<np; i++){
      p_loc[i]=r[i]+beta*p_loc[i];
    }
    
    MPI_Allgather(p_loc , np, MPI_DOUBLE, p, np, MPI_DOUBLE, MPI_COMM_WORLD);
    iter ++;
  }//END OF WHILE
  printf("Process %d  counter=%d, iter=%d, maxnorm=%.01f\n", rank, counter, iter, maxnorm);

  MPI_Gather(x, np, MPI_DOUBLE, xout, np, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);



  printf("total matvec time %f\n", j_total);

  delete []x;
  delete []r;
  delete []p;
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

  int    *colind=NULL; //new int [annz];
  int    *rbegin=NULL; //new int [m+1];
  double *avrr  =NULL;  
  
  double *brr   =NULL;
  double *b     =NULL;              //local b
  
  double *xrr   =NULL;    //*&xrr   = brr;
  int m, annz, mrhs;
  int mp;                 //m/p -> mm

  FILE *f1 = NULL;
  FILE *f2 = NULL;

  if(rank==FILEA_READER_ID)
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

  }

  if(rank==FILEB_READER_ID)
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

    if(rank==FILEB_READER_ID) {
      for(int i=0; i<mrhs; i++) {
        fscanf(f2, "%lf\n", brr+i);
      }
    }
  }


  MPI_Bcast(&m,    1, MPI_INT, FILEA_READER_ID, MPI_COMM_WORLD);
  MPI_Bcast(&mrhs, 1, MPI_INT, FILEB_READER_ID, MPI_COMM_WORLD);
  MPI_Bcast(&annz, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);

  mp = m/nproc;
  if(m!=mrhs) exit(3);

  if(avrr==NULL  )  avrr = new double[annz];
  if(   b==NULL  )     b = new double[mp  ];
  if(colind==NULL) colind= new int [annz];
  if(rbegin==NULL) rbegin= new int [m+1 ];
  if(xrr ==NULL  )    xrr= new double[m];


  MPI_Bcast  (avrr,   annz, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
  MPI_Bcast  (colind, annz, MPI_INT,    ROOT_ID, MPI_COMM_WORLD);
  MPI_Bcast  (rbegin, m+1 , MPI_INT,    ROOT_ID, MPI_COMM_WORLD);
  
  MPI_Scatter(brr, mp, MPI_DOUBLE, b, mp, MPI_DOUBLE, FILEB_READER_ID, MPI_COMM_WORLD);

  if(rank==ROOT_ID){
    rtime = omp_get_wtime();
    fprintf(stdout,"start count\n");
  }


  myCG(m, avrr, colind, rbegin, b, xrr, rank, nproc);

  if(rank==ROOT_ID){
    fprintf(stdout,"time %f\n",-(rtime-=omp_get_wtime()));
  }

  //if(rank==ROOT_ID){
  //  printf("the result are:\n");
  //  for(int i=0; i<m; i++)
  //    fprintf(stdout, "%f ", xrr[i]);
  //  printf("\n"); 
  //}


  MPI_Finalize();

  delete[] colind;
  delete[] rbegin;
  delete[] avrr;
  delete[] brr;

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













