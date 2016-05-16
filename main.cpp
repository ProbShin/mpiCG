/*************************************************************************
    > File Name: main.cpp
    > Author: Shin Cheng
    > Descriptions: 
    > Created Time: Mon May 16 12:22:57 2016
 ************************************************************************/

#include<iostream>
#include<cstdio>     /* stderr */
#include<cstdlib>    /* exit */

using namespace std;



#define ACRED     "\x1b[31m"
#define ACGREEN   "\x1b[32m"
#define ACYELLOW  "\x1b[33m"
#define ACBLUE    "\x1b[34m"
#define ACMAGENTA "\x1b[35m"
#define ACCYAN    "\x1b[36m"
#define ACRESET   "\x1b[0m"


// =========================================
// solve Ax=b with CG 
// A input
// b input
// b will be overwrite by x
// =========================================
int main(int argc, char const *argv[]){
  if( argc!=3 ){
    fprintf(stderr,  ACRED"usage: ""%s "ACGREEN"<inA.mtx> <inb.mtx>"ACRED" <outx.mtx>"ACRESET"\n", argv[0]);
    exit(1);
  }


  int    *airr, *ajrr;
  double *avrr;  
  double *brr;
  double *&xrr = brr;
  int m, p;
  const int base =1;      //pain in ass mm is 1 based
  fstream f1(argv[1],"r");
  fstream f2(argv[2],"r");
  {
    int tn, int tnnz;
    if(!f1 ||
       check_head_sparse(f1, &m, &tn, &tnnz) ||
       m!=tn ||
       tnnz==0 
      ) {
      fprintf(stderr, "read "ACGRED"%s!"ACRESET" failed\n", argv[1]);
      exit(2);
    }
    airr = new double[tnnz];
    ajrr = new double[tnnz];
    avrr = new double[tnnz];
  
  }
  {
    int tm, tn;
    if(!f1  ||
        check_head_dense(f2, &tm, &tn) ||
        tn  !=1  ||
        rhs!=m) {
      fprintf(stderr, "read "ACGRED"%s!"ACRESET" failed\n", argv[2]);
      exit(2);
    }
    brr = new double[m];
  }


  for(int i=0; i<tnnz;i++){
    fscanf(f1, "%d %d %lf\n", airr+i, ajrr+i, avrr+i);
  }
  



  delete[] airr;
  delete[] ajrr;
  delete[] avrr;
  delete[] brr;

}
