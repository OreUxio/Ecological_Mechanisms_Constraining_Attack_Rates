//$Id: linpack_eigen.cc 2057 2011-01-22 12:32:41Z axel $
// the gsl_eigen_symmv seems to be broken!! dont use it 
#include <string.h> //for memcpy
#include <stdio.h>
#include <stdlib.h>
#include "linpack_eigen.h"
#include "error.h"

typedef int integer;
typedef double doublereal;
typedef char character;
typedef float real;

extern "C" {
  void dsyev_(character *,character *,integer *, doublereal *, 
	      integer *,doublereal *, doublereal *,integer *,integer *);
  void dsyevr_(character* JOBZ,character* RANGE,character* UPLO,  
	       integer* N, doublereal* A,integer*  LDA,doublereal* VL,
	       doublereal* VU,integer*  IL,integer*  IU,doublereal* ABSTOL,
	       integer* M,doublereal*  W,doublereal*  Z,integer* LDZ,
	       integer* ISUPPZ,doublereal* WORK, 
	       integer* LWORK,integer* IWORK,integer* LIWORK,
	       integer* INFO );
  void dgeev_(character * JOBVL,character * JOBVR,integer *N, doublereal *A, 
	      integer *LDA,doublereal *WR,doublereal *WI,doublereal *VL,
	      integer *LDVL,doublereal *VR,integer *LDVR,doublereal *WORK,
	      integer *LWORK, integer *INFO);
  void dpotrf_(character* UPLO,integer* N,doublereal* A,
	       integer* LDA,integer* INFO );
  void spotrf_(character* UPLO,integer* N,real* A,
	       integer* LDA,integer* INFO );
}

// get eigenvalues of real mat into ddr and ddi:
void eigen(const CLHEP::HepMatrix & mat, CLHEP::HepVector & ddr, CLHEP::HepVector & ddi){
  
//       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,  LDVR,
//                          WORK, LWORK, INFO )
//            CHARACTER     JOBVL, JOBVR
//            INTEGER       INFO, LDA, LDVL, LDVR, LWORK, N
//            DOUBLE        PRECISION  A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
//                          WI( * ), WORK( * ), WR( * )


  int n=mat.num_row();
  if(n!=mat.num_col()) FATAL_ERROR("");

  ddr=CLHEP::HepVector(n);
  ddi=CLHEP::HepVector(n);
  double *dr=&(ddr[0]);
  double *di=&(ddi[0]);

  character jobvl='N',jobvr='N'; // don't compute eigenvectors

  integer info,lda=n,lwork;
  static int largest_n=1;
  static doublereal * a = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  static doublereal * work = (doublereal *) malloc(sizeof(doublereal));
  static int workspacesize=0;

  //query optimal workspacesize
  lwork=-1;
  dgeev_(&jobvl,&jobvr,&n, a,&lda,dr,di,0,&lda,0,&lda,work,&lwork,&info);

  if(info!=0){
    REPORT(info);
    WARNING("error calling lapack routine dgeev");
  }
  
  if(*work>workspacesize){
    workspacesize=int(*work+0.5);
    work=(doublereal *)realloc(work,sizeof(doublereal)*workspacesize);
    if(work==0){
      FATAL_ERROR("memory allocation problem");
    }
  }

  if(n>largest_n){
    largest_n=n;
    a=(doublereal *)realloc(a,sizeof(doublereal)*largest_n*largest_n);
    if(a==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  for(int i=0;i<n;i++){
    memcpy(a+i*lda,&(mat[i][0]),n*sizeof(double));
  }

  lwork=workspacesize;
  dgeev_(&jobvl,&jobvr,&n, a,&lda,dr,di,0,&lda,0,&lda,work,&lwork,&info);
  if(info!=0){
    REPORT(info);
    WARNING("error calling lapack routine dgeev");
  }

  return;
}

// get eigenvalues of real mat into ddr and ddi and its left and right
// eigenvectors (following LAPACK convention to split complex matrices
// into successive real and imaginary parts), into vecs_l and vecs_r,
// respectively.
void eigen(const CLHEP::HepMatrix & mat, CLHEP::HepVector & ddr, CLHEP::HepVector & ddi, CLHEP::HepMatrix & vecs_l, CLHEP::HepMatrix & vecs_r){
  
//       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,  LDVR,
//                          WORK, LWORK, INFO )
//            CHARACTER     JOBVL, JOBVR
//            INTEGER       INFO, LDA, LDVL, LDVR, LWORK, N
//            DOUBLE        PRECISION  A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
//                          WI( * ), WORK( * ), WR( * )


  int n=mat.num_row();
  if(n!=mat.num_col()) FATAL_ERROR("");

  ddr=CLHEP::HepVector(n);
  ddi=CLHEP::HepVector(n);
  double *dr=&(ddr[0]);
  double *di=&(ddi[0]);

  character jobvl='V',jobvr='V'; // don't compute eigenvectors

  integer info,lda=n,lwork;
  static int largest_n=1;
  static doublereal * a = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  static doublereal * vl = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  static doublereal * vr = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  static doublereal * work = (doublereal *) malloc(sizeof(doublereal));
  static int workspacesize=0;

  //query optimal workspacesize
  lwork=-1;
  dgeev_(&jobvl,&jobvr,&n, a,&lda,dr,di,0,&lda,0,&lda,work,&lwork,&info);

  if(info!=0){
    REPORT(info);
    WARNING("error calling lapack routine dgeev");
  }
  
  if(*work>workspacesize){
    workspacesize=int(*work+0.5);
    work=(doublereal *)realloc(work,sizeof(doublereal)*workspacesize);
    if(work==0){
      FATAL_ERROR("memory allocation problem");
    }
  }

  if(n>largest_n){
    largest_n=n;
    a=(doublereal *)realloc(a,sizeof(doublereal)*largest_n*largest_n);
    if(a==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
    vr=(doublereal *)realloc(vr,sizeof(doublereal)*largest_n*largest_n);
    if(vr==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
    vl=(doublereal *)realloc(vl,sizeof(doublereal)*largest_n*largest_n);
    if(vl==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  for(int i=0;i<n;i++){
    memcpy(a+i*lda,&(mat[i][0]),n*sizeof(double));
  }

  lwork=workspacesize;
  dgeev_(&jobvl,&jobvr,&n, a,&lda,dr,di,vl,&lda,vr,&lda,work,&lwork,&info);
  if(info!=0){
    REPORT(info);
    WARNING("error calling lapack routine dgeev");
  }

  // we do not transpose here (!) just exchange r and l
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      vecs_r[i][j]=vl[i*lda+j];
      vecs_l[i][j]=vr[i*lda+j];
    }
  }
  return;
}

// get eigenvalues of real symmetric mat into dd
void eigen(CLHEP::HepMatrix & mat, CLHEP::HepVector & dd){
  
//   SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
//            CHARACTER     JOBZ, UPLO
//            INTEGER       INFO, LDA, LWORK, N
//            DOUBLE        PRECISION A( LDA, * ), W( * ), WORK( * )

  int n=mat.num_row();
  if(n!=mat.num_col()) FATAL_ERROR("non-square matrix");

  dd=CLHEP::HepVector(n);
  double *d=&(dd[0]);

  character jobz='V',uplo='U'; // the matrix mat is lower, but in
			       // FORTRAN rows and columns are
			       // exchanced
  integer info,lda=n,lwork;
  static int largest_n=1;
  static doublereal * a = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  static doublereal * work = (doublereal *) malloc(sizeof(doublereal));
  static int workspacesize=0;

  //query optimal workspacesize
  lwork=-1;
  dsyev_(&jobz,&uplo,&n,a,&lda,d,work,&lwork,&info);
  if(info<0){
    REPORT(info);
    WARNING("error calling lapack routine dsyev");
  }else if(info>0){
    WARNING("eigensystem computation did not converge");
  }
  
  if(*work>workspacesize){
    workspacesize=int(*work+0.5);
    work=(doublereal *)realloc(work,sizeof(doublereal)*workspacesize);
    if(work==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  if(n>largest_n){
    largest_n=n;
    a=(doublereal *)realloc(a,sizeof(doublereal)*largest_n*largest_n);
    if(a==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  int i,j;
  
  for(i=0;i<n;i++){
    memcpy(a+i*lda,&(mat[i][0]),n*sizeof(double));
  }

  lwork=workspacesize;
  dsyev_(&jobz,&uplo,&n,a,&lda,d,work,&lwork,&info);
  if(info<0){
    REPORT(info);
    WARNING("error calling lapack routine dsyev");
  } if(info>0){
    WARNING("eigensystem computation did not converge");
  }
  


  // we have to transpose here:
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      mat[j][i]=a[i*lda+j];
    }
  }
  return;
}


// Get the m largest eigenvalues of real symmetric mat into dd???
void eigen(CLHEP::HepMatrix & mat, CLHEP::HepVector & dd,double accuracy,int m){
  
//    SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO,  N,  A,  LDA,  VL,  VU,  IL,  IU,
//                       ABSTOL,  M,  W,  Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
//                       LIWORK, INFO )
//        CHARACTER      JOBZ, RANGE, UPLO
//        INTEGER        IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
//        DOUBLE         PRECISION ABSTOL, VL, VU
//        INTEGER        ISUPPZ( * ), IWORK( * )
//        DOUBLE         PRECISION A( LDA, * ), W( * ), WORK( * ), Z( LDZ,  *)


  int n=mat.num_row();
  if(n!=mat.num_col()) FATAL_ERROR("non-square matrix");

  dd=CLHEP::HepVector(n);
  double *d=&(dd[0]);

  character jobz='V',uplo='U'; // the matrix mat is lower, but in
			       // FORTRAN rows and columns are
			       // exchanced
  character range=(m?'I':'A');
  integer info,lda=n,lwork,liwork,n_eigenvalues;
  integer ilower=n-m+1,iupper=n;
  doublereal rdummy,abstol=accuracy;
  static int largest_n=1;
  static doublereal * a = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  static doublereal * z = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  static doublereal * work = (doublereal *) malloc(sizeof(doublereal));
  static integer * iwork = (integer *) malloc(sizeof(integer));
  static int workspacesize=0;
  static int iworkspacesize=0;
  integer isuppz[2*n];//range of non-zero entries in eigenvectors, NOT USED.

  //query optimal workspacesize
  lwork=-1;
  liwork=-1;
  dsyevr_(&jobz,&range,&uplo,&n,a,&lda,&rdummy,&rdummy,&ilower,&iupper,&abstol,
	  &n_eigenvalues,d,z,&lda,isuppz,work,&lwork,iwork,&liwork,&info);
  if(info<0){
    REPORT(info);
    WARNING("error calling lapack routine dsyevr");
  }else if(info>0){
    WARNING("eigensystem computation did not converge");
  }
  
  if(*work>workspacesize){
    workspacesize=int(*work+0.5);
    work=(doublereal *)realloc(work,sizeof(doublereal)*workspacesize);
    if(work==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  if(*iwork>iworkspacesize){
    iworkspacesize=*iwork;
    iwork=(integer *)realloc(iwork,sizeof(integer)*iworkspacesize);
    if(iwork==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  if(n>largest_n){
    largest_n=n;
    a=(doublereal *)realloc(a,sizeof(doublereal)*largest_n*largest_n);
    if(a==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
    z=(doublereal *)realloc(z,sizeof(doublereal)*largest_n*largest_n);
    if(z==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  int i,j;
  
  for(i=0;i<n;i++){
    memcpy(a+i*lda,&(mat[i][0]),n*sizeof(double));
  }

  lwork=workspacesize;
  liwork=iworkspacesize;
  dsyevr_(&jobz,&range,&uplo,&n,a,&lda,&rdummy,&rdummy,&ilower,&iupper,&abstol,
	  &n_eigenvalues,d,z,&lda,isuppz,work,&lwork,iwork,&liwork,&info);
  if(m) ALWAYS_ASSERT(n_eigenvalues==m);
  m=n_eigenvalues;
  if(info<0){
    REPORT(info);
    WARNING("error calling lapack routine dsyevr");
    exit(1);
  } if(info>0){
    WARNING("lapack internal error");
  }
  
  // clear superfluous entries:
  for(int i=m;i<n;i++){
    d[i]=0;
  }

  // we have to transpose here:
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      mat[j][i]=(i<m?z[i*lda+j]:0);
    }
  }
  return;
}

int CholeskyL(CLHEP::HepMatrix & mat){

//   SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
//     CHARACTER      UPLO
//     INTEGER        INFO, LDA, N
//     DOUBLE         PRECISION A( LDA, * )

  int n=mat.num_row();
  if(n!=mat.num_col()) FATAL_ERROR("non-square matrix");

  character uplo='L';

  integer info,lda=n;
  static int largest_n=1;
  static doublereal * a = (doublereal *) malloc(sizeof(doublereal)*largest_n*largest_n);
  if(n>largest_n){
    largest_n=n;
    a=(doublereal *)realloc(a,sizeof(doublereal)*largest_n*largest_n);
    if(a==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  int i,j;
  
  for(i=0;i<n;i++){
    memcpy(a+i*lda,&(mat[i][0]),n*sizeof(double));
  }

  dpotrf_(&uplo,&n,a,&lda,&info);
  if(info<0){
    REPORT(info);
    WARNING("error calling lapack routine dpotrf");
  } if(info>0){
    WARNING("could not do this decomposition");
  }

  // we have to transpose here:
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      mat[j][i]=(i>j?0:a[i*lda+j]);
    }
  }
  return info;
}

int rCholeskyL(CLHEP::HepMatrix & mat){

//   SUBROUTINE SPOTRF( UPLO, N, A, LDA, INFO )
//     CHARACTER      UPLO
//     INTEGER        INFO, LDA, N
//     REAL           A( LDA, * )

  int n=mat.num_row();
  if(n!=mat.num_col()) FATAL_ERROR("non-square matrix");

  character uplo='L';

  integer info,lda=n;
  static int largest_n=1;
  static real * a = (real *) malloc(sizeof(real)*largest_n*largest_n);
  if(n>largest_n){
    largest_n=n;
    a=(real *)realloc(a,sizeof(real)*largest_n*largest_n);
    if(a==0){
      fprintf(stderr,"memory allocation problem\n");
      exit(1);
    }
  }

  int i,j;
  
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      a[i*lda+j]=mat[j][i];
    }
  }

  spotrf_(&uplo,&n,a,&lda,&info);
  if(info<0){
    REPORT(info);
    WARNING("error calling lapack routine spotrf");
  } if(info>0){
    WARNING("could not do this decomposition");
  }

  // we have to transpose here:
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      mat[j][i]=(i>j?0:a[i*lda+j]);
    }
  }
  return info;
}


// Multipliation of HepMatrices using BLAS:
extern "C" {
  void dgemm_(const char*,const char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
}
// recall that, as compared to FORTAN, C matrices are transposed
CLHEP::HepMatrix operator^(CLHEP::HepMatrix B, CLHEP::HepMatrix A){
  int m=A.num_col() ,n=B.num_row(),k=A.num_row();
  CLHEP::HepMatrix C(n,m);
  static double alpha=1,beta=0;
  if(B.num_col()!=k){
    FATAL_ERROR("matrix dimensions don't fit");
  }
  dgemm_("n","n",&m,&n,&k,&alpha,&A[0][0],&m,&B[0][0],&k,&beta,&C[0][0],&m);
  return C;
}
