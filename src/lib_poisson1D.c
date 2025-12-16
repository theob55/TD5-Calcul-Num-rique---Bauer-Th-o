/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int h){
//   // TODO: Fill AB with the tridiagonal Poisson operator
//   //LAPACKE_dgbsv(LAPACK_COL_MAJOR, n, &kl, &ku, &NRHS, AB, &lab, ipiv, b, ldb);

//   double values[] = {0, -1, 2, -1};

//   for (int j=0; j<(*la); ++j) {
//     for (int i=0; i<(*lab); ++i) {
//       AB[i + j*(*lab)] = values[i];
//     }
//   }

//   AB[1] = 0;
//   AB[(*la)*(*lab)-1] = 0;
// }

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator

  for(size_t i = 0; i <*la*(*lab); ++i) {
      AB[i] = 0.0; //initialisation de toute la matrice à 0;
  }

  for(size_t i = 1; i <*la; ++i) {
    AB[i * (*lab) +*kv] = -1.0;
  }

  for(size_t i = 0; i < *la; ++i) {
    AB[i*(*lab) + (*kv + 1)] = 2.0;
  }

  for(size_t i = 0; i < (*la - 1); ++i) {
    AB[i*(*lab) + (*lab - 1)] = -1.0;
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
  for (int i=0; i<*la; ++i) {
    for (int j=0; j<*lab; ++j) {
      AB[i*(*la)+j] = 0;
    }
    AB[i**kv] = 1;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
  double h2 = 1/((*la+1)*(*la+1));
  for (int i=0 ; i<*la; ++i) {
    RHS[i] = 0;
  }
  RHS[0] = *BC0*h2;
  RHS[*la-1] = *BC1*h2;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
  
  for (int i=0; i<*la; ++i) {
    EX_SOL[i] = -0.5*X[i]*X[i] + 0.5*X[i] + *BC0;
  }
  for (int i=0; i<*la; ++i) {
    EX_SOL[i] += *BC1-EX_SOL[*la-1];
  }
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
  x[0]=0;
  double la_inv = 1/(*la-1); 
  for (int i=1; i<*la; ++i) {
    x[i] = x[i-1] + la_inv;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  double* vec = malloc(*la*sizeof(double));
  for (int i = 0; i < *la; i++) {
     vec[i] = x[i];
  }
  cblas_daxpy(*la, -1.0, y, 1, vec, 1);
  double norm_up = cblas_dnrm2(*la, vec, 1);
  double norm_down = cblas_dnrm2(*la, x, 1);
  return norm_up/norm_down;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return (i-j+(*lab-1)) + j*(*lab);
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return *info;
}
