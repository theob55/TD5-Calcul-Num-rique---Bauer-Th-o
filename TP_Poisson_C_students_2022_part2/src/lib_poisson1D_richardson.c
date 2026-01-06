/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
  // TODO: Compute all eigenvalues for the 1D Poisson operator
  for (int k=0; k<(*la); k++) {
    eigval[k] = 2-2*cos((k+1)*M_PI/((*la)+1));
  }
}

double eigmax_poisson1D(int *la){
  // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
  return 2-2*cos((*la)*M_PI/((*la)+1));
}

double eigmin_poisson1D(int *la){
  // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
  return 2-2*cos(M_PI/((*la)+1));
}

double richardson_alpha_opt(int *la){
  // TODO: Compute alpha_opt
  double* eigval = (double *) malloc(sizeof(double)*(*la));
  eig_poisson1D(eigval, la);
  double eigmax = eigmax_poisson1D(la);
  double eigmin = eigmin_poisson1D(la);
  free(eigval);
  return 2/(eigmin+eigmax);
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat
  if (cblas_dnrm2(*la, RHS, 1) == 0) {
    perror("||RHS|| = 0");
    return;
  }

  double* r=(double *) malloc(sizeof(double)*(*la));
  double my_norm = 100;
  double r_norm = 100;
  double nrm_RHS = cblas_dnrm2(*la, RHS, 1);
  while (my_norm > *tol && *nbite < *maxit) {
    cblas_daxpy(*la, *alpha_rich, r, 1, X, 1);
    cblas_dcopy(*la, RHS, 1, r, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1);
    r_norm = cblas_dnrm2(*la, r, 1);
    my_norm = r_norm / nrm_RHS;
    
    resvec[*nbite] = my_norm;
    (*nbite)++;
  }
  free(r);
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
  for (int i=0; i<*la; ++i) {
    MB[i*(*lab)] = 0;
    MB[1+i*(*lab)] = AB[1+i*(*lab)];
    MB[2+i*(*lab)] = 0;
  }
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
    for (int i=0; i<*la; ++i) {
      MB[i*(*lab)] = 0;
      MB[1+i*(*lab)] = AB[1+i*(*lab)];
      MB[2+i*(*lab)] = AB[2+i*(*lab)];
    }
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iterative method
  if (cblas_dnrm2(*la, RHS, 1) == 0) {
    perror("||RHS|| = 0");
    return;
  }

  double* r=(double *) malloc(sizeof(double)*(*la));
  double my_norm = 100;
  double r_norm = 100;
  double nrm_RHS = cblas_dnrm2(*la, RHS, 1);

  if (MB[2+(*lab)] == 0 && MB[2+(*la-1)*(*lab)] == 0) { //Jacobi.
    double* MB_inv = (double *) malloc(sizeof(double)*(*lab)*(*la));
    for (int i=0; i<(*la); i++) {
      MB_inv[i*(*lab)] = 0;
      MB_inv[1+i*(*lab)] = 1/(MB[1+i*(*lab)]); //main diag
      MB_inv[2+i*(*lab)] = 0;
    }
    while (my_norm > *tol && *nbite < *maxit) {
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, MB_inv, *lab, r, 1, 1.0, X, 1); //X += (MB)^-1 * r
      cblas_dcopy(*la, RHS, 1, r, 1); //r=RHS
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1); //r=RHS-A*X
      r_norm = cblas_dnrm2(*la, r, 1);
      my_norm = r_norm / nrm_RHS;
      
      resvec[*nbite] = my_norm;
      (*nbite)++;
    }
    free(MB_inv);
  }
  else { //Seidel
    int* ipiv = (int *) calloc(*la, sizeof(int));
    int NRHS = 1;
    int info=0;
    dgbtrftridiag(la, la, kl, ku, MB, lab, ipiv, &info); //factorize MB
    //MB(x^(k+1) - x^(k)) = (b - A*x^(k))
    while (my_norm > *tol && *nbite < *maxit) {
      cblas_dcopy(*la, RHS, 1, r, 1); //r=RHS
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1); //r=RHS-A*X
      dgbtrs_("N", la, kl, ku, &NRHS, MB, lab, ipiv, r, la, &info); //solves MB*y = r. stored in r
      cblas_daxpy(*la, 1.0, r, 1, X, 1); //x_k + x_{k+1} - x_k = x_{k+1}. stored in X
      cblas_dcopy(*la, RHS, 1, r, 1); //r=RHS
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1); //r=RHS-A*X

      r_norm = cblas_dnrm2(*la, r, 1);
      my_norm = r_norm / nrm_RHS;
      
      resvec[*nbite] = my_norm;
      (*nbite)++;
    }
    free(ipiv);
  }
  free(r);
}