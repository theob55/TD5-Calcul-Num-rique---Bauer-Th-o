/******************************************/
/* tp_env.c                               */
/* This file contains a main function to  */
/* test the environment of compilation    */
/* and verify BLAS/LAPACK installation    */
/******************************************/
#include "tp_env.h"

/**
 * Main function to test the compilation environment and computed functions.
 * Prints various system constants and tests BLAS functions.
 * Prints result of all tests for each implemented function
 * 
 * @param argc: Number of command-line arguments (unused)
 * @param argv: Array of argument strings (unused)
 * @return 0 on success
 */
int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  printf("--------- Test environment of execution for Practical exercises of Numerical Algorithmics ---------\n\n");
  
  /* Display mathematical and precision constants */
  printf("The exponantial value is e = %f \n",M_E);
  printf("The maximum single precision value from values.h is maxfloat = %e \n",MAXFLOAT);
  printf("The maximum single precision value from float.h is flt_max = %e \n",FLT_MAX);
  printf("The maximum double precision value from float.h is dbl_max = %e \n",DBL_MAX);
  printf("The epsilon in single precision value from float.h is flt_epsilon = %e \n",FLT_EPSILON);
  printf("The epsilon in double precision value from float.h is dbl_epsilon = %e \n",DBL_EPSILON);

  printf("\n\n Test of ATLAS (BLAS/LAPACK) environment \n");

  /* Test BLAS functionality */
  double x[5], y[5];  /* Test vectors */
  int ii;
  
  /* Initialize test vectors */
  for (ii=0;ii<5;ii++){
    x[ii]=ii+1; y[ii]=ii+6;
    printf("x[%d] = %lf, y[%d] = %lf\n",ii,x[ii],ii,y[ii]);
  }

  /* Test cblas_dcopy: copies vector x into vector y */
  printf("\nTest DCOPY y <- x \n");
  cblas_dcopy(5,x,1,y,1);  /* Copy 5 elements from x to y with stride 1 */
  for (ii=0;ii<5;ii++){
    printf("y[%d] = %lf\n",ii,y[ii]);
  }

  printf("\n\nTest of functions in lib_poisson1D.c \n");

  int NRHS=1;           /* Solving Ax=b with one right-hand side */
  int nbpoints=10;      /* Total number of discretization points (including boundaries) */
  int la=nbpoints-2;    /* Number of interior points (excluding boundaries) */
  double T0=-5.0;          /* Dirichlet boundary condition at x=0 */
  double T1=5.0;           /* Dirichlet boundary condition at x=1 */
  double h = 1.0/(la+1);

  /* Allocate memory for vectors */
  double *RHS=(double *) malloc(sizeof(double)*la);      /* Right-hand side vector */
  double *EX_SOL=(double *) malloc(sizeof(double)*la);   /* Analyticalexact solution */
  double *X=(double *) malloc(sizeof(double)*la);        /* Grid points */
  double *APPROX_SOL=(double *) malloc(sizeof(double)*la);   /* Approximated solution */
  double *res_ex=(double *) malloc(sizeof(double)*la);   /* residu calculé avec EX_SOL */

  for (int i=0; i<la; ++i) {
    res_ex[i] = 0;
    X[i] = (i+1)*h;
    EX_SOL[i] = T0 + X[i]*(T1 - T0); //Calcul théorique de la solution exacte
  }
  RHS[0] += T0/(h*h);
  RHS[la-1] += T1/(h*h);
  
  /* Set up band storage parameters for tridiagonal matrix */
  int kv=1;             /* Number of superdiagonals */
  int ku=1;             /* Number of superdiagonals in original matrix */
  int kl=1;             /* Number of subdiagonals */
  int lab=2*kl+ku+1;   /* Leading dimension of band storage */

  /* Allocate and initialize the coefficient matrix */
  double *AB1 = (double *) malloc(sizeof(double)*lab*la);

  int *ipiv = (int *) calloc(la, sizeof(int));  /* Pivot indices for LU factorization */

  set_GB_operator_colMajor_poisson1D(AB1, &lab, &la, &kv); //Sets AB1 in ColMajor
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB1, lab, EX_SOL, 1, -1.0, res_ex, 1); //res_ex = AB1*EX_SOL
  
  double res_ex_nrm = cblas_dnrm2(la, res_ex, 1); // res_ex_nrm = ||res_ex||

  printf("Exercice 4 : ");
  if (res_ex_nrm < 1e-12) {
    printf("SUCCESS\n"); //Methode de validation : exercice 4
  }
  else {printf("FAIL\n");}

  //Exercice 5 :
  //Nous avons deja set AB1 ColMajor
  int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB1, lab, ipiv, RHS, la); //info:code d'erreur ; RHS = solution approchée

  if (info != 0) {
    printf("ERROR : info of dgbsv != 0 \n");
    return 1;
  }
  printf("Exercice 5 :\n");
  for (int i=0; i<la; ++i) {
    printf("sol_approx[%d] = %f\n", i, RHS[i]);
  }

  printf("----------- freeing ---------\n");
  free(RHS); free(EX_SOL); free(X); free(APPROX_SOL); free(res_ex); free(AB1); free(ipiv);
  printf("\n\n--------- End -----------\n");
}

/*
runs with :
docker build -f docker/Dockerfile --progress=plain -t tp-cn:latest .
docker run -it tp-cn:latest
make run
*/
