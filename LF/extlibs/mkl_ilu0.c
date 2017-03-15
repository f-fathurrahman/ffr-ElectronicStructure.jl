#include <stdio.h>
#include "mkl.h"

void mkl_ilu0( long long int n, long long int nnz, double *values, long long int *rowIndex,
               long long int *columns, double *bilu0 )
{
  long long int ipar[128];
  double dpar[128];
  long long int ierr;

  ipar[30] = 0;
  dpar[30] = 1.0e-10;
  //printf("n = %lli\n", n);
  //printf("nnz = %lli\n", nnz);

  //long long int i;
  //for( i = 0; i<nnz; i++){
  //  printf("%8lli %8lli %18.10f\n", i+1, columns[i], values[i]);
  //}

  dcsrilu0( &n , values, rowIndex, columns, bilu0, ipar, dpar, &ierr );

  //printf("mkl_ilu0: ierr = %lli\n", ierr);
}


void apply_prec_ILU0( long long int N, double *bilu0, long long int *ia, long long int *ja,
                      double *v, double *pv )
{
  char cvar1, cvar, cvar2;
  double *tmp = (double*)malloc( N*sizeof(double) );

  //printf("N = %lld\n", N);

  cvar1 = 'L';
  cvar  = 'N';
  cvar2 = 'U';
  mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &N, bilu0, ia, ja, v, tmp);

  cvar1 = 'U';
  cvar = 'N';
  cvar2 = 'N';
  mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &N, bilu0, ia, ja, tmp, pv);

  free(tmp);

  return;
}

