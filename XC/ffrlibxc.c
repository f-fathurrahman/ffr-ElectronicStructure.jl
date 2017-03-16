#include <stdio.h>

#include <xc.h>

void ffrlibxc( int func_id, int N, double *rho, double *sigma, double *exc )
{
  xc_func_type func;
  int i, vmajor, vminor, vmicro;

  //xc_version(&vmajor, &vminor, &vmicro);
  //printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);

  xc_func_init(&func, func_id, XC_UNPOLARIZED);

  //printf("The functional is '%s'\n", func.info->name);

  switch(func.info->family)
  {
    case XC_FAMILY_LDA:
      xc_lda_exc(&func, N, rho, exc);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc(&func, N, rho, sigma, exc);
      break;
  }

  //for(i=0; i<N; i+=1){
  //  printf("C: %lf %lf\n", rho[i], exc[i]);
  //}

  xc_func_end(&func);
}
