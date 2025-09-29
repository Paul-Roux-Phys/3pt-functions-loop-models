#include "common.h"
#include "hash.h"
#include "transfer.h"

void main()
{
  int i,j;
  int k1,k2,k3;
  double Z123,Z101,Z303,C123,C220;
  double Z123_exponent,Z101_exponent,Z303_exponent;
  double lambda_over_pi,lambda;

  init_hash();
  // data_position = 0;

#ifdef RESUME
  lambda_over_pi = RESUME;
#else
  lambda_over_pi = 0.13;
#endif

  for(; lambda_over_pi < 0.4951; lambda_over_pi += 0.005)
    {
      /* Weights of Warnaar-Nienhuis-Seaton */
      lambda = PI*lambda_over_pi;
      n_loop = -2*cos(4*lambda);
      rho[0] = 0.0; // Not used
      rho[1] = 1 + sin(lambda) + sin(3*lambda) - sin(5*lambda);
      rho[2] = 2*sin(2*lambda)*sin((6*lambda+PI)/4);
      rho[3] = rho[2];
      rho[4] = rho[2];
      rho[5] = rho[2];
      rho[6] = 1 + sin(3*lambda);
      rho[7] = rho[6];
      rho[8] = sin(lambda) + cos(2*lambda);
      rho[9] = rho[8];

      k1 = kk1;
      k2 = 0;
      k3 = kk3;
      
      /* CASE 1: Z(k1,k2,k3) */
      Z123 = compute_Z_special(k1,k2,k3);
      Z123_exponent = exponent;
      // printf("Z123 %.16f %.16f\n",Z123,Z123_exponent);

      /* CASE 2: Z(k1,0,k1) */
      Z101 = compute_Z(k1,0,k1);
      Z101_exponent = exponent;
      // printf("Z101 %.16f %.16f\n",Z101,Z101_exponent);
      
      /* CASE 3: Z(k3,0,k3) */
      Z303 = compute_Z(k3,0,k3);
      Z303_exponent = exponent;
      // printf("Z303 %.16f %.16f\n",Z303,Z303_exponent);
      
      exponent = Z123_exponent - 0.5*(Z101_exponent + Z303_exponent);
      C123 = Z123 / sqrt(Z101 * Z303) * exp(exponent);
      
      k1 = 0;
      k2 = 0;
      k3 = 0;
      
      /* CASE 1: Z(k1,k2,k3) */
      Z123 = compute_Z_seam(1,0);
      Z123_exponent = exponent;
      // printf("Z220 %.16f %.16f\n",Z123,Z123_exponent);
      
      /* CASE 2: Z(k1,0,k1) */
      Z101 = compute_Z_seam(1,1);
      Z101_exponent = exponent;
      // printf("Z202 %.16f %.16f\n",Z101,Z101_exponent);
      
      /* CASE 3: Z(k3,0,k3) */
      Z303 = compute_Z_seam(0,0);
      Z303_exponent = exponent;
      // printf("Z000 %.16f %.16f\n",Z303,Z303_exponent);
      
      exponent = Z123_exponent - 0.5*(Z101_exponent + Z303_exponent);
      C220 = Z123 / sqrt(Z101 * Z303) * exp(exponent);
      
      printf("%2d %.6f %.10f %2d %2d %.10f  %.16f \n",n,lambda_over_pi,n_loop,kk1,kk3,n2,C123/C220);
      fflush(NULL);
    }
}
