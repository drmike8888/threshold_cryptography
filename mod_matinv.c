/*************************************************************************************
 *                                                                                   *
 *      Simple Gauss-Jordan matrix inversion.                                        *
 *  Enter with pointer to multiprecion data in row*n + column format, dimension n    *
 *  and place to store result.                                                       *
 *  Returns +1 if ok, -1 if singular.                                                *
 *                                                                                   *
 *                        Author = Mike Rosing                                       *
 *                        Date = Feb. 16, 2026                                       *
 *                                                                                   *
 ************************************************************************************/

#include "modulo.h"

/* inv should already be initialized to n x n mpz_t values.
   This assumes matrix has different prime than elliptic curves,
   so prime value is also an input.
*/

int mod_matinv(mpz_t *inv, long n, mpz_t *mat, mpz_t prime)
{
  mpz_t *work, rowmod, prevprm, tmp, pvt;
  long i, j, k;

  work = (mpz_t*)malloc(sizeof(mpz_t)*2*n*n);

/* set up modular routines for input prime modulus */
  
  mget(prevprm);
  mset(prime);
  
/*  intitialize matrix  */

  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
    {
      mpz_init(work[2*i*n + j]);
      mpz_set(work[i*2*n + j], mat[i*n + j]);
      mpz_init(work[2*i*n + n + j]);
      mpz_set_ui(work[i*2*n + n + j], 0);
      if(i == j) mpz_set_ui(work[i*2*n + n + j], 1);
    }
  }

/* pivot is diagonal on each row */
  
  mpz_inits(pvt, rowmod, tmp, NULL);
  for(k=0; k<n; k++)
  {
  
/*  normalize row with diagonal  */

    mpz_set(pvt, work[2*n*k + k]);
    if(!mpz_cmp_ui(pvt, 0))
    {
      for(i=0; i<n; i++)
	for(j=0; j<2*n; j++)
	  mpz_clear(work[2*i*n + j]);
      free(work);
      mset(prevprm);
      mpz_clears(pvt, rowmod, prevprm, tmp, NULL);
      return -1;
    }
    for(j=0; j<2*n; j++)
      mdiv(work[2*n*k + j], work[2*n*k + j], pvt);

  /*   printf("normalize k: %ld\n", k); */
  /* for(i=0; i<n; i++) */
  /* { */
  /*   for(j=0; j<2*n; j++) */
  /*     gmp_printf("%Zd ", work[i*2*n + j]); */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */
/*  eliminate all rows  */

    for(i=0; i<n; i++)
    {
      if(i == k) continue;
      mpz_set(rowmod, work[i*2*n + k]);
      for(j=0; j<2*n; j++)
      {
	mmul(tmp, rowmod, work[k*2*n + j]);
	msub(work[i*2*n + j], work[i*2*n + j], tmp);
      }
    }
  /*   printf("eliminate k: %ld\n", k); */
  /* for(i=0; i<n; i++) */
  /* { */
  /*   for(j=0; j<2*n; j++) */
  /*     gmp_printf("%Zd ", work[i*2*n + j]); */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */
  }   // end k loop

/*  use row and col vectors to unscramble result.
    row is from, col is to.
*/

  for(k=0; k<n; k++)
  {
    for(j=0; j<n; j++)
      mpz_set(inv[k*n + j], work[k*2*n + n + j]);
  }
  for(i=0; i<n; i++)
    for(j=0; j<2*n; j++)
      mpz_clear(work[2*i*n + j]);
  free(work);
  mset(prevprm);
  mpz_clears(pvt, rowmod, prevprm, tmp, NULL);
  return 1;
}

/* compute matrix multiply c(m x s) = a(m x n) * b(r x s)
   if n != r die horribly.
   assumes c is already initialized.
*/

int mod_matmul(mpz_t *c, mpz_t *a, long m, long n,
	       mpz_t *b, long r, long s, mpz_t prime)
{
  mpz_t prevprm, tmp;
  long i, j, k;

  if(n != r)
  {
    printf("incompatible matricies for modular multiplication\n");
    exit(-17);
  }

/* set up modular routines for input prime modulus */
  
  mget(prevprm);
  mset(prime);

  mpz_init(tmp);
  for(i=0; i<m; i++)
  {
    for(j=0; j<s; j++)
    {
      mpz_set_ui(c[i*s + j], 0);
      for(k=0; k<n; k++)
      {
	mmul(tmp, a[i*n + k], b[k*s + j]);
	madd(c[i*s + j], c[i*s + j], tmp);
      }
    }
  }

  mset(prevprm);
  mpz_clears(tmp, prevprm, NULL);
}
