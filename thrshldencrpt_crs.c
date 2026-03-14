/************************************************************
 *                                                          *
 *   Subroutines specific to threshold encryption algorithm *
 *   which assumes polynomial basis elliptic curves exist.  *
 *                                                          *
 *                       Author = Mike Rosing               *
 *                        Date = 17 Feb. 2026               *
 *                                                          *
 ***********************************************************/

#include "poly_eliptic.h"
#include "eliptic.h"
#include "thrshldencrpt_setup.h"

/* create share generating matrix.
   Enter with 2N-1 x N matrix initialized to hold values
   and value N (which is assumed to be 2^n) as well as
   modulus (assumed to be torsion of elliptic curve).
*/

void genshare(mpz_t *share, long N, mpz_t mod)
{
  long i, j, n1;
  mpz_t x;

  n1 = 2*N - 1;
  mpz_init(x);
  for(i=0; i<n1; i++)
  {
    mpz_set_ui(x, i + 1);
    for(j=0; j<N; j++)
      mpz_powm_ui(share[i*N + j], x, j, mod);
  }
  mpz_clear(x);
}

/* compute inner product lth row with secret vector.
   input share matrix, row l, s_vec and N. result ms
   should already be initialized.
*/

void mdots(mpz_t ms, long row, mpz_t *s_vec, mpz_t *shrmtrx, long N, mpz_t mod)
{
  mpz_t prevprm, tmp;
  long i, rn;
  
/* set up modular routines for input prime modulus */
  
  mget(prevprm);
  mset(mod);
  mpz_init(tmp);
  rn = row*N;
  mpz_set_ui(ms, 0);
  for(i=0; i<N; i++)
  {
    mmul(tmp, s_vec[i], shrmtrx[rn + i]);
    madd(ms, ms, tmp);
  }

  mset(prevprm);
  mpz_clears(prevprm, tmp, NULL);
}

/*  compute zhat_0 and v_l,0 vector from share matrix and ghat.
    assumes elliptic curve already set up in CRS. Input CRS,
    share matrix, secrets c and s_vec,  and number of users N.
    In addition to computing values inside CRS, returns c_vec
    for use in computing y_hat and tau_hat values. c_vec is
    large: from -2N + 2 to 5N - 1 is 7N - 3 entries.
*/

void zvcalc(mpz_t *cvec, CRS *crs, mpz_t *shrmtrx, mpz_t c, mpz_t *s_vec, long N)
{
  long i, ldx, n1, lofst, n2;
  mpz_t *mdsvec, tmp, coef;

/* start with m . s for each row of share matrix */
  
  n1 = 2*N - 1;
  mdsvec = (mpz_t*)malloc(sizeof(mpz_t)*n1);
  for(i=0; i<n1; i++)
  {
    mpz_init(mdsvec[i]);
    mdots(mdsvec[i], i, s_vec, shrmtrx, N, crs->grp->tor);
  }

/* create c^j for all j in [-2N+1 ... 4N-2]  */

  mpz_init(cvec[0]);
  mpz_powm_ui(cvec[0], c, n1, crs->grp->tor);   // negative n1 failed??
  mpz_invert(cvec[0], cvec[0], crs->grp->tor);
  for(i=1; i<6*N; i++)
  {
    mpz_init(cvec[i]);
    mod_mul(cvec[i], c, cvec[i - 1], crs->grp->tor);
  }

/* compute z0hat using cvec from 1 to 2N-1 and all of m.s */

  printf("computing z0 hat\n");
  n2 = 2*N;
  mpz_inits(tmp, coef, NULL);
  for(i=0; i<n1; i++)
  {
    mod_mul(tmp, cvec[i + n2], mdsvec[i], crs->grp->tor);
    mod_add(coef, coef, tmp, crs->grp->tor);
  }
  poly_point_init(&crs->z0hat);
  poly_elptic_mul(&crs->z0hat, crs->ghat, coef, crs->Ex);

/* compute v_l,0 for each l in [1 ... 2N-1] */

  for(ldx=0; ldx<n1; ldx++)
  {
    printf("computing v_%ld,0 hat\r", ldx);
    fflush(stdout);
    poly_point_init(&crs->vl0[ldx]);
    mpz_set_ui(coef, 0);
    for(i=0; i<n1; i++)
    {
      if(i == ldx) continue;
      lofst = i - ldx + n1;  // 2N-1 index is 0 exponent
      mod_mul(tmp, cvec[lofst], mdsvec[i], crs->grp->tor);
      mod_add(coef, coef, tmp, crs->grp->tor);
    }
    poly_elptic_mul(&crs->vl0[ldx], crs->ghat, coef, crs->Ex);
  }
  printf("\n");
  
/* clear out temporary values */

  mpz_clears(tmp, coef, NULL);
  for(i=0; i<n1; i++)
    mpz_clear(mdsvec[i]);
  free(mdsvec);
}

/* compute y_hat vector and tau_hat matrix. Use cvec from 
   zvcalc. Enter with CRS, cvec, n and N. 
   NOTE: gamma values are toxic waste - use and destroy 
   securely.  */

void ytaucalc(CRS *crs, mpz_t *cvec, long n, long N)
{
  FILE *dbg;
  long *X, i, j, k, n01, n3, n4;
  mpz_t tmp, coef, *gamma;

  n01 = N - 1;
  
/* create index list for each bit */
  
  X = (long*)malloc(sizeof(long)*(n + 1));
  X[0] = 1;
  for(i=1; i<=n; i++)
    X[i] = 2*X[i - 1];
  
/* create set of random values */

  gamma = (mpz_t*)malloc(sizeof(mpz_t)*n01);
  dbg = fopen("debug.gamma", "w");
  for(i=0; i<n01; i++)
  {
    mpz_init(gamma[i]);
    mod_rand(gamma[i], crs->grp->tor);
    mpz_out_raw(dbg, gamma[i]);
  }
  fclose(dbg);
/* compute y_hat for each j < n */

  n3 = 3*N - 1;              // gets to exponent X_m + N + index
  mpz_inits(tmp, coef, NULL);
  for(j=0; j<n; j++)
  {
    printf("computing y_hat[%ld]\r", j);
    fflush(stdout);
    poly_point_init(&crs->yhat[j]);
    i = X[j];
    mpz_set_ui(coef, 0);
    while(i < X[j + 1])
    {
      mod_mul(tmp, gamma[i - 1], cvec[i + n3], crs->grp->tor);
      mod_add(coef, coef, tmp, crs->grp->tor);
      i++;
    }
    poly_elptic_mul(&crs->yhat[j], crs->ghat, coef, crs->Ex);
  }
  printf("\n");
  
/* compute tau_hat matrix for each j < n and i in [-2N+1 .. 2N-1] */

  n3 = 5*N - 1;
  n4 = 4*N - 1;
  for(j=0; j<n; j++)
  {
    for(i=0; i<n4; i++)
    {
      printf("computing tau[%ld, %ld]\r", j, i);
      fflush(stdout);
      if((i < 3*N - 1 + X[j]) || (i >= 3*N - 1 + X[j + 1]))
      {
	poly_point_init(&crs->tauhat[j*n4 + i]);
	k = X[j] - 1;
	mpz_set_ui(coef, 0);
	while(k < X[j + 1] - 1)
        {
	  mod_mul(tmp, gamma[k], cvec[k - i + n3], crs->grp->tor);
	  mod_add(coef, coef, tmp, crs->grp->tor);
	  k++;
	}
	poly_elptic_mul(&crs->tauhat[j*n4 + i], crs->ghat, coef, crs->Ex);
      }
      else // exclude X_j values in this row
	poly_point_init(&crs->tauhat[j*n4 + i]); // point at infinity
    }
  }
  printf("\n");
  mpz_clears(tmp, coef, NULL);
  for(i=0; i<n01; i++)
    mpz_clear(gamma[i]);
  free(gamma);
}

/* compute g_hat*c_vec[i] and g*c_vec[i] for i in -2N+1 to 2N-1 */

void cgcalc(CRS *crs, mpz_t *cvec, long N)
{
  long i, n1;

  n1 = 4*N - 1;
  for(i=0; i<n1; i++)
  {
    printf("computing cg[%ld]\r", i);
    fflush(stdout);
    poly_point_init(&crs->cghat[i]);
    poly_elptic_mul(&crs->cghat[i], crs->ghat, cvec[i], crs->Ex);
    point_init(&crs->cg[i]);
    elptic_mul(&crs->cg[i], crs->g, cvec[i], crs->grp->E);
  }
  printf("\n");
}

/*  convert G1 point to G2 point.  
    Enter with space for G2 initalized.
*/

void tog2(POLY_POINT *G2, POINT G1)
{
  G2->x.deg = 0;
  mpz_set(G2->x.coef[0], G1.x);
  G2->y.deg = 0;
  mpz_set(G2->y.coef[0], G1.y);
}

