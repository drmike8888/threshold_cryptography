/**********************************************************
 *                                                        *
 *    Deal with CRS structure. create space and read in   *
 *    data from file and free all space when done.        *
 *                                                        *
 *                 Author = Mike Rosing                   *
 *                  Date = 19 Feb. 2026                   *
 *                                                        *
 *********************************************************/

#include <stdio.h>
#include "mpz_raw.h"
#include "thrshldencrpt_setup.h"

/* read in CRS data from file. Returns -1 if file not found,
   +1 otherwise.  initializes GMP.
*/

int crs_read(CRS *crs, char *filename)
{
  FILE *cdat;
  long i, j, m, n, N;

  cdat = fopen(filename, "r");
  if(!cdat)
  {
    printf("can't find CRS file %s\n", filename);
    return -1;
  }
  fread(&crs->n, sizeof(long), 1, cdat);
  fread(&crs->N, sizeof(long), 1, cdat);
  N = crs->N;
  n = crs->n;
  crs->grp = (GROUP*)malloc(sizeof(GROUP));
  fread(&crs->grp->degree, sizeof(long), 1, cdat);
  mpz_inits(crs->grp->prm, crs->grp->tor, crs->grp->t, crs->grp->cardE,
	    crs->grp->cardEx, crs->grp->cobse, crs->grp->coxtd, NULL);
  mpz_inp_raw(crs->grp->prm, cdat);
  mpz_inp_raw(crs->grp->tor, cdat);
  mpz_inp_raw(crs->grp->t, cdat);
  mpz_inp_raw(crs->grp->cardE, cdat);
  mpz_inp_raw(crs->grp->cardEx, cdat);
  mpz_inp_raw(crs->grp->cobse, cdat);
  mpz_inp_raw(crs->grp->coxtd, cdat);
  curve_init(&crs->grp->E);
  curve_read(&crs->grp->E, cdat);
  poly_point_init(&crs->grp->S);
  poly_point_read(&crs->grp->S, cdat);
  poly_curve_init(&crs->Ex);
  poly_curve_read(&crs->Ex, cdat);
  point_init(&crs->g);
  point_read(&crs->g, cdat);
  poly_point_init(&crs->ghat);
  poly_point_read(&crs->ghat, cdat);
  poly_init(&crs->B);
  poly_read(&crs->B, cdat);
  poly_point_init(&crs->uhat);
  poly_point_read(&crs->uhat, cdat);
  poly_point_init(&crs->hhat);
  poly_point_read(&crs->hhat, cdat);
  m = 4*N - 1;
  crs->cghat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*m);
  for(i=0; i<m; i++)
  {
    poly_point_init(&crs->cghat[i]);
    poly_point_read(&crs->cghat[i], cdat);
  }
  crs->cg = (POINT*)malloc(sizeof(POINT)*m);
  for(i=0; i<m; i++)
  {
    point_init(&crs->cg[i]);
    point_read(&crs->cg[i], cdat);
  }
  poly_point_init(&crs->z0hat);
  poly_point_read(&crs->z0hat, cdat);
  crs->vl0 = (POLY_POINT*)malloc(sizeof(POLY_POINT)*(2*N-1));
  for(i=0; i<2*N-1; i++)
  {
    poly_point_init(&crs->vl0[i]);
    poly_point_read(&crs->vl0[i], cdat);
  }
  crs->yhat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n);
  for(j=0; j<n; j++)
  {
    poly_point_init(&crs->yhat[j]);
    poly_point_read(&crs->yhat[j], cdat);
  }
  crs->tauhat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n*m);
  for(j=0; j<n; j++)
  {
    for(i=0; i<m; i++)
    {
      poly_point_init(&crs->tauhat[j*m + i]);
      poly_point_read(&crs->tauhat[j*m + i], cdat);
    }
  }
  fclose(cdat);
  return 1;
}

/* clear out all CRS data from memory */

void crs_clear(CRS *crs)
{
  long i, j, m, n, N;

  n = crs->n;
  N = crs->N;
  m = 4*N - 1;
  mpz_clears(crs->grp->prm, crs->grp->tor, crs->grp->t, crs->grp->cardE,
	    crs->grp->cardEx, crs->grp->cobse, crs->grp->coxtd, NULL);
  curve_clear(&crs->grp->E);
  poly_point_clear(&crs->grp->S);
  free(crs->grp);
  poly_curve_clear(&crs->Ex);
  point_clear(&crs->g);
  poly_point_clear(&crs->ghat);
  poly_clear(&crs->B);
  poly_point_clear(&crs->uhat);
  poly_point_clear(&crs->hhat);
  for(i=0; i<m; i++)
    poly_point_clear(&crs->cghat[i]);
  free(crs->cghat);
  for(i=0; i<m; i++)
    point_clear(&crs->cg[i]);
  free(crs->cg);
  for(i=0; i<2*N-1; i++)
    poly_point_clear(&crs->vl0[i]);
  free(crs->vl0);
  for(j=0; j<n; j++)
    poly_point_clear(&crs->yhat[j]);
  free(crs->yhat);
  for(j=0; j<n; j++)
    for(i=0; i<m; i++)
      poly_point_clear(&crs->tauhat[j*m + i]);
  free(crs->tauhat);
}

/* initialize GMP math system, irreducible polynomials
   and random number generator.  */

void mathinit(unsigned long seed, CRS crs)
{
  POLY irrd;
  
  mseed(seed);                        // these two routines
  minit(crs.grp->prm);                // a change from book code
  poly_init(&irrd);
  if(poly_irreducible(&irrd, crs.grp->degree))
    poly_printf("Found irreducible polynomial:\n", irrd);
  else
    printf("no irreducible polynomial found...\n");
  poly_irrd_set(irrd);
  poly_mulprep(irrd);
}

/*  input public key file data from disk. */

int key_read(POLY **A, POLY_POINT **vprm, char *filename)
{
  FILE *keys;
  int j, k, N, m;

  keys = fopen(filename, "r");
  if(!keys)
  {
    printf("can't find public key file %s\n", filename);
    return -1;
  }
  fread(&N, sizeof(long), 1, keys);
  *A = (POLY*)malloc(sizeof(POLY)*N);
  m = 4*N - 1;
  *vprm = (POLY_POINT*)malloc(sizeof(POLY_POINT)*N*m);
  for(j=0; j<N; j++)
  {
    poly_init(&(*A)[j]);
    poly_read(&(*A)[j], keys);
    for(k=0; k<m; k++)
    {
      poly_point_init(&(*vprm)[j*m + k]);
      poly_point_read(&(*vprm)[j*m + k], keys);
    }
  }
  fclose(keys);
  return 1;
}

/* clear memory of public key data */

void key_clear(POLY **A, POLY_POINT **vprm, long N)
{
  long j, k, m;

  m = 4*N - 1;
  for(j=0; j<N; j++)
  {
    poly_clear(&(*A)[j]);
    for(k=0; k<m; k++)
      poly_point_clear(&(*vprm)[j*m + k]);
  }
  free(*A);
  free(*vprm);
}

/* read in encryption key. 1 is success, -1 is failure */

int encrptkey_read(long *N, long *k, long **L, POLY_POINT *zhat,
		   POLY_POINT **vhat, char *filename)
{
  FILE *keys;
  long j, n1;
  
  keys = fopen(filename, "r");
  if(!keys)
  {
    printf("can't find encryption key file %s\n", filename);
    return -1;
  }
  fread(N, sizeof(long), 1, keys);
  fread(k, sizeof(long), 1, keys);
  *L = (long*)malloc(sizeof(long)*(*k));
  fread(*L, sizeof(long), *k, keys);
  poly_point_init(zhat);
  poly_point_read(zhat, keys);
  n1 = 2*(*N) - 1;
  *vhat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n1);
  for(j=0; j<n1; j++)
  {
    poly_point_init(&(*vhat)[j]);
    poly_point_read(&(*vhat)[j], keys);
  }
  fclose(keys);
  return 1;
}

void encrptkey_clear(long N, long **L, POLY_POINT *zhat, POLY_POINT **vhat)
{
  long j, n1;

  free(*L);
  n1 = 2*N - 1;
  poly_point_clear(zhat);
  for(j=0; j<n1; j++)
    poly_point_clear(&(*vhat)[j]);
  free(*vhat);
}

/*  read in cipher text file 
    return +1 on success, -1 on failure  */

int cipher_read(mpz_t tag, long *T, POLY *C1, POINT *c2, POLY_POINT *c3hat,
		POLY_POINT *c4hat, char *filename)
{
  FILE *cphr;

  cphr = fopen(filename, "r");
  if(!cphr)
  {
    printf("can't find cipher text file %s\n", filename);
    return -1;
  }
  mpz_init(tag);
  mpz_inp_raw(tag, cphr);
  fread(T, sizeof(long), 1, cphr);
  poly_init(C1);
  poly_read(C1, cphr);
  point_init(c2);
  point_read(c2, cphr);
  poly_point_init(c3hat);
  poly_point_read(c3hat, cphr);
  poly_point_init(c4hat);
  poly_point_read(c4hat, cphr);
  return 1;
}

void cipher_clear(POLY *C1, POINT *c2, POLY_POINT *c3hat,
		  POLY_POINT *c4hat)
{
  poly_clear(C1);
  point_clear(c2);
  poly_point_clear(c3hat);
  poly_point_clear(c4hat);
}
