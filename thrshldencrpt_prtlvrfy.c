/*******************************************************
 *                                                     *
 *    Perform partial verify simulation. In reality,   *
 *    this is done while aggregating partial decrypt   *
 *    from each individual.  This just reads every     *
 *    set of data from one file.  Input CRS, public    *
 *    key file, cipher file, partial.sig file and      *
 *    S list file. does math check as described in     *
 *    paper.                                           *
 *                                                     *
 *                 Author = Mike Rosing                *
 *                  Date = 25 Feb. 2026                *
 *                                                     *
 ******************************************************/

#include "mpz_raw.h"
#include "eliptic.h"
#include "thrshldencrpt_setup.h"
#include "thrshldencrpt_crsin.h"
#include "pairing.h"

int main(int argc, char *argv[])
{
  FILE *key, *cphr, *slst;
  long i, N;
  long T, *S, slen;
  CRS crs;
  POLY *A;                      // vector
  POLY_POINT *vprm;             // matrix
  mpz_t tag, *alpha, *rl;
  POLY C1, left, right;
  POINT c2;
  POLY_POINT c3hat, c4hat, Tmp;
  POINT *sigma1;                 // vector
  POLY_POINT *sigma2;            // vector
  POLY_POINT G2;
  
/* make sure enough files listed on command line */

  if(argc < 5)
  {
    printf("Use: ./thrshldencrpt_prtlvrfy <CRS> <public key> <cipher data> <partial sig>\n");
    exit(-1);
  }
  
/* read in CRS data (or die) */

  if(crs_read(&crs, argv[1]) < 0)
    exit(-2);
  printf("CRS data read in\n");
  N = crs.N;
  mathinit(261793148, crs);

/* read in public keys (or die) */

  if(key_read(&A, &vprm, argv[2]) < 0)
    exit(-3);
  printf("public keys read in\n");

/* read in cipher text file 
   Each user should input their own tag which is compared to the
   cipher tag. This simulation ignores that step. */

  if(cipher_read(tag, &T, &C1, &c2, &c3hat, &c4hat, argv[3]) < 0)
    exit(-4);
  printf("cipher text file input\n");

/* read in partial sig data */

  cphr = fopen(argv[4], "r");
  if(!cphr)
  {
    printf("can't find cipher data file %s\n", argv[4]);
    exit(-5);
  }
  fread(&slen, sizeof(long), 1, cphr);
  S = (long*)malloc(sizeof(long)*slen);
  fread(S, sizeof(long), slen, cphr);
  sigma1 = (POINT*)malloc(sizeof(POINT)*slen);
  sigma2 = (POLY_POINT*)malloc(sizeof(POLY_POINT)*slen);
  for(i=0; i<slen; i++)
  {
    point_init(&sigma1[i]);
    point_read(&sigma1[i], cphr);
    poly_point_init(&sigma2[i]);
    poly_point_read(&sigma2[i], cphr);
  }
  fclose(cphr);

/* compute left hand side of verify test */

  poly_point_init(&Tmp);
  poly_elptic_mul(&Tmp, crs.uhat, tag, crs.Ex);
  poly_elptic_sum(&Tmp, Tmp, crs.hhat, crs.Ex);
  poly_init(&left);
  poly_init(&right);
  poly_point_init(&G2);
  for(i=0; i<slen; i++)
  {
    printf("compute left side formula for %ld\n", S[i]+1);
    tog2(&G2, sigma1[i]);
    weil(&left, G2, Tmp, crs.grp->S, crs.grp->tor, crs.Ex);
    poly_mul(&left, left, A[S[i]]);
    tog2(&G2, crs.g);
    weil(&right, G2, sigma2[i], crs.grp->S, crs.grp->tor, crs.Ex);
    if(poly_cmp(left, right))
      printf("Signer %ld verifies\n", S[i]+1);
    else
      printf("Signer %ld FAILS!\n", S[i]+1);
  }

/* clean up memory */

  poly_clear(&left);
  poly_clear(&right);
  poly_point_clear(&Tmp);
  for(i=0; i<slen; i++)
  {
    point_clear(&sigma1[i]);
    poly_point_clear(&sigma2[i]);
  }
  free(sigma1);
  free(sigma2);
  cipher_clear(&C1, &c2, &c3hat, &c4hat);
  key_clear(&A, &vprm, N);
  crs_clear(&crs);
}
