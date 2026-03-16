/***************************************************************
 *                                                             *
 *   Threshold encryption decryption final step. Enter with    *
 *   CRS, encryption key, cipher data, S list, and partial     *
 *   signatures. Outputs computation that should give message  *
 *   polynomial value.                                         *
 *                                                             *
 *                      Author = Mike Rosing                   *
 *                       Date = 26 Feb. 2026                   *
 *                                                             *
 **************************************************************/

#include "mpz_raw.h"
#include "eliptic.h"
#include "thrshldencrpt_setup.h"
#include "thrshldencrpt_crsin.h"
#include "pairing.h"

int main(int argc, char *argv[])
{
  FILE *key, *cphr, *slst;
  long i, j, k, *L, N, n, lmt, n1;
  long T, *S, slen, *X, m, mask, n4;
  long mdex, sdex, ldex, *Spad, *Tdex;
  CRS crs;
  POLY *A;                      // vector
  POLY_POINT *vprm;             // matrix
  POLY_POINT zhat;
  POLY_POINT *vhat;             // vector
  mpz_t tag, *share, *Mprm;     // share & Mprm are matricies
  mpz_t *inv;                   // inverse of M'
  mpz_t *omegat;                // vector
  POLY C1, weiltop, weilbot1, weilbot2;
  POINT c2;
  POLY_POINT c3hat, c4hat, Tmp;
  POINT *sigma1;                 // vector
  POLY_POINT *sigma2;            // vector
  POLY_POINT G2;
  POINT sigagg1, sigagg3, tmul;
  POLY_POINT sigagg21, sigagg22, sigagg23, sigagg2;
  POLY rcrvdmsg;
  
/* make sure enough files listed on command line */

  if(argc < 5)
  {
    printf("Use: ./thrshldencrpt_dcrypt <CRS> <encryption key> <cipher data> <partial sig>\n");
    exit(-1);
  }
  
/* read in CRS data (or die) */

  if(crs_read(&crs, argv[1]) < 0)
    exit(-2);
  printf("CRS data read in\n");
  N = crs.N;
  n = crs.n;
  mathinit(261793148, crs);

/* read in encryption key (k is length of L) */

  if(encrptkey_read(&N, &k, &L, &zhat, &vhat, argv[2]) < 0)
    exit(-3);
  printf("encryption key read in\n");
  
/* read in cipher text file  */

  if(cipher_read(tag, &T, &C1, &c2, &c3hat, &c4hat, argv[3]) < 0)
    exit(-4);
  printf("cipher text file input\n");

/* read in partial sig data */

  cphr = fopen(argv[4], "r");
  if(!cphr)
  {
    printf("can't find cipher data file %s\n", argv[4]);
    exit(-7);
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
  printf("partial sig data file read in\n");
  
/* create values used to compute decryption function */  

  lmt = k - T;
  printf("lmt: %ld\n", lmt);
  if(lmt < 0)
  {
    printf("threshold %ld larger than encryptor list %ld\n", T, k);
    exit(-8);
  }
  n1 = 2*N - 1;
  share = (mpz_t*)malloc(sizeof(mpz_t)*N*n1);
  for(i=0; i<n1; i++)
    for(j=0; j<N; j++)
      mpz_init(share[i*N + j]);
  genshare(share, N, crs.grp->tor);
  printf("share matrix computed\n");
	 
/* create index list for each bit in L - T */
  
  X = (long*)malloc(sizeof(long)*(n + 1));
  X[0] = 1;
  for(i=1; i<=n; i++)
    X[i] = 2*X[i - 1];
  
/* pull out N rows from share matrix that are used for decryption */

  Mprm = (mpz_t*)malloc(sizeof(mpz_t)*N*N);
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      mpz_init(Mprm[i*N + j]);

/* THIS ASSUMES S LIST AND L LIST ARE IN NUMERICAL ORDER! */
  
  mdex = 0;
  sdex = 0;
  ldex = 0;
  Spad = (long*)malloc(sizeof(long)*N);
  Tdex = (long*)malloc(sizeof(long)*slen);
  for(i=0; i<N; i++)
  {
    if((i == S[sdex]) && (sdex < T))
    {
      Spad[mdex] = i;
      Tdex[sdex] = mdex;
      for(j=0; j<N; j++)
	mpz_set(Mprm[mdex*N + j], share[i*N + j]);
      if(S[sdex] == L[ldex])
	ldex++;
      mdex++;
      sdex++;
    }
    else if(i != L[ldex])
    {
      Spad[mdex] = i;
      for(j=0; j<N; j++)
	mpz_set(Mprm[mdex*N + j], share[i*N + j]);
      mdex++;
    }
    else
      ldex++;
  }
  printf("mdex: %ld lmt: %ld\n", mdex, lmt);
  mask = 1;
  for(j=0; j<n; j++)
  {
    if(mask & lmt)
    {
      i = X[j];
      while(i < X[j + 1])
      {
	Spad[mdex] = N + i - 1;
	for(m=0; m<N; m++)
	  mpz_set(Mprm[mdex*N + m], share[Spad[mdex]*N + m]);
	mdex++;
	i++;
      }
    }
    mask <<= 1;
  }
  printf("check mdex %ld = N %ld\n", mdex, N);
  for(i=0; i<slen; i++)
    printf("Tdex[%ld] = %ld\n", i, Tdex[i]);
  for(i=0; i<N; i++)
    printf("Spad[%ld] = %ld\n", i, Spad[i]);

/* compute matrix inverse of square M' matrix (modulo torsion value) */

  inv = (mpz_t*)malloc(sizeof(mpz_t)*N*N);
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      mpz_init(inv[i*N + j]);
  if(mod_matinv(inv, N, Mprm, crs.grp->tor) < 0)
  {
    printf("M' matrix not invertable\n");
    exit(-7);
  }

/* first row of inverse matrix is reconstruction vector */
  
  omegat =(mpz_t*)malloc(sizeof(mpz_t)*N);
  for(i=0; i<N; i++)
    mpz_init_set(omegat[i], inv[i]);

/* compute decryption components
   first is point sigma_agg_1. S list order, not user number.  */

  printf("computing sigma agg 1\n");
  point_init(&sigagg1);
  point_init(&tmul);
  for(i=0; i<T; i++)
  {
    elptic_mul(&tmul, sigma1[i], omegat[Tdex[i]], crs.grp->E);
    elptic_sum(&sigagg1, sigagg1, tmul, crs.grp->E);
  }
//  point_printf("sigagg1:\n", sigagg1);
  
  printf("computing sigma agg 2 part 1\n");
  poly_point_init(&sigagg21);
  for(i=0; i<T; i++)
  {
    poly_elptic_mul(&Tmp, sigma2[i], omegat[Tdex[i]], crs.Ex);
    poly_elptic_sum(&sigagg21, sigagg21, Tmp, crs.Ex);
  }
  printf("computing sigma agg 2 part 2\n");
  poly_point_init(&sigagg22);
  for(i=0; i<N; i++)
  {
    poly_elptic_mul(&Tmp, vhat[Spad[i]], omegat[i], crs.Ex);
    poly_elptic_sum(&sigagg22, sigagg22, Tmp, crs.Ex);
  }
  printf("computing sigma agg 2 part 3\n");
  poly_point_init(&sigagg23);
  mask = 1;
  n4 = 4*N - 1;
  m = 2*N;
  for(j=0; j<n; j++)
  {
    if(!(mask & lmt))      // only clear bits allowed in sum
    {
      for(i=0; i<N; i++)
      {
	poly_elptic_mul(&Tmp, crs.tauhat[j*n4 + m + Spad[i]], omegat[i], crs.Ex);
	poly_elptic_sum(&sigagg23, sigagg23, Tmp, crs.Ex);
      }
    }
    mask <<= 1;
  }
  printf("computing sigma agg 2\n");
  poly_point_init(&sigagg2);
  poly_elptic_sum(&sigagg2, sigagg21, sigagg22, crs.Ex);
  poly_elptic_sum(&sigagg2, sigagg2, sigagg23, crs.Ex);
//  poly_point_printf("sigagg2:\n", sigagg2);

  printf("computing sigma agg 3\n");
  point_init(&sigagg3);
  m = 2*N - 2;
  for(i=0; i<N; i++)
  {
    elptic_mul(&tmul, crs.cg[m - Spad[i]], omegat[i], crs.grp->E);
    elptic_sum(&sigagg3, sigagg3, tmul, crs.grp->E);
  }
  
/* compute numerator of output */

  printf("computing numerator\n");
  poly_init(&weiltop);
  poly_point_init(&G2);
  tog2(&G2, c2);
  weil(&weiltop, G2, sigagg2, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_mul(&weiltop, weiltop, C1);

/* and denominator */

  printf("computing denominator\n");
  poly_init(&weilbot1);
  tog2(&G2, sigagg1);
  weil(&weilbot1, G2, c3hat, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_init(&weilbot2);
  tog2(&G2, sigagg3);
  weil(&weilbot2, G2, c4hat, crs.grp->S, crs.grp->tor, crs.Ex);
  
/* recover message data */

  poly_init(&rcrvdmsg);
  poly_div(&rcrvdmsg, weiltop, weilbot1);
  poly_div(&rcrvdmsg, rcrvdmsg, weilbot2);

  poly_printf("recovered message: \n", rcrvdmsg);
  
/* clear out memory */

  crs_clear(&crs);
  encrptkey_clear(N, &L, &zhat, &vhat);
  cipher_clear(&C1, &c2, &c3hat, &c4hat);
  free(S);
  for(i=0; i<slen; i++)
  {
    point_clear(&sigma1[i]);
    poly_point_clear(&sigma2[i]);
  }
  free(sigma1);
  free(sigma2);
  for(i=0; i<n1; i++)
    for(j=0; j<N; j++)
      mpz_clear(share[i*N + j]);
  free(share);
  free(X);
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      mpz_clear(Mprm[i*N + j]);
  free(Mprm);
  free(Spad);
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      mpz_clear(inv[i*N + j]);
  free(inv);
  for(i=0; i<N; i++)
    mpz_clear(omegat[i]);
  free(omegat);
  point_clear(&sigagg1);
  point_clear(&tmul);
  point_clear(&sigagg3);
  poly_point_clear(&sigagg21);
  poly_point_clear(&sigagg22);
  poly_point_clear(&sigagg23);
  poly_point_clear(&sigagg2);
  poly_clear(&weiltop);
  poly_point_clear(&G2);
  poly_clear(&weilbot1);
  poly_clear(&weilbot2);
  poly_clear(&rcrvdmsg);
}  
  
