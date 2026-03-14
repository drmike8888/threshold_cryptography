/*****************************************************************
 *                                                               *
 *   Use Weil pairing to compute threshold encryption base on    *
 *   paper by Waters & Wu. Section 5.2                           *
 *                                                               *
 ****************************************************************/

#include "mpz_raw.h"
#include "pairing.h"
#include "thrshldencrpt_setup.h"
#include <string.h>

/* N = 2^n is the maximum number of users */

/* #define n (long)6 */
/* #define N (long)64 */
/* #define n (long)2 */
/* #define N (long)4 */
#define n (long)4
#define N (long) 16

/* parse parameter buffer. Skip comments and look for * = <value>.
   Set values only - can not compute anything until initializations
   completed. */

void get_group(GROUP *grp, char *bufr)
{
  int i;
  char *pos;

/* read in degree value */
  
  pos = strstr(bufr, "k = ");
  if(!pos)
  {
    printf("no degree in file\n");
    exit(-3);
  }
  pos += 4;
  sscanf(pos, "%ld", &grp->degree);

/* read in prime field value */
  
  pos = strstr(bufr, "q = ");
  if(!pos)
  {
    printf("no prime in file\n");
    exit(-4);
  }
  pos += 4;
  i = 20;
  while(pos[i] != '\n') i++;
  pos[i] = 0;
  mpz_init_set_str(grp->prm, pos, 10);
  pos[i] = '\n';
  
/* read in torsion value (group order) */
  
  pos = strstr(bufr, "r = ");
  if(!pos)
  {
    printf("no torsion in file\n");
    exit(-5);
  }
  pos += 4;
  i = 20;
  while(pos[i] != '\n') i++;
  pos[i] = 0;
  mpz_init_set_str(grp->tor, pos, 10);
  pos[i] = '\n';

/* read in t value (trace of Frobenius) */
  
  pos = strstr(bufr, "t = ");
  if(!pos)
  {
    printf("no t value in file\n");
    exit(-9);
  }
  pos += 4;
  i = 2;
  while(pos[i] != '\n') i++;
  pos[i] = 0;
  mpz_init_set_str(grp->t, pos, 10);
  pos[i] = '\n';

/* read in base curve cardinality */
  
  pos = strstr(bufr, "#E = ");
  if(!pos)
  {
    printf("no cardinality in file\n");
    exit(-6);
  }
  pos += 5;
  i = 20;
  while(pos[i] != '\n') i++;
  pos[i] = 0;
  mpz_init_set_str(grp->cardE, pos, 10);
  pos[i] = '\n';

/* read in base curve a4 and a6 parameters */
  
  pos = strstr(bufr, "a4 = ");
  if(!pos)
  {
    printf("no curve a4 parameter in file\n");
    exit(-7);
  }
  pos += 5;
  i = 20;
  while(pos[i] != '\n') i++;
  pos[i] = 0;
  mpz_init_set_str(grp->E.a4, pos, 10);
  pos[i] = '\n';
  pos = strstr(bufr, "a6 = ");
  if(!pos)
  {
    printf("no curve a6 parameter in file\n");
    exit(-8);
  }
  pos += 5;
  i = 20;
  while(pos[i] != '\n') i++;
  pos[i] = 0;
  mpz_init_set_str(grp->E.a6, pos, 10);
}

int main(int argc, char *argv[])
{
  FILE *prmdat, *crsdat;
  char *parmbuf;
  GROUP thrshldgroup;
  CRS crs;
  POLY irrd;
  mpz_t prm, c, *s_vec;
  long i, j, k, m, found, n1;
  POINT Tmp;
  POLY_POINT Tmp2;
  mpz_t *sharemtx, *cvec;
  
/* see if parameter file exits, read it in */
  
  if(argc < 2)
  {
    printf("Use: ./threshod_encrypt <file with parameters>\n");
    exit(-1);
  }
  prmdat = fopen(argv[1], "r");
  if(!prm)
  {
    printf("can't find file %s\n", argv[1]);
    exit(-2);
  }
  parmbuf = (char*)malloc(2048);
  i = 0;
  while(!feof(prmdat))
  {
    parmbuf[i] = fgetc(prmdat);
    i++;
  }
  fclose(prmdat);
  
/*  initialize base prime and polynomial  */

  get_group(&thrshldgroup, parmbuf);
  gmp_printf("q: %Zd\n", thrshldgroup.prm);
  gmp_printf("r: %Zd\n", thrshldgroup.tor);
  printf("k: %ld\n", thrshldgroup.degree);
  gmp_printf("t: %Zd\n", thrshldgroup.t);
  gmp_printf("#E: %Zd\n", thrshldgroup.cardE);
  gmp_printf("a4: %Zd\n", thrshldgroup.E.a4);
  gmp_printf("a6: %Zd\n", thrshldgroup.E.a6);

/* initialize random generator */

  mseed(594612);                      // these two routines
  minit(thrshldgroup.prm);            // a change from book code
  poly_init(&irrd);
  if(poly_irreducible(&irrd, thrshldgroup.degree))
    poly_printf("Found irreducible polynomial:\n", irrd);
  else
    printf("no irreducible polynomial found...\n");
  poly_irrd_set(irrd);
  poly_mulprep(irrd);

  
  poly_curve_init(&crs.Ex);    // this sets .deg to 0
  mpz_set(crs.Ex.a4.coef[0], thrshldgroup.E.a4);
  mpz_set(crs.Ex.a6.coef[0], thrshldgroup.E.a6);

/* compute cofactor of base curve  */
  
  mpz_inits(thrshldgroup.cardEx, thrshldgroup.cobse, thrshldgroup.coxtd, NULL);
  mpz_div(thrshldgroup.cobse, thrshldgroup.cardE, thrshldgroup.tor);
  
/*  compute #E extension curves  */

  cardinality(thrshldgroup.cardEx, thrshldgroup.t, thrshldgroup.degree);
//  gmp_printf("extension has %Zd points\n", thrshldgroup.cardEx);
  mpz_div(thrshldgroup.coxtd, thrshldgroup.cardEx, thrshldgroup.tor);
  mpz_div(thrshldgroup.coxtd, thrshldgroup.coxtd, thrshldgroup.tor);   // remove tor^2 for cofactor

/*  create generator points from random points and co-factors  */

  point_init(&crs.g);
  point_rand(&crs.g, thrshldgroup.E);
  elptic_mul(&crs.g, crs.g, thrshldgroup.cobse, thrshldgroup.E);
  point_printf("G1 generator: ", crs.g);

  printf("computing G2 generator\n");
  poly_point_init(&crs.ghat);
  poly_point_rand(&crs.ghat, crs.Ex);
  poly_elptic_mul(&crs.ghat, crs.ghat, thrshldgroup.coxtd, crs.Ex);
//  poly_point_printf("crs.ghat generator:\n", crs.ghat);

/* create reference point for Weil curve computation */

  printf("computing reference point S\n");
  poly_point_init(&thrshldgroup.S);
  found = 0;
  while(!found)
  {
    poly_point_rand(&thrshldgroup.S, crs.Ex);
    poly_elptic_mul(&thrshldgroup.S, thrshldgroup.S, thrshldgroup.tor, crs.Ex);
    if(poly_test_point(thrshldgroup.S))
      found = 0;
    else
      found = 1;
  }
//  poly_point_printf("S reference:\n", thrshldgroup.S);

/* begin common reference string creation */

  printf("computing u hat\n");
  crs.grp = &thrshldgroup;                    // point to pairing curve
  poly_point_init(&crs.uhat);
  poly_point_rand(&crs.uhat, crs.Ex);
  poly_elptic_mul(&crs.uhat, crs.uhat, thrshldgroup.coxtd, crs.Ex);
  printf("computing h hat\n");
  poly_point_init(&crs.hhat);
  poly_point_rand(&crs.hhat, crs.Ex);
  poly_elptic_mul(&crs.hhat, crs.hhat, thrshldgroup.coxtd, crs.Ex);

/* the values c and s_vec are "toxic waste", 
   should be securely generated and lost forever */

  mpz_init(c);
  mod_rand(c, thrshldgroup.tor);              // new function not in book
  s_vec = (mpz_t*)malloc(sizeof(mpz_t)*N);
  for(i=0; i<N; i++)
  {
    mpz_init(s_vec[i]);
    mod_rand(s_vec[i], thrshldgroup.tor);
  }

  prmdat = fopen("debug.svec", "w");
  mpz_out_raw(prmdat, c);
  for(i=0; i<N; i++)
    mpz_out_raw(prmdat, s_vec[i]);
  fclose(prmdat);
  
/* B = w(G1, G2)^s[0] = w(s[0]*G1, G2) */

  printf("computing B\n");
  point_init(&Tmp);
  elptic_mul(&Tmp, crs.g, s_vec[0], thrshldgroup.E);
  poly_point_init(&Tmp2);
  tog2(&Tmp2, Tmp);
  poly_init(&crs.B);
  weil(&crs.B, Tmp2, crs.ghat, thrshldgroup.S, thrshldgroup.tor, crs.Ex);

/* create share generation matrix */

  n1 = 2*N - 1;
  sharemtx = (mpz_t*)malloc(sizeof(mpz_t)*n1*N);
  for(i=0; i<n1; i++)
    for(j=0; j<N; j++)
      mpz_init(sharemtx[i*N + j]);
  genshare(sharemtx, N, thrshldgroup.tor);
  /* printf("source matrix: \n"); */
  /* for(i=0; i<n1; i++) */
  /* { */
  /*   for(j=0; j<N; j++) */
  /*     gmp_printf("%Zd ", sharemtx[i*N + j]); */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */

/* compute z0_hat and vl0_hat points of CRS */

  cvec = (mpz_t*)malloc(sizeof(mpz_t)*6*N);
  crs.vl0 = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n1);
  zvcalc(cvec, &crs, sharemtx, c, s_vec, N);

/* compute y_hat and tau_hat vector and matrix points of CRS */

  crs.yhat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n);
  crs.tauhat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n*(2*n1 + 1));
  ytaucalc(&crs, cvec, n, N);

/* compute c_vec[i]*g_hat for i in -2N+1 ... 2N-1 */

  crs.cghat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*4*N);
  crs.cg = (POINT*)malloc(sizeof(POINT)*4*N);
  cgcalc(&crs, cvec, N);

/* save CRS to disk */

  crsdat = fopen("threshold_setup.bin", "w");
  j = n;
  fwrite(&j, sizeof(long), 1,  crsdat);
  j = N;
  fwrite(&j, sizeof(long), 1,  crsdat);
  fwrite(&crs.grp->degree, sizeof(long), 1, crsdat);  // group data in order
  mpz_out_raw(crsdat, crs.grp->prm);
  mpz_out_raw(crsdat, crs.grp->tor);
  mpz_out_raw(crsdat, crs.grp->t);
  mpz_out_raw(crsdat, crs.grp->cardE);
  mpz_out_raw(crsdat, crs.grp->cardEx);
  mpz_out_raw(crsdat, crs.grp->cobse);
  mpz_out_raw(crsdat, crs.grp->coxtd);
  curve_write(&crs.grp->E, crsdat);
  poly_point_write(&crs.grp->S, crsdat);
  poly_curve_write(&crs.Ex, crsdat);                 // then rest of CRS
  point_write(&crs.g, crsdat);
  poly_point_write(&crs.ghat, crsdat);
  poly_write(&crs.B, crsdat);
  poly_point_write(&crs.uhat, crsdat);
  poly_point_write(&crs.hhat, crsdat);
  for(i=0; i<4*N-1; i++)
    poly_point_write(&crs.cghat[i], crsdat);
  for(i=0; i<4*N-1; i++)
    point_write(&crs.cg[i], crsdat);
  poly_point_write(&crs.z0hat, crsdat);
  for(i=0; i<n1; i++)
    poly_point_write(&crs.vl0[i], crsdat);
  for(i=0; i<n; i++)
    poly_point_write(&crs.yhat[i], crsdat);
  m = 4*N - 1;
  for(j=0; j<n; j++)
    for(i=0; i<m; i++)
      poly_point_write(&crs.tauhat[j*m + i], crsdat);
  fclose(crsdat);
  
/* clean up all variables to make sure memory was ok */

  mpz_clears(thrshldgroup.prm, thrshldgroup.tor, thrshldgroup.t, NULL);
  mpz_clears(thrshldgroup.cardE, thrshldgroup.E.a4, thrshldgroup.E.a6, NULL);
  mpz_clears(thrshldgroup.cardEx, thrshldgroup.cobse, thrshldgroup.coxtd, NULL);
  mpz_clear(c);
  for(i=0; i<N; i++)
    mpz_clear(s_vec[i]);
  free(s_vec);
  for(i=0; i<n1; i++)
    for(j=0; j<N; j++)
      mpz_clear(sharemtx[i*N + j]);
  free(sharemtx);
  poly_clear(&irrd);
  poly_clear(&crs.B);
  poly_point_clear(&crs.ghat);
  poly_point_clear(&thrshldgroup.S);
  poly_point_clear(&crs.uhat);
  poly_point_clear(&crs.hhat);
  poly_curve_clear(&crs.Ex);
  for(i=0; i<4*N; i++)
    poly_point_clear(&crs.cghat[i]);
  free(crs.cghat);
  point_clear(&Tmp);
  poly_point_clear(&Tmp2);
  for(i=0; i<6*N; i++)
    mpz_clear(cvec[i]);
  free(cvec);
  for(i=0; i<n1; i++)
    poly_point_clear(&crs.vl0[i]);
  free(crs.vl0);
  m = 2*n1 + 1;
  for(j=0; j<n; j++)
  {
    poly_point_clear(&crs.yhat[j]);
    for(i=0; i<m; i++)
      poly_point_clear(&crs.tauhat[j*m + i]);
  }
  free(crs.yhat);
  free(crs.tauhat);
}
