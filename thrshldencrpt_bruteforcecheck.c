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
 *  26 Mar. 4 modified to use "secret" values directly to      *
 *  verify the math does what I think it is supposed to do.    *
 **************************************************************/

#include "mpz_raw.h"
#include "eliptic.h"
#include "thrshldencrpt_setup.h"
#include "thrshldencrpt_crsin.h"
#include "pairing.h"

int mod_matinv(mpz_t *inv, long n, mpz_t *mat, mpz_t prime);
int mod_matmul(mpz_t *c, mpz_t *a, long m, long n,
	       mpz_t *b, long r, long s, mpz_t prime);

void mdots(mpz_t ms, long row, mpz_t *s_vec, mpz_t *shrmtrx, long N, mpz_t mod);

int main(int argc, char *argv[])
{
  FILE *key, *cphr, *slst;
  long i, j, k, *L, N, n, lmt, n1;
  long T, *S, slen, *X, m, mask, n4;
  long mdex, sdex, ldex, *Spad;//, *Srow;
  CRS crs;
  POLY *A;                      // vector
  POLY_POINT *vprm;             // matrix
  POLY_POINT zhat;
  POLY_POINT *vhat;             // vector
  mpz_t tag, *share, *Mprm;     // share & Mprm are matricies
  mpz_t *inv;                   // inverse of M'
  mpz_t *omegat;                // vector
  POLY C1, tatetop, tatebot1, tatebot2;
  POINT c2;
  POLY_POINT c3hat, c4hat, Tmp;
  POINT *sigma1, *sig1chk;                 // vector
  POLY_POINT *sigma2, *sig2chk;            // vector
  POLY_POINT G2;
  POINT sigagg1, sigagg3, tmul;
  POLY_POINT sigagg21, sigagg22, sigagg23, sigagg2;
  POLY rcrvdmsg;
  // DEBUGGING
  mpz_t *rl, smtst1, smtst2, *svec, c, cpow, cpm;
  POLY tstrght, tstlft, taugg;
  mpz_t sigchk, *gma, t, *msvec, *alpha, ta, tb, tc, td;
  mpz_t tb1, tb2, tb3, tc1, tc2;
  POINT omc;
  POLY_POINT c4chk, tauuh, *vl0chk, *vprmchk, *tauchk;
  long n3, xpnt;
  
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
  poly_init(&tstrght);
  poly_point_init(&G2);
  tog2(&G2, crs.g);
  weil(&tstrght, G2, crs.ghat, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_printf("weil g, ghat):\n", tstrght);
  exit(0);
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
  
/* get t value */
  
  key = fopen("debug.t", "r");
  mpz_init(t);
  mpz_inp_raw(t, key);
  fclose(key);

/* get rl values */
  
  key = fopen("debug.rl", "r");
  rl = (mpz_t*)malloc(sizeof(mpz_t)*slen);
  for(i=0; i<slen; i++)
  {
    mpz_init(rl[i]);
    mpz_inp_raw(rl[i], key);
  }
  fclose(key);

/* get gamma_k values */
  
  key = fopen("debug.gamma", "r");
  gma = (mpz_t*)malloc(sizeof(mpz_t)*N);
  for(i=0; i<N-1; i++)
  {
    mpz_init(gma[i]);
    mpz_inp_raw(gma[i], key);
  }
  fclose(key);

/* get s vector values */
  
  key = fopen("debug.svec", "r");
  mpz_init(c);
  mpz_inp_raw(c, key);
  svec = (mpz_t*)malloc(sizeof(mpz_t)*N);
  for(i=0; i<N; i++)
  {
    mpz_init(svec[i]);
    mpz_inp_raw(svec[i], key);
  }
  fclose(key);

/* read in alpha key values */
  
  key = fopen("secret_keys.bin", "r");
  fread(&i, sizeof(long), 1, key);  // ignore, already known
  alpha = (mpz_t*)malloc(sizeof(mpz_t)*N);
  for(i=0; i<N; i++)
  {
    mpz_init(alpha[i]);
    mpz_inp_raw(alpha[i], key);
  }
  fclose(key);

/* read in public keys (or die) */

  if(key_read(&A, &vprm, "public_keys.bin") < 0)
    exit(-3);
  printf("public keys read in\n");

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

  for(i=0; i<k; i++)
    printf("L[%ld]: %ld\n", i, L[i]);
  for(i=0; i<slen; i++)
    printf("T[%ld]: %ld\n", i, S[i]);
  mdex = 0;
  sdex = 0;
  ldex = 0;
  Spad = (long*)malloc(sizeof(long)*N);
  /* Srow = (long*)malloc(sizeof(long)*slen); */
  for(i=0; i<N; i++)
  {
    if((i == S[sdex]) && (sdex < T))
    {
      Spad[mdex] = i;
      /* Srow[sdex] = i; */
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

/********************************************************************************/
/* brute force compute saved terms from secret values
   to check they were computed correctly. */

  /* poly_point_init(&tauuh); */
  /* poly_elptic_mul(&tauuh, crs.uhat, tag, crs.Ex); */
  /* poly_elptic_sum(&tauuh, tauuh, crs.hhat, crs.Ex); */
  /* sig1chk = (POINT*)malloc(sizeof(POINT)*slen); */
  /* sig2chk = (POLY_POINT*)malloc(sizeof(POLY_POINT)*slen); */
  poly_point_init(&Tmp);
  /* for(i=0; i<slen; i++) */
  /* { */
  /*   point_init(&sig1chk[i]); */
  /*   elptic_mul(&sig1chk[i], crs.g, rl[i], crs.grp->E); */
  /*   poly_point_init(&sig2chk[i]); */
  /*   poly_elptic_mul(&sig2chk[i], tauuh, rl[i], crs.Ex); */
  /*   poly_elptic_mul(&Tmp, crs.ghat, alpha[i+1], crs.Ex); */
  /*   poly_elptic_sum(&sig2chk[i], sig2chk[i], Tmp, crs.Ex); */
  /* } */
  
  /* for(i=0; i<slen; i++) */
  /* { */
  /*   printf("file point %ld", i); */
  /*   point_printf(" \n", sigma1[i]); */
  /*   printf("test point %ld", i);        sigma1, sigma2 checks*/
  /*   point_printf(" \n", sig1chk[i]); */
  /* } */
  /* for(i=0; i<slen; i++) */
  /* { */
  /*   printf("file point %ld", i); */
  /*   poly_point_printf("\n", sigma2[i]); */
  /*   printf("test point %ld", i); */
  /*   poly_point_printf("\n", sig2chk[i]); */
  /* } */

// z0 hat - checks OK

  msvec = (mpz_t*)malloc(sizeof(mpz_t)*n1);
  mpz_inits(smtst1, smtst2, cpow, sigchk, cpm, NULL);
  mpz_set(cpow, c);
  mpz_set_ui(smtst1, 0);
  for(i=0; i<n1; i++)
  {                         // compute z0hat terms
    mpz_init(msvec[i]);
    mdots(msvec[i], i, svec, share, N, crs.grp->tor);
    mod_mul(sigchk, msvec[i], cpow, crs.grp->tor);
    mod_add(smtst1, smtst1, sigchk, crs.grp->tor);
    mod_mul(cpow, cpow, c, crs.grp->tor);
  }
//  gmp_printf("z0: %Zd\n", smtst1);
// smtst1 has z0_hat coefficient
  /* poly_elptic_mul(&Tmp, crs.ghat, smtst1, crs.Ex); */
  /* poly_point_printf("z0 from file:\n", crs.z0hat); */
  /* poly_point_printf("z0 computed:\n", Tmp); */

  /* mpz_set(cpow, c); */  // z_hat checks ok
  /* mpz_set_ui(smtst2, 0); */
  /* for(i=0; i<k; i++)        // compute v' terms */
  /* {                           */
  /*   mod_mul(sigchk, alpha[i], cpow, crs.grp->tor); */
  /*   mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   mod_add(smtst2, smtst2, sigchk, crs.grp->tor); */
  /* } */
  // smtst2 has sum v' coefficient
  /* mod_add(sigchk, smtst1, smtst2, crs.grp->tor); */
  /* poly_elptic_mul(&Tmp, crs.ghat, sigchk, crs.Ex); */
  /* poly_point_printf("z hat from file:\n", zhat); */
  /* poly_point_printf("z hat computed:\n", Tmp); */

// check y hat terms OK

  /* mpz_powm_ui(cpow, c, 5, crs.grp->tor); */
  /* mod_mul(sigchk, gma[0], cpow, crs.grp->tor); */
  /* poly_elptic_mul(&Tmp, crs.ghat, sigchk, crs.Ex); */
  /* poly_point_printf("y0 from file:\n", crs.yhat[0]); */
  /* poly_point_printf("y0 computed:\n", Tmp); */
  /* mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* mod_mul(smtst1, gma[1], cpow, crs.grp->tor); */
  /* mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* mod_mul(sigchk, gma[2], cpow, crs.grp->tor); */
  /* mod_add(sigchk, sigchk, smtst1, crs.grp->tor); */
  /* poly_elptic_mul(&Tmp, crs.ghat, sigchk, crs.Ex); */
  /* poly_point_printf("y1 from file:\n", crs.yhat[1]); */
  /* poly_point_printf("y1 computed:\n", Tmp); */

// c3, c4 hat check  OK

  /* mpz_inits(tb1, tb2, tb3, NULL); */
  /* mpz_set(cpow, c); */
  /* mpz_set_ui(smtst1, 0); */
  /* for(i=0; i<n1; i++) */
  /* { */
  /*   mod_mul(sigchk, msvec[i], cpow, crs.grp->tor); */
  /*   mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*   mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* } */
  /* mpz_set(tb1, smtst1); */
  /* gmp_printf("b1: %Zd\n", smtst1); */
  /* mpz_set(cpow, c); */
  /* mpz_set_ui(smtst2, 0); */
  /* for(i=0; i<k; i++)        // compute v' terms */
  /* { */
  /*   mod_mul(sigchk, alpha[i], cpow, crs.grp->tor); */
  /*   mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   mod_add(smtst2, smtst2, sigchk, crs.grp->tor); */
  /* } */
  /* mpz_set(tb2, smtst2); */
  /* mod_add(smtst1, smtst1, smtst2, crs.grp->tor); */
  /* mpz_powm_ui(cpow, c, N+2, crs.grp->tor); */
  /* mod_mul(smtst2, gma[1], cpow, crs.grp->tor); */
  /* mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* mod_mul(sigchk, gma[2], cpow, crs.grp->tor); */
  /* mod_add(smtst2, smtst2, sigchk, crs.grp->tor); */
  /* mpz_set(tb3, smtst2); */
  /* mod_add(smtst1, smtst1, smtst2, crs.grp->tor); */
  /* mpz_inits(ta, tb, tc, td, NULL); */
  /* mod_mul(tb, smtst1, t, crs.grp->tor); */
  
  /* poly_point_init(&c4chk); */
  /* poly_elptic_mul(&c4chk, crs.ghat, smtst1, crs.Ex); */
  /* poly_point_printf("c4 from file:\n", c4hat); */
  /* poly_point_printf("c4 computed:\n", c4chk); */
  /* exit(0); */
// check vl0 OK

  /* vl0chk = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n1); */
  /* for(m=0; m<n1; m++) */
  /* { */
  /*   mpz_set_ui(smtst1, 0); */
  /*   mpz_powm_ui(cpm, c, m, crs.grp->tor); */
  /*   for(i=0; i<n1; i++) */
  /*   { */
  /*     if(i == m) */
  /* 	continue; */
  /*     mpz_powm_ui(cpow, c, i, crs.grp->tor); */
  /*     mod_div(cpow, cpow, cpm, crs.grp->tor); */
  /*     mod_mul(sigchk, msvec[m], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*   } */
  /*   poly_point_init(&vl0chk[m]); */
  /*   poly_elptic_mul(&vl0chk[m], crs.ghat, smtst1, crs.Ex); */
  /* } */

  /* for(m=0; m<n1; m++) */
  /* { */
  /*   printf("vl0 %ld ", m); */
  /*   if(poly_point_cmp(vl0chk[m], crs.vl0[m])) */
  /*     printf("matches\n"); */
  /*   else */
  /*     printf("fails\n"); */
  /* } */

// check v' matrix OK

  n4 = 4*N - 1;
  /* vprmchk = (POLY_POINT*)malloc(sizeof(POLY_POINT)*N*n4); */
  /* for(i=0; i<N; i++) */
  /* { */
  /*   mpz_powm_ui(cpow, c, 2*N-1, crs.grp->tor); */
  /*   mpz_invert(cpow, cpow, crs.grp->tor); */
  /*   for(m=0; m<j; m++) */
  /*   { */
  /*     printf("v prm (%ld, %ld)\r", i, m); */
  /*     fflush(stdout); */
  /*     mod_mul(sigchk, alpha[i], cpow, crs.grp->tor); */
  /*     poly_point_init(&vprmchk[j*i + m]); */
  /*     poly_elptic_mul(&vprmchk[j*i + m], crs.ghat, sigchk, crs.Ex); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   } */
  /* } */
  /* printf("\n"); */
  /* for(i=0; i<N; i++) */
  /* { */
  /*   for(m=0; m<j; m++) */
  /*   { */
  /*     printf("%ld %ld: ", i, m); */
  /*     fflush(stdout); */
  /*     if(poly_point_cmp(vprmchk[j*i + m], vprm[j*i + m])) */
  /* 	printf("matches\n"); */
  /*     else */
  /* 	printf("fails!\n"); */
  /*   } */
  /* } */

// check v hat  OK

  /* vl0chk = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n1); */
  /* for(m=0; m<n1; m++) */
  /* { */
  /*   mpz_set_ui(smtst1, 0); */
  /*   mpz_powm_ui(cpm, c, m, crs.grp->tor); */
  /*   for(i=0; i<k; i++) */
  /*   { */
  /*     if(i == m) */
  /* 	continue; */
  /*     mpz_powm_ui(cpow, c, i, crs.grp->tor); */
  /*     mod_div(cpow, cpm, cpow, crs.grp->tor); */
  /*     mod_mul(sigchk, alpha[i], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*   } */
  /*   poly_point_init(&vl0chk[m]); // extra term in v-hat */
  /*   poly_elptic_mul(&vl0chk[m], crs.ghat, smtst1, crs.Ex); */
  /*   poly_elptic_sum(&vl0chk[m], vl0chk[m], crs.vl0[m], crs.Ex); */
  /* } */
  /* for(m=0; m<n1; m++) */
  /* { */
  /*   printf("vl0 %ld ", m); */
  /*   if(poly_point_cmp(vl0chk[m], crs.vl0[m])) */
  /*     printf("matches\n"); */
  /*   else */
  /*     printf("fails\n"); */
  /* } */
  
  /* poly_point_printf("computed vhat[0]:\n", vl0chk[0]); */
  /* mod_div(sigchk, alpha[1], c, crs.grp->tor); */
  /* mod_div(smtst1, alpha[2], c, crs.grp->tor); */
  /* mod_div(smtst1, smtst1, c, crs.grp->tor); */
  /* mod_add(smtst1, sigchk, smtst1, crs.grp->tor); */
  /*   poly_elptic_mul(&vl0chk[m], crs.ghat, smtst1, crs.Ex); */
  /*   poly_elptic_sum(&vl0chk[m], vl0chk[m], crs.vl0[m], crs.Ex); */
  /* poly_point_printf("brute force vhat[0]:\n", vl0chk[0]); */
  /* exit(0); */
  
  /* for(m=0; m<n1; m++) */
  /* { */
  /*   printf("vhat index %ld\n", m); */
  /*   poly_point_printf("from file:\n", vhat[m]); */
  /*   poly_point_printf("computed:\n", vl0chk[m]); */
  /* } */

// check tau hat matrix  OK
  /* tauchk = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n*n4); */
  /* for(j=0; j<n; j++) */
  /* { */
  /*   i = X[j]; */
  /*   n3 = N + i; */
  /*   mpz_powm_ui(cpow, c, n3, crs.grp->tor); */
  /*   mpz_set_ui(smtst1, 0); */
  /*   while(i < X[j + 1]) */
  /*   { */
  /*     mod_mul(sigchk, cpow, gma[i - 1], crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*     i++; */
  /*   } */
  /*   n3 = 2*N - 1; */
  /*   mpz_powm_ui(cpow, c, n3, crs.grp->tor); */
  /*   for(m=0; m<n4; m++) */
  /*   { */
  /*     printf("tau %ld %ld\r", j, m); */
  /*     fflush(stdout); */
  /*     if((m < 3*N + X[j] - 1) || (m >= 3*N + X[j+1] - 1)) */
  /*     { */
  /* 	mod_mul(sigchk, smtst1, cpow, crs.grp->tor); */
  /* 	poly_point_init(&tauchk[n4*j + m]); */
  /* 	poly_elptic_mul(&tauchk[n4*j + m], crs.ghat, sigchk, crs.Ex); */
  /* 	mod_div(cpow, cpow, c, crs.grp->tor); */
  /*     } */
  /*     else */
  /*     { */
  /* 	poly_point_init(&tauchk[n4*j + m]); */
  /* 	mod_div(cpow, cpow, c, crs.grp->tor); */
  /*     } */
  /*   } */
  /* } */
  /* printf("taucheck\n"); */
  /* for(j=0; j<n; j++) */
  /* { */
  /*   for(m=0; m<n4; m++) */
  /*   { */
  /*     if(poly_test_point(crs.tauhat[j*n4 + m])) */
  /* 	printf("at infinity "); */
  /*     printf("%ld %ld: ", j, m); */
  /*     fflush(stdout); */
  /*     if(poly_point_cmp(tauchk[j*n4 + m], crs.tauhat[j*n4 + m])) */
  /* 	printf("matches\n"); */
  /*     else */
  /* 	printf("fails!\n"); */
  /*     poly_point_printf("\n", tauchk[j*n4 + m]); */
  /*   } */
  /* } */

// check sigagg2 part 2 (get msvec from above)

  mpz_inits(tc1, tc2, NULL);
  /* mpz_set_ui(tc, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_powm_ui(sigchk, c, Spad[m], crs.grp->tor); */
  /*   mpz_invert(cpow, sigchk, crs.grp->tor); */
  /*   mpz_set_ui(smtst1, 0); */
  /*   for(i=0; i<n1; i++) */
  /*   { */
  /*     if(i == Spad[m]) */
  /*     { */
  /* 	mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* 	continue; */
  /*     } */
  /*     mod_mul(sigchk, msvec[i], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   } */
  /*   mod_mul(smtst1, smtst1, omegat[m], crs.grp->tor); */
  /*   mod_add(tc1, tc1, smtst1, crs.grp->tor); */
    
  /*   /\* poly_elptic_mul(&Tmp, crs.ghat, smtst1, crs.Ex); *\/ */
  /*   /\* printf("m: %ld vl0[%ld] ", m, Spad[m]); *\/ */
  /*   /\* if(poly_point_cmp(Tmp, crs.vl0[Spad[m]])) *\/ */
  /*   /\*   printf("matches\n"); *\/ */
  /*   /\* else *\/ */
  /*   /\*   printf("fails\n"); *\/ */
       
  /*   mpz_powm_ui(cpow, c, Spad[m], crs.grp->tor); */
  /*   mpz_invert(cpow, cpow, crs.grp->tor); */
  /*   mpz_set_ui(smtst2, 0); */
  /*   for(i=0; i<k; i++) */
  /*   { */
  /*     if(i == Spad[m]) */
  /*     { */
  /* 	mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* 	continue; */
  /*     } */
  /*     mod_mul(sigchk, alpha[i], cpow, crs.grp->tor); */
  /*     mod_add(smtst2, smtst2, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*     /\* poly_elptic_mul(&Tmp, crs.ghat, sigchk, crs.Ex); *\/ */
  /*     /\* printf("i: %ld vprm[%ld] ", i, Spad[m]); *\/ */
  /*     /\* if(poly_point_cmp(Tmp, vprm[4*N*i + n1 - Spad[m]])) *\/ */
  /*     /\* 	printf("matches\n"); *\/ */
  /*     /\* else *\/ */
  /*     /\* 	printf("fails\n"); *\/ */
  /*   } */
  /*   mod_mul(smtst2, smtst2, omegat[m], crs.grp->tor); */
  /*   mod_add(tc2, tc2, smtst2, crs.grp->tor); */
    
  /*   mod_add(smtst1, smtst1, smtst2, crs.grp->tor); */
  /*   /\* poly_elptic_mul(&Tmp, crs.ghat, smtst1, crs.Ex); *\/ */
  /*   /\* printf("m: %ld vhat[%ld] ", m, Spad[m]); *\/ */
  /*   /\* if(poly_point_cmp(Tmp, vhat[Spad[m]])) *\/ */
  /*   /\*   printf("matches\n"); *\/ */
  /*   /\* else *\/ */
  /*   /\*   printf("fails\n"); *\/ */
  /*   /\* mod_mul(smtst1, smtst1, omegat[m], crs.grp->tor); *\/ */
  /*   mod_add(tc, tc, smtst1, crs.grp->tor); */
  /* } */
  /* poly_elptic_mul(&Tmp, crs.ghat, cpm, crs.Ex); */
  /* poly_point_printf("brute force sigma 2 part 2:\n", Tmp); */

// check sigagg 2 part 3

  /* mpz_set_ui(td, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_set_ui(smtst1, 0); */
  /*   i = X[1]; */
  /*   mpz_powm_ui(cpow, c, N + i - 1 - Spad[m], crs.grp->tor); */
  /*   xpnt = N + i - 1 - Spad[m]; */
  /*   while(i < X[2]) */
  /*   { */
  /*     mod_mul(sigchk, gma[i - 1], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*     xpnt++; */
  /*     i++; */
  /*   } */
  /*   mod_mul(smtst1, smtst1, omegat[m], crs.grp->tor); */
  /*   mod_add(td, td, smtst1, crs.grp->tor); */
  /* } */
  /* poly_elptic_mul(&Tmp, crs.ghat, smtst2, crs.Ex); */
  /* poly_point_printf("brute force sigma 2 part 3:\n", Tmp); */

// sig agg 3 check

  /* mpz_set_ui(ta, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_powm_ui(cpow, c, Spad[m]+1, crs.grp->tor); */
  /*   mod_div(sigchk, omegat[m], cpow, crs.grp->tor); */
  /*   mod_add(ta, ta, sigchk, crs.grp->tor); */
  /* } */
  /* point_init(&omc); */
  /* elptic_mul(&omc, crs.g, smtst1, crs.grp->E); */
  /* point_printf("computed sigagg3:\n", omc); */
    
// compare individual terms from every component of
// t*(sigma22+sigma23) - sigagg3*c4hat
// first term matches
  /* mpz_set_ui(cpm, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_powm_ui(sigchk, c, Spad[m], crs.grp->tor); */
  /*   mpz_invert(cpow, sigchk, crs.grp->tor); */
  /*   mpz_set_ui(smtst1, 0); */
  /*   for(i=0; i<n1; i++) */
  /*   { */
  /*     if(i == Spad[m]) */
  /*     { */
  /* 	mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* 	continue; */
  /*     } */
  /*     mod_mul(sigchk, msvec[Spad[m]], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   } */
  /*   mod_mul(smtst1, smtst1, omegat[m], crs.grp->tor); */
  /*   mod_add(cpm, cpm, smtst1, crs.grp->tor); */
  /* } */
  /* mod_mul(cpm, cpm, t, crs.grp->tor); */
  /* mpz_set_ui(smtst2, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_set(cpow, c); */
  /*   mpz_set_ui(smtst1, 0); */
  /*   for(i=0; i<n1; i++) */
  /*   { */
  /*     mod_mul(sigchk, msvec[Spad[m]], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   } */
  /*   mpz_powm_ui(cpow, c, Spad[m]+1, crs.grp->tor); */
  /*   mod_div(sigchk, smtst1, cpow, crs.grp->tor); */
  /*   /\* mpz_invert(cpow, sigchk, crs.grp->tor); *\/ */
  /*   /\* mod_mul(sigchk, cpow, smtst1, crs.grp->tor); *\/ */
  /*   mod_mul(smtst1, sigchk, omegat[m], crs.grp->tor); */
  /*   mod_add(smtst2, smtst2, smtst1, crs.grp->tor); */
  /* } */
  /* mod_mul(smtst2, smtst2, t, crs.grp->tor); */
  /* mod_sub(smtst1, cpm, smtst2, crs.grp->tor); */
  /* gmp_printf("first block difference: %Zd\n", smtst1); */
  /* mod_mul(smtst2, svec[0], t, crs.grp->tor); */
  /* mod_neg(smtst1, smtst2, crs.grp->tor); */
  /* gmp_printf("compares to: %Zd\n", smtst1); */

// second term matches
  /* mpz_set_ui(cpm, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_powm_ui(sigchk, c, Spad[m], crs.grp->tor); */
  /*   mpz_invert(cpow, sigchk, crs.grp->tor); */
  /*   mpz_set_ui(smtst1, 0); */
  /*   for(i=0; i<k; i++) */
  /*   { */
  /*     if(i == Spad[m]) */
  /*     { */
  /* 	mod_mul(cpow, cpow, c, crs.grp->tor); */
  /* 	continue; */
  /*     } */
  /*     mod_mul(sigchk, alpha[i], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   } */
  /*   mod_mul(smtst1, smtst1, omegat[m], crs.grp->tor); */
  /*   mod_add(cpm, cpm, smtst1, crs.grp->tor); */
  /* } */
  /* mod_mul(cpm, cpm, t, crs.grp->tor); */
  /* mpz_set_ui(smtst2, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_set(cpow, c); */
  /*   mpz_set_ui(smtst1, 0); */
  /*   for(i=0; i<k; i++) */
  /*   { */
  /*     mod_mul(sigchk, alpha[i], cpow, crs.grp->tor); */
  /*     mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*     mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   } */
  /*   mpz_powm_ui(cpow, c, Spad[m]+1, crs.grp->tor); */
  /*   mod_div(sigchk, smtst1, cpow, crs.grp->tor); */
  /*   mod_mul(smtst1, sigchk, omegat[m], crs.grp->tor); */
  /*   mod_add(smtst2, smtst2, smtst1, crs.grp->tor); */
  /* } */
  /* mod_mul(smtst2, smtst2, t, crs.grp->tor); */
  /* mod_sub(smtst1, cpm, smtst2, crs.grp->tor); */
  /* gmp_printf("first block difference: %Zd\n", smtst1); */
  /* mpz_set_ui(smtst1, 0); */
  /* for(i=0; i<T; i++) */
  /* { */
  /*   mod_mul(sigchk, omegat[i], alpha[Spad[i]], crs.grp->tor); */
  /*   mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /* } */
  /* mod_mul(smtst2, smtst1, t, crs.grp->tor); */
  /* mod_neg(smtst1, smtst2, crs.grp->tor); */
  /* gmp_printf("compares to: %Zd\n", smtst1); */

// third term checks
  /* mpz_set_ui(smtst2, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_powm_ui(cpow, c, N + 1 - Spad[m], crs.grp->tor); */
  /*   mod_mul(smtst1, gma[1], cpow, crs.grp->tor); */
  /*   mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   mod_mul(sigchk, gma[2], cpow, crs.grp->tor); */
  /*   mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*   mod_mul(smtst1, smtst1, omegat[m], crs.grp->tor); */
  /*   mod_add(smtst2, smtst2, smtst1, crs.grp->tor); */
  /* } */
  /* mod_mul(smtst2, smtst2, t, crs.grp->tor); */
  /* mpz_set_ui(cpm, 0); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mpz_powm_ui(cpow, c, N + 2, crs.grp->tor); */
  /*   mod_mul(smtst1, gma[1], cpow, crs.grp->tor); */
  /*   mod_mul(cpow, cpow, c, crs.grp->tor); */
  /*   mod_mul(sigchk, gma[2], cpow, crs.grp->tor); */
  /*   mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /*   mpz_powm_ui(cpow, c, Spad[m] + 1, crs.grp->tor); */
  /*   mod_div(smtst1, smtst1, cpow, crs.grp->tor); */
  /*   mod_mul(smtst1, smtst1, omegat[m], crs.grp->tor); */
  /*   mod_add(cpm, cpm, smtst1, crs.grp->tor); */
  /* } */
  /* mod_mul(cpm, cpm, t, crs.grp->tor); */
  /* mod_sub(smtst1, smtst2, cpm, crs.grp->tor); */
  /* gmp_printf("third term: %Zd\n", smtst1); */

// compare terms that should be divided -> subtracted in exponent

  /* mod_mul(sigchk, ta, tb1, crs.grp->tor); */
  /* mod_sub(smtst1, sigchk, tc1, crs.grp->tor); */
  /* gmp_printf("a b1 - c1: %Zd\n", smtst1); */
  /* gmp_printf("s_0: %Zd\n", svec[0]); */
  /* printf("\n"); */
  /* mod_mul(sigchk, ta, tb2, crs.grp->tor); */
  /* mod_sub(smtst1, sigchk, tc2, crs.grp->tor); */
  /* gmp_printf("a b2 - c2: %Zd\n", smtst1); */
  /* mpz_set_ui(smtst1, 0); */
  /* for(i=0; i<T; i++) */
  /* {            */
  /*   mod_mul(sigchk, omegat[i], alpha[Spad[i]], crs.grp->tor); */
  /*   mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /* } */
  /* gmp_printf("w0 a1 + w1 a2: %Zd\n", smtst1); */
  /* printf("\n"); */
  /* mod_mul(sigchk, ta, tb3, crs.grp->tor); */
  /* mod_sub(smtst1, sigchk, td, crs.grp->tor); */
  /* gmp_printf("a b3 - d: %Zd\n", smtst1); */
  
  /* exit(0); */
/********************************************************************************/

/* compute decryption components
   first is point sigma_agg_1. S list order, not user number.  */
#if 1
  printf("computing sigma agg 1\n");
  point_init(&sigagg1);
  point_init(&tmul);
  for(i=0; i<T; i++)
  {
    elptic_mul(&tmul, sigma1[i], omegat[i], crs.grp->E);
    elptic_sum(&sigagg1, sigagg1, tmul, crs.grp->E);
  }
  poly_init(&tstrght);
  poly_point_init(&G2);
  tog2(&G2, sigagg1);
  tate(&tstrght, G2, c3hat, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_printf("first bottom term (tate(sigagg1, c3hat)): \n", tstrght);
//#endif
  printf("computing sigma agg 3\n");
  point_init(&sigagg3);
  m = 2*N - 2;
  for(i=0; i<N; i++)
  {
    elptic_mul(&tmul, crs.cg[m - Spad[i]], omegat[i], crs.grp->E);
    elptic_sum(&sigagg3, sigagg3, tmul, crs.grp->E);
  }
//  point_printf("sigagg3:\n", sigagg3);
  
  printf("compute bottom second term\n");
  poly_init(&tatebot2);
  tog2(&G2, sigagg3);
  tate(&tatebot2, G2, c4hat, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_printf("tate(sig3, c4hat):\n", tatebot2);
#endif

// compute from exponent result
#if 0
  poly_init(&taugg);
  tog2(&G2, crs.g);
  tate(&taugg, G2, crs.ghat, crs.grp->S, crs.grp->tor, crs.Ex);
  mod_mul(sigchk, ta, tb, crs.grp->tor);
  poly_init(&tatebot2);
  poly_pow(&tatebot2, taugg, sigchk);
//  poly_printf("tate(g, ghat)^sig3*c4hat:\n", tatebot2);

#endif

  poly_init(&tatebot1);
  poly_mul(&tatebot1, tatebot2, tstrght);
  poly_printf("bottom term:\n", tatebot1);
#if 1
  printf("computing sigma agg 2 part 1\n");
  poly_point_init(&sigagg21);
  for(i=0; i<T; i++)
  {
    poly_elptic_mul(&Tmp, sigma2[i], omegat[i], crs.Ex);
    poly_elptic_sum(&sigagg21, sigagg21, Tmp, crs.Ex);
  }
#endif
#if 0
  poly_init(&tatetop);
  tog2(&G2, c2);
  tate(&tatetop, G2, sigagg21, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_init(&rcrvdmsg);
  poly_div(&rcrvdmsg, tatetop, tstrght);
  poly_printf("first term top / first term bottom:\n", rcrvdmsg);
  mpz_inits(smtst1, smtst2, cpow, sigchk, cpm, NULL);
  mpz_set_ui(smtst1, 0);
  for(i=0; i<T; i++)
  {           
    mod_mul(sigchk, omegat[i], alpha[Spad[i]], crs.grp->tor);
    mod_add(smtst1, smtst1, sigchk, crs.grp->tor);
  }
  poly_elptic_mul(&Tmp, crs.ghat, smtst1, crs.Ex);
  poly_init(&tatebot2);
  tate(&tatebot2, G2, Tmp, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_printf("what left over should be:\n", tatebot2);
  exit(0);
#endif       
#if 1
  printf("computing sigma agg 2 part 2\n");
  poly_point_init(&sigagg22);
  for(i=0; i<N; i++)
  {
    poly_elptic_mul(&Tmp, vhat[Spad[i]], omegat[i], crs.Ex);
    poly_elptic_sum(&sigagg22, sigagg22, Tmp, crs.Ex);
  }
  /* poly_point_printf("sigma agg 2 part 2 from tables:\n", sigagg22); */
  /* exit(0); */
//#endif  
  printf("computing sigma agg 2 part 3\n");
  poly_point_init(&sigagg23);
  mask = 1;
  n4 = 4*N - 1;
  m = 2*N;
  for(j=0; j<n; j++)
  {
    if(!(mask & lmt))      // only clear bits allowed in sum
    {
      printf("using j %ld\n", j);
      for(i=0; i<N; i++)
      {
	poly_elptic_mul(&Tmp, crs.tauhat[j*n4 + m + Spad[i]], omegat[i], crs.Ex);
	poly_elptic_sum(&sigagg23, sigagg23, Tmp, crs.Ex);
      }
    }
    mask <<= 1;
  }
//  poly_point_printf("sigma agg 2 part 3 table:\n", sigagg23);
  
  printf("compute top test tate(c2, sig21+sig22+sig23)\n");
  poly_init(&tstlft);
  poly_point_init(&sigagg2);
  poly_elptic_sum(&sigagg2, sigagg23, sigagg22, crs.Ex);
  poly_elptic_sum(&sigagg2, sigagg2, sigagg21, crs.Ex);
  tog2(&G2, c2);
  tate(&tstlft, G2, sigagg2, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_printf("tate(c2, sig21+sig22+sig23):\n", tstlft);

  poly_div(&tstrght, tstlft, tatebot1);
  poly_init(&rcrvdmsg);
  poly_mul(&rcrvdmsg, tstrght, C1);
  poly_printf("recovered message:\n", rcrvdmsg);
  exit(0);
#endif
// compute from exponent result
#if 0
  poly_init(&taugg);
  tog2(&G2, crs.g);
  tate(&taugg, G2, crs.ghat, crs.grp->S, crs.grp->tor, crs.Ex);
  mod_add(sigchk, tc, td, crs.grp->tor);
  mod_mul(sigchk, sigchk, t, crs.grp->tor);
  poly_init(&tatetop);
  poly_pow(&tatetop, taugg, sigchk);
//  poly_printf("tate(g, ghat)^t(sig22+sig23):\n", tatetop);
  poly_init(&tstlft);
  poly_div(&tstlft, tatetop, tatebot2);
  poly_printf("top / bottom:\n", tstlft);
#endif
#if 0
  mod_mul(smtst1, ta, tb, crs.grp->tor);
  mod_add(sigchk, tc, td, crs.grp->tor);
  mod_mul(sigchk, sigchk, t, crs.grp->tor);
  mod_sub(smtst2, sigchk, smtst1, crs.grp->tor);
  gmp_printf("exponent top - bottom: %Zd\n", smtst2);
#endif
#if 0
  poly_pow(&tstlft, taugg, smtst2);
  poly_printf("tate(g, ghat)^t(C+D)-AB:\n", tstlft);
  exit(0);
#endif
/* check reconstruction vector */

  /* mpz_inits(smtst1, smtst2, cpow, sigchk, cpm, NULL); */
  /* for(m=0; m<N; m++) */
  /* { */
  /*   mod_mul(sigchk, omegat[m], msvec[Spad[m]], crs.grp->tor); */
  /*   mod_add(smtst1, smtst1, sigchk, crs.grp->tor); */
  /* } */
  /* gmp_printf("msvec*omega: %Zd\n", smtst1); */
  /* gmp_printf("svec[0]: %Zd\n", svec[0]); */
  /* exit(0); */
// repeat test from above, but add in s_0

  mpz_set_ui(smtst1, 0);
  for(i=0; i<T; i++)
  {           
    mod_mul(sigchk, omegat[i], alpha[Spad[i]], crs.grp->tor);
    mod_add(smtst1, smtst1, sigchk, crs.grp->tor);
  }

//  mod_div(smtst1, smtst1, c, crs.grp->tor);

  mod_div(sigchk, svec[0], c, crs.grp->tor);
//  mod_add(smtst1, smtst1, svec[0], crs.grp->tor);
  mod_add(smtst1, smtst1, sigchk, crs.grp->tor);
  
  mod_mul(smtst1, smtst1, t, crs.grp->tor);  // t(s_0 + omga_0*a_1 + omga_1*a_2)
  mod_neg(sigchk, smtst1, crs.grp->tor);
  gmp_printf("left over terms: %Zd\n", sigchk);
  exit(0);
  
  poly_init(&tatetop);
  tog2(&G2, crs.g);
  tate(&tatetop, G2, crs.ghat, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_pow(&tatetop, tatetop, smtst1);
  poly_invert(&tatetop, tatetop);
//  poly_invert(&tatebot2, tatetop);
  
  poly_printf("what left over should be:\n", tatetop);

  exit(0);

  printf("computing sigma agg 2\n");
  poly_point_init(&sigagg2);
  poly_elptic_sum(&sigagg2, sigagg21, sigagg22, crs.Ex);
  poly_elptic_sum(&sigagg2, sigagg2, sigagg23, crs.Ex);

/* compute numerator of output */

  printf("computing numerator\n");
  tog2(&G2, c2);
  tate(&tatetop, G2, sigagg2, crs.grp->S, crs.grp->tor, crs.Ex);
  poly_div(&tatebot2, tatetop, tatebot1);
  poly_printf("top / bottom: \n", tatebot2);
  /* poly_init(&rcrvdmsg); */
  /* poly_mul(&rcrvdmsg, tatebot2, C1); */
  /* poly_printf("recovered message:\n", rcrvdmsg); */

  poly_pow(&tatetop, crs.B, t);
  poly_printf("B^t value:\n", tatetop);
  poly_invert(&tatebot1, tatetop);
  poly_printf("B^-t value:\n", tatebot1);
  /* exit(0); */
  /* poly_mul(&tatetop, tatetop, C1); */

/* and denominator */

  /* printf("computing denominator\n"); */
  /* poly_init(&tatebot1); */
  /* tog2(&G2, sigagg1); */
  /* tate(&tatebot1, G2, c3hat, crs.grp->S, crs.grp->tor, crs.Ex); */
  /* poly_init(&tatebot2); */
  /* tog2(&G2, sigagg3); */
  /* tate(&tatebot2, G2, c4hat, crs.grp->S, crs.grp->tor, crs.Ex); */
  
/* recover message data */

  /* poly_copy(&rcrvdmsg, tatetop); */
  /* poly_div(&rcrvdmsg, rcrvdmsg, tatebot1); */
  /* poly_div(&rcrvdmsg, rcrvdmsg, tatebot2); */
  /* poly_printf("recovered message: \n", rcrvdmsg); */
  
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
  /* free(Srow); */
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
  poly_clear(&tatetop);
  poly_point_clear(&G2);
  poly_clear(&tatebot1);
  poly_clear(&tatebot2);
  poly_clear(&rcrvdmsg);
}  
  
