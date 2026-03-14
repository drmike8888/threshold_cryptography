/***************************************************************
 *                                                             *
 *   Simulate partial decryption step. Reality requires way    *
 *   more complexity as each person does their own decryption  *
 *   and the decryption signatures are combined later. For     *
 *   simulation the private keys are in one file and the       *
 *   public keys are in another file so the algorithm can be   *
 *   examined. Private keys should only be in volitile memory  *
 *   and never saved.                                          *
 *                                                             *
 *                     Author = Mike Rosing                    *
 *                      Date = 25 Feb. 2026                    *
 *                                                             *
 **************************************************************/

#include "mpz_raw.h"
#include "eliptic.h"
#include "thrshldencrpt_setup.h"
#include "thrshldencrpt_crsin.h"

int main(int argc, char *argv[])
{
  FILE *key, *cphr, *slst;
  long i, N;
  long T, *S, slen;
  CRS crs;
  POLY *A;                      // vector
  POLY_POINT *vprm;             // matrix
  mpz_t tag, *alpha, *rl;
  POLY C1;
  POINT c2;
  POLY_POINT c3hat, c4hat, Tmp;
  POINT *sigma1;                 // vector
  POLY_POINT *sigma2;            // vector
  
/* input order is CRS, public key, private key, cipher text (binary!), and S list */

  if(argc < 6)
  {
    printf("Use: ./thrshldencrpt_prtldcrpt <crs> <public key> <private key> <cipher data> <S list>\n");
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

/* read in private key file
   THIS IS A SIMULATION - DON'T EVER DO THIS!
*/

  key = fopen(argv[3], "r");
  if(!key)
  {
    printf("can't find private key file %s\n", argv[3]);
    exit(-4);
  }
  fread(&i, sizeof(long), 1, key);
  if(i != N)
  {
    printf("inconsistent simulation data %ld vs %ld\n", i, N);
    exit(-5);
  }
  alpha = (mpz_t*)malloc(sizeof(mpz_t)*N);
  for(i=0; i<N; i++)
  {
    mpz_init(alpha[i]);
    mpz_inp_raw(alpha[i], key);
  }
  fclose(key);
  printf("private keys read in\n");
  
/* read in cipher text file */

  if(cipher_read(tag, &T, &C1, &c2, &c3hat, &c4hat, argv[4]) < 0)
    exit(-6);
  printf("cipher text file input\n");

/* read in S list of decryptors */

  slst = fopen(argv[5], "r");
  if(!slst)
  {
    printf("can't find decryptor list file %s\n", argv[5]);
    exit(-7);
  }
  fscanf(slst, "%ld", &slen);
  if(slen < T)
  {
    printf("number of decryptors %ld below threshold %ld\n", slen, T);
    exit(-8);
  }
  S = (long*)malloc(sizeof(long)*slen);
  for(i=0; i<slen; i++)
  {
    fscanf(slst, "%ld", &S[i]);
    S[i]--;                       // convert from user number to index
  }
  fclose(slst);
  printf("decryptor list read in\n");

/* create signature values for each person in S list */

  sigma1 = (POINT*)malloc(sizeof(POINT)*slen);
  sigma2 = (POLY_POINT*)malloc(sizeof(POLY_POINT)*slen);
  rl = (mpz_t*)malloc(sizeof(mpz_t)*slen);
  poly_point_init(&Tmp);
  key = fopen("debug.rl", "w");
  for(i=0; i<slen; i++)
  {
    printf("computing user %ld signatures\n", S[i]+1);
    mpz_init(rl[i]);
    mod_rand(rl[i], crs.grp->tor);
    mpz_out_raw(key, rl[i]);
    point_init(&sigma1[i]);
    elptic_mul(&sigma1[i], crs.g, rl[i], crs.grp->E);
    poly_point_init(&sigma2[i]);
    poly_elptic_mul(&Tmp, crs.ghat, alpha[S[i]], crs.Ex);
    poly_elptic_mul(&sigma2[i], crs.uhat, tag, crs.Ex);
    poly_elptic_sum(&sigma2[i], sigma2[i], crs.hhat, crs.Ex);
    poly_elptic_mul(&sigma2[i], sigma2[i], rl[i], crs.Ex);
    poly_elptic_sum(&sigma2[i], sigma2[i], Tmp, crs.Ex);
  }
  fclose(key);
  
/* save decryption values to disk */

  key = fopen("partial.sig", "w");
  fwrite(&slen, sizeof(long), 1, key);
  fwrite(S, sizeof(long), slen, key);
  for(i=0; i<slen; i++)
  {
    point_write(&sigma1[i], key);
    poly_point_write(&sigma2[i], key);
  }
  fclose(key);

/* clean up ram */

  crs_clear(&crs);
  key_clear(&A, &vprm, N);
  cipher_clear(&C1, &c2, &c3hat, &c4hat);
  for(i=0; i<N; i++)
    mpz_clear(alpha[i]);
  free(alpha);
  free(S);
  for(i=0; i<slen; i++)
  {
    point_clear(&sigma1[i]);
    poly_point_clear(&sigma2[i]);
    mpz_clear(rl[i]);
  }
  free(sigma1);
  free(sigma2);
  free(rl);
}

