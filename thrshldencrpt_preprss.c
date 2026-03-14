/******************************************************************
 *                                                                *
 *   compute an encryption key for a specified quorum of users.   *
 *   Input CRS, list of public keys and list of indexes into the  *
 *   the public keys that make up the quorum of users doing the   *
 *   encryption.                                                  *
 *                                                                *
 *                     Author = Mike Rosing                       *
 *                      Date = 20 Feb. 2026                       *
 *                                                                *
 *****************************************************************/

#include "mpz_raw.h"
#include "eliptic.h"
#include "thrshldencrpt_setup.h"
#include "thrshldencrpt_crsin.h"

int main(int argc, char *argv[])
{
  FILE *qrm, *enct;
  long i, j, k, n, N, *L, m, n1, n2;
  CRS crs;
  POLY *A;                      // vector
  POLY_POINT *vprm;             // matrix
  POLY_POINT zhat;
  POLY_POINT *vhat;             // vector

  if(argc < 4)
  {
    printf("Use: ./thrshldencrpt_preprss <CRS data> <public key block> <index list of users>\n");
    exit(-1);
  }

/* read in CRS data (or die) */

  if(crs_read(&crs, argv[1]) < 0)
    exit(-2);
  printf("CRS data read in\n");
  n = crs.n;
  N = crs.N;
  mathinit(261793148, crs);

/* read in public keys (or die) */

  if(key_read(&A, &vprm, argv[2]) < 0)
    exit(-3);
  printf("public keys read in\n");

/* read in quorum list of encryptors */
  
  L = (long*)malloc(sizeof(long)*N);
  qrm = fopen(argv[3], "r");
  if(!qrm)
  {
    printf("can't find quorum list file %s\n", argv[3]);
    crs_clear(&crs);
    key_clear(&A, &vprm, N);
    exit(-4);
  }
  fscanf(qrm, "%ld", &k);
  if(k < 1)
  {
    printf("quorum too small? %ld\n", k);
    crs_clear(&crs);
    key_clear(&A, &vprm, N);
    exit(-5);
  }
  for(i=0; i<k; i++)
  {
    fscanf(qrm, "%ld", &L[i]);
    L[i]--;                        // convert user number to index
  }
  fclose(qrm);

/* compute z_hat point */

  printf("computing z hat\n");
  m = 4*N - 1;
  n2 = 2*N;
  poly_point_init(&zhat);
  for(i=0; i<k; i++)
    poly_elptic_sum(&zhat, zhat, vprm[L[i]*m + L[i] + n2], crs.Ex);
  poly_elptic_sum(&zhat, zhat, crs.z0hat, crs.Ex);

/* compute v_hat points. 2N - 1 components but only sum
   over L of them. */

  n1 = 2*N - 1;
  vhat = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n1);
  for(j=0; j<n1; j++)
  {
    printf("computing v hat %ld\r", j);
    fflush(stdout);
    poly_point_init(&vhat[j]);
    for(i=0; i<k; i++)
    {
      if(L[i] == j) continue;
      poly_elptic_sum(&vhat[j], vhat[j], vprm[L[i]*m + L[i] - j + n1], crs.Ex);
    }
    poly_elptic_sum(&vhat[j], vhat[j], crs.vl0[j], crs.Ex);
  }
  printf("\n");
  
/* save encryption keys to disk */

  enct = fopen("encryption_key.bin", "w");
  fwrite(&N, sizeof(long), 1, enct);
  fwrite(&k, sizeof(long), 1, enct);    // length of L
  fwrite(L, sizeof(long), k, enct);     // L table
  poly_point_write(&zhat, enct);
  for(j=0; j<n1; j++)
    poly_point_write(&vhat[j], enct);
  fclose(enct);

/* clear out memory */

  for(j=0; j<n1; j++)
    poly_point_clear(&vhat[j]);
  free(vhat);
  poly_point_clear(&zhat);
  crs_clear(&crs);
  key_clear(&A, &vprm, N);
}

