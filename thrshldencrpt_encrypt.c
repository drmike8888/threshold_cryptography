/************************************************************
 *                                                          *
 *   Encryption step for threshold encryption algorithm.    *
 *   Input CRS data, encryption key, tag file, message file *
 *   and threshold for decryption. threshold must be less   *
 *   than or equal to L in the encryption key.              *
 *   For this excercise the tag file is converted to the    *
 *   torsion of the elliptic curve. If the message file is  *
 *   larger than the polynomial field size (rounded down to *
 *   bytes) it will be reduced by XOR with itself to fit.   *
 *                                                          *
 *                   Author = Mike Rosing                   *
 *                    Date = 23 Feb. 2026                   *
 *                                                          *
 ***********************************************************/

#include "mpz_raw.h"
#include "eliptic.h"
#include "thrshldencrpt_setup.h"
#include "thrshldencrpt_crsin.h"
#include <string.h>

int main(int argc, char *argv[])
{
  FILE *tagfl, *msg;
  long i, j, k, n, N, *L, m, n1, n2;
  long taglen, msglen, T, torlen, lmt;
  CRS crs;
  POLY *A;                      // vector
  POLY_POINT *vprm;             // matrix
  POLY_POINT zhat;
  POLY_POINT *vhat;             // vector
  char filename[256], *tgbf, *msgbf, *numbf;
  mpz_t tag, t;
  POLY msgpoly, C1, Inv;
  POINT c2;
  POLY_POINT c3hat, c4hat;
  long deg, prtl, mask;
  
/* get threshold number */

  if(argc >= 2)
  {
    T = atoi(argv[1]);
    if(T < 1)
    {
      printf("T must be > 1.\n");
      exit(-1);
    }
  }
  else
  {
    printf("Threshold for decryption: ");
    scanf("%ld", &T);
  }
  
/* read in CRS data (or die) */

  if(argc >= 3)
  {
    if(crs_read(&crs, argv[2]) < 0)
      exit(-2);
  }
  else
  {
    printf("CRS data file: ");
    scanf("%s", filename);
    if(crs_read(&crs, filename) < 0)
      exit(-2);
  }
  printf("CRS data read in\n");
  n = crs.n;
  N = crs.N;
  mathinit(261793148, crs);

/* read in encryption key */

  if(argc >= 4)
  {
    if(encrptkey_read(&N, &k, &L, &zhat, &vhat, argv[3]) < 0)
      exit(-3);
  }
  else
  {
    printf("encryption key file: ");
    scanf("%s", filename);
    if(encrptkey_read(&N, &k, &L, &zhat, &vhat, filename) < 0)
      exit(-3);
  }
  printf("encryption key read in\n");
  if(T > k)
  {
    printf("Threshold is larger than encryption key allows.\n");
    exit(-6);
  }
  
/* read in tag file */

  tgbf = (char*)malloc(1024);
  if(argc >= 5)
  {
    tagfl = fopen(argv[4], "r");
    if(!tagfl)
    {
      printf("can't find tag file %s\n", argv[4]);
      exit(-4);
    }
    i = 0;
    while((!feof(tagfl)) && (i < 1024))
    {
      tgbf[i] = fgetc(tagfl);
      i++;
    }
    fclose(tagfl);
    taglen = i - 1;
    printf("tag file read in\n");
  }
  else
  {
    while((getchar()) != '\n');
    printf("enter tag: ");
    tgbf[0] = 0;
    i = 0;
    while(tgbf[i] != '\n')
    {
      i++;
      tgbf[i] = getchar();
    }
    tgbf[i] = 0;
    if(strlen(tgbf) > 1022)
    {
      printf("buffer overflow!!\n");
      exit(-4);
    }
    taglen = strlen(tgbf);
  }
  printf("\n");

/* convert tag line into gmp number */

  mpz_init(tag);
  torlen = mpz_sizeinbase(crs.grp->tor, 16);
  torlen /= 2;
  if(taglen <= torlen)
    mpz_import(tag, taglen, 1, 1, 0, 0, (void*)tgbf);
  else
  {
    numbf = (char*)malloc(torlen + 1);
    for(i=0; i<torlen; i++)
      numbf[i] = 0;
    i = 0;
    while(i < taglen)
    {
      j = 0;
      while((j < torlen) && (i < taglen))
      {
	numbf[j] ^= tgbf[i];
	j++;
	i++;
      }
    }
    mpz_import(tag, torlen, 1, 1, 0, 0, (void*)numbf);
    free(numbf);
  }
  gmp_printf("tag: %Zd\n", tag); 
    
/* read in message file */

  if(argc >= 6)
  {
    msg = fopen(argv[5], "r");
    if(!msg)
    {
      printf("can't find message file %s\n", argv[5]);
      exit(-5);
    }
  }
  else
  {
    while((getchar()) != '\n');
    printf("enter message filename: ");
    scanf("%s", filename);
    msg = fopen(filename, "r");
    if(!msg)
    {
      printf("can't find message file %s\n", filename);
      exit(-5);
    }
  }

  msgbf = (char*)malloc(1024*1024);
  i = 0;
  while((!feof(msg)) && (i < 1024*1024))
  {
    msgbf[i] = fgetc(msg);
    i++;
  }
  fclose(msg);
  msglen = i - 1;
  printf("message buffer read in %ld\n", msglen);
  
/* convert message to polynomial value */

  poly_init(&msgpoly);
  torlen = mpz_sizeinbase(crs.grp->prm, 16);
  torlen /= 2;
  printf("torlen: %ld degree: %ld\n", torlen, crs.grp->degree);
  if(msglen < torlen*crs.grp->degree)
  {
    deg = msglen/torlen;
    prtl = msglen % torlen;
    if(prtl) deg++;
    for(i=0; i<deg; i++)
    {
      j = msglen - torlen*(i + 1);
      if(j < 0)
	mpz_import(msgpoly.coef[i], prtl, 1, 1, 0, 0, (void*)msgbf);
      else
	mpz_import(msgpoly.coef[i], torlen, 1, 1, 0, 0, (void*)&msgbf[j]);
    }
    msgpoly.deg = deg;
  }
  else
  {
    j = torlen*crs.grp->degree;
    numbf = (char*)malloc(j + 1);
    for(i=0; i<j; i++)
      numbf[i] = 0;
    i =  0;
    while(i < msglen)
    {
      j = 0;
      while(j < torlen*crs.grp->degree)
      {
	numbf[j] ^= msgbf[i];
	j++;
	i++;
      }
    }
    for(i=0; i<crs.grp->degree; i++)
    {
      j = torlen*(crs.grp->degree - i - 1);
      mpz_import(msgpoly.coef[i], torlen, 1, 1, 0, 0, (void*)&numbf[j]);
    }
    msgpoly.deg = crs.grp->degree;
    free(numbf);
  }
  poly_printf("message:\n", msgpoly);

/* compute cipher text terms */

  printf("computing C1\n");
  mpz_init(t);
  mod_rand(t, crs.grp->tor);
  poly_init(&C1);
  poly_pow(&C1, crs.B, t);
/***************************************************/
  /* tagfl = fopen("debug.t", "w"); */
  /* mpz_out_raw(tagfl, t); */
  /* fclose(tagfl); */
  
  /* poly_printf("B^t:\n", C1);    // DEBUGING! */
  /* poly_init(&Inv); */
  /* poly_invert(&Inv, C1); */
  /* poly_printf("1/B^t:\n", Inv); */
  /* poly_clear(&Inv); */
/***************************************************/
  
  poly_mul(&C1, C1, msgpoly);
  printf("computing c2\n");
  point_init(&c2);
  elptic_mul(&c2, crs.g, t, crs.grp->E);
  printf("computing c3 hat\n");
  poly_point_init(&c3hat);
  poly_elptic_mul(&c3hat, crs.uhat, tag, crs.Ex);
  poly_elptic_sum(&c3hat, c3hat, crs.hhat, crs.Ex);
  poly_elptic_mul(&c3hat, c3hat, t, crs.Ex);
  printf("computing c4 hat\n");
  poly_point_init(&c4hat);
  lmt = k - T;
  poly_point_copy(&c4hat, zhat);
  mask = 1;
  for(j=0; j<n; j++)
  {
    if(!(lmt & mask))
      poly_elptic_sum(&c4hat, c4hat, crs.yhat[j], crs.Ex);
    mask <<= 1;
  }
  poly_elptic_mul(&c4hat, c4hat, t, crs.Ex);
//  poly_point_printf("c4hat original:\n", c4hat);
/* save cipher text to disk */

  msg = fopen("cipher_text.bin", "w");
  mpz_out_raw(msg, tag);
  fwrite(&T, sizeof(long), 1, msg);
  poly_write(&C1, msg);
  point_write(&c2, msg);
  poly_point_write(&c3hat, msg);
  poly_point_write(&c4hat, msg);
  fclose(msg);
  
  encrptkey_clear(N, &L, &zhat, &vhat);
  crs_clear(&crs);
  mpz_clear(t);
  poly_clear(&C1);
  point_clear(&c2);
  poly_point_clear(&c3hat);
  poly_point_clear(&c4hat);
  free(tgbf);
  free(msgbf);
}
