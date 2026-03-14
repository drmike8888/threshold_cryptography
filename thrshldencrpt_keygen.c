/***********************************************************
 *                                                         *
 *   Use pass phases for each of the N users to create     *
 *   private and public keys (called secret key in paper)  *
 *   to create the list of keys for use with threshold     *
 *   encryption algorithm.  private keys should never be   *
 *   stored - users should input pass phrase to create the *
 *   secret value so memory can be wiped clean.            *
 *                                                         *
 *                    Author = Mike Rosing                 *
 *                     Date = 19 Feb. 2026                 *
 *                                                         *
 **********************************************************/

#include "mpz_raw.h"
#include "eliptic.h"
#include "pairing.h"
#include "thrshldencrpt_setup.h"
#include "thrshldencrpt_crsin.h"
#include "k12.header/KangarooTwelve.h"

/*  add number of bits of security to bits in prime to find number of
    bytes required from extensible hash. Used in hash0 and hash1.
*/

long secbyte(mpz_t p)
{
  int m, k, b;
  
/*  find security level and add to bit size of prime  */
  
  m = mpz_sizeinbase(p, 2);
  if(m < 208)
    k = 80;
  else if(m < 320)
    k = 128;
  else if(m < 448)
    k = 192;
  else
    k = 256;
  b = m + k + 7;
  b >>= 3;
  return b;
}

/*  Hash to GF(r) using KangarooTwelve.
    Enter with pointer to data as bytes, number of bytes and torsion
    value r.
    Returns value mod r.
    NOTE: we assume p is twice the protection level.  Protection
    levels are 80, 128, 192 and 256 bits, so we add 10, 16, 24 or
    32 bytes to ensure result mod p is actually random.  See
    IETF irtf-cfrg-hash-to-curve-11 for details section 5.
*/

void hash1(mpz_t hsh, unsigned char *dat, long len, mpz_t r)
{
  unsigned char *outp, *dst;
  mpz_t zero;
  long b;
  
  dst = (char*)malloc(24);
  sprintf(dst, "Hash_1 pring&sig");

/*  find security level and add to bit size of prime  */

  b = secbyte(r);
  outp = (unsigned char*)malloc(b + 2);
  KangarooTwelve(dat, len, outp, b, dst, 16);

/*  convert hash bytes into an integer  */

  mpz_import(hsh, b, -1, 1, 0, 0, outp);
  mpz_init(zero);
  mod_add(hsh, hsh, zero, r);      // force to be mod r
  mpz_clear(zero);
  free(outp);
  free(dst);
}

int main(int argc, char *argv[])
{
  FILE *crsdat, *phrs, *keys;
  long i, j, k, n, N, n1, m, *phstrt;
  CRS crs;
  char *phrsbfr;
  mpz_t *alpha;         // vector
  POLY_POINT *vprm;     // matrix
  POLY *A;              // vector
  POINT Tmp;
  POLY_POINT Tmp2;
  
/* are there enough items on command line? */

  if(argc < 3)
  {
    printf("Use: ./threshold_encrypt_keygen <CRS file> <phrase file>\n");
    exit(-1);
  }

  phrs = fopen(argv[2], "r");
  if(!phrs)
  {
    printf("can't find phrase file %s\n", argv[2]);
    exit(-3);
  }

/* read in CRS data (or die) */

  if(crs_read(&crs, argv[1]) < 0)
    exit(-2);
  printf("data read in\n");
  n = crs.n;
  N = crs.N;
  mathinit(239553881, crs);

/* read in all phrases for every user */
  
  phrsbfr = (char*)malloc(1024*1024);
  i = 0;
  while(!feof(phrs))
  {
    phrsbfr[i] = fgetc(phrs);
    i++;
  }
  i--;
  fclose(phrs);

/* make sure there are enough phrases for all possible users.
   There are many possible ways to create and modify the public
   key data - most likely one person at a time, not in bulk as 
   in this example. */

  j = 0;
  k = 0;
  phstrt = (long*)malloc(sizeof(long)*(N + 1));
  while((j < i) && (k < N))
  {
    phstrt[k] = j;
    while(phrsbfr[j] != '\n') j++;
    k++;
  }
  phstrt[k] = j;
  printf("read in %ld phrases\n", k);
  if(k < N)
  {
    printf("not enough phrases in phrase file: %ld out of %ld\n", k, N);
    exit(-4);
  }

/* convert each phrase into an integer via hashing to get each user's private key */

  alpha = (mpz_t*)malloc(sizeof(mpz_t)*N);
  for(k=0; k<N; k++)
  {
    m = phstrt[k + 1] - phstrt[k];
    mpz_init(alpha[k]);
    hash1(alpha[k], &phrsbfr[phstrt[k]], m, crs.grp->tor);
  }

/* compute each users public key */

  A = (POLY*)malloc(sizeof(POLY)*N);
  point_init(&Tmp);
  poly_point_init(&Tmp2);
  for(k=0; k<N; k++)
  {
    printf("computing public key %ld\r", k);
    fflush(stdout);
    poly_init(&A[k]);
    elptic_mul(&Tmp, crs.g, alpha[k], crs.grp->E);
    tog2(&Tmp2, Tmp);
    weil(&A[k], Tmp2, crs.ghat, crs.grp->S, crs.grp->tor, crs.Ex);
  }
  printf("\n");

/* compute hint matrix - one row for each user */

  m = 4*N - 1;
  vprm = (POLY_POINT*)malloc(sizeof(POLY_POINT)*N*m);
  for(k=0; k<N; k++)
  {
    for(i=0; i<m; i++)
    {
      printf("computing hint %ld, %ld\r", k, i);
      fflush(stdout);
      poly_point_init(&vprm[k*m + i]);
      poly_elptic_mul(&vprm[k*m + i], crs.cghat[i], alpha[k], crs.Ex);
    }
  }
  printf("\n");

/* save public key block to disk. */

  keys = fopen("public_keys.bin", "w");
  fwrite(&N, sizeof(long), 1, keys);
  for(k=0; k<N; k++)
  {
    poly_write(&A[k], keys);
    for(i=0; i<m; i++)
      poly_point_write(&vprm[k*m + i], keys);
  }
  fclose(keys);

/* save secret keys to disk.
   Never actually do this - secret keys should always be
   created before each use and erased from memory carefully.
*/
  keys = fopen("secret_keys.bin", "w");
  fwrite(&N, sizeof(long), 1, keys);
  for(k=0; k<N; k++)
    mpz_out_raw(keys, alpha[k]);
  fclose(keys);

  
  crs_clear(&crs);
  free(phstrt);
  for(k=0; k<N; k++)
  {
    mpz_clear(alpha[k]);
    poly_clear(&A[k]);
  }
  free(alpha);
  free(A);
}
