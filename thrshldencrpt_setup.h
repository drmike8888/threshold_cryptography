/************************************************************************
 *                                                                      *
 *   Threshold specific functions as well as modular inverse routines.  *
 *                                                                      *
 ***********************************************************************/

/* create structure to hold main parameters of an
   elliptic curve group for threshold pairings. */

typedef struct
{
  long degree;
  mpz_t prm;   // coefficient prime
  mpz_t tor;   // elliptic curve torsion 
  mpz_t t;     // elliptic curve Frobenius
  mpz_t cardE; // base curve cardinality
  mpz_t cardEx;// extended curve cardinality
  mpz_t cobse; // base curve cofactor
  mpz_t coxtd; // extension curve cofactor
  CURVE E;     // base curve a4, a6 values
  POLY_POINT S;// reference point for Weil pairing
} GROUP;

/* this structure holds the common reference string values */

typedef struct
{
  long n;             // log_2(N)
  long N;             // max number of users (2^n)
  GROUP *grp;         // pointer
  POLY_CURVE Ex;
  POINT g;
  POLY_POINT ghat;
  POLY B;
  POLY_POINT uhat;
  POLY_POINT hhat;
  POLY_POINT *cghat;  // vector
  POINT      *cg;     // vector
  POLY_POINT z0hat;
  POLY_POINT *vl0;    // vector
  POLY_POINT *yhat;   // vector
  POLY_POINT *tauhat; // matrix
}CRS;

int mod_matinv(mpz_t *inv, long n, mpz_t *mat, mpz_t prime);
int mod_matmul(mpz_t *c, mpz_t *a, long m, long n,
	       mpz_t *b, long r, long s, mpz_t prime);
void genshare(mpz_t *share, long N, mpz_t mod);
void mdots(mpz_t ms, long row, mpz_t *s_vec, mpz_t *shrmtrx, long N, mpz_t mod);
void zvcalc(mpz_t *cvec, CRS *crs, mpz_t *shrmtrx, mpz_t c, mpz_t *s_vec, long N);
void ytaucalc(CRS *crs, mpz_t *cvec, long n, long N);
void cgcalc(CRS *crs, mpz_t *cvec, long N);
void tog2(POLY_POINT *G2, POINT G1);

