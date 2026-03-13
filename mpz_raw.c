#include "mpz_raw.h"

/* Write and read points, polynomials and poly_points to
   a disk file. Enter with pointer to object, and FILE pointer.
   On reads, object must be pre-initialized.
*/

void point_write(POINT *P, FILE *f)
{
  mpz_out_raw(f, P->x);
  mpz_out_raw(f, P->y);
}

void point_read(POINT *P, FILE *f)
{
  mpz_inp_raw(P->x, f);
  mpz_inp_raw(P->y, f);
}

void curve_write(CURVE *E, FILE *f)
{
  mpz_out_raw(f, E->a4);
  mpz_out_raw(f, E->a6);
}

void curve_read(CURVE *E, FILE *f)
{
  mpz_inp_raw(E->a4, f);
  mpz_inp_raw(E->a6, f);
}

void poly_write(POLY *p, FILE *f)
{
  int i, k;

  k = p->deg;
  fwrite(&k, sizeof(long), 1, f);
  for(i=0; i<=k; i++)
    mpz_out_raw(f, p->coef[i]);
}

void poly_read(POLY *p, FILE *f)
{
  int i, k;

  fread(&k, sizeof(long), 1, f);
  p->deg = k;
  for(i=0; i<=k; i++)
    mpz_inp_raw(p->coef[i], f);
}

void poly_point_write(POLY_POINT *P, FILE *f)
{
  poly_write(&P->x, f);
  poly_write(&P->y, f);
}

void poly_point_read(POLY_POINT *P, FILE *f)
{
  poly_read(&P->x, f);
  poly_read(&P->y, f);
}

void poly_curve_write(POLY_CURVE *E, FILE *f)
{
  poly_write(&E->a4, f);
  poly_write(&E->a6, f);
}

void poly_curve_read(POLY_CURVE *E, FILE *f)
{
  poly_read(&E->a4, f);
  poly_read(&E->a6, f);
}
