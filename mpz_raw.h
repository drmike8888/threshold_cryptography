/*********************************************************************
 *                                                                   *
 *   I/O functions for points and curves on base and polynomial      *
 *   basis in binary format.                                         *
 *                                                                   *
 ********************************************************************/

#include "poly_eliptic.h"
#include "eliptic.h"

void point_write(POINT *P, FILE *f);
void point_read(POINT *P, FILE *f);
void curve_write(CURVE *E, FILE *f);
void curve_read(CURVE *E, FILE *f);
void poly_write(POLY *p, FILE *f);
void poly_read(POLY *p, FILE *f);
void poly_point_write(POLY_POINT *P, FILE *f);
void poly_point_read(POLY_POINT *P, FILE *f);
void poly_curve_write(POLY_CURVE *E, FILE *f);
void poly_curve_read(POLY_CURVE *E, FILE *f);
