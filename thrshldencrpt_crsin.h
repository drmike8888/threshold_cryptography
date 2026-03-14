/**************************************************************
 *                                                            *
 *    Read in CRS from disk and then free space when done.    *
 *                                                            *
 *************************************************************/

int crs_read(CRS *crs, char *filename);
void crs_clear(CRS *crs);
void mathinit(unsigned long seed, CRS crs);
int key_read(POLY **A, POLY_POINT **vprm, char *filename);
void key_clear(POLY **A, POLY_POINT **vprm, long N);
int encrptkey_read(long *N, long *k, long **L, POLY_POINT *zhat, POLY_POINT **vhat, char *filename);
void encrptkey_clear(long N, long **L, POLY_POINT * zhat, POLY_POINT **vhat);
int cipher_read(mpz_t tag, long *T, POLY *C1, POINT *c2, POLY_POINT *c3hat,
		POLY_POINT *c4hat, char *filename);
void cipher_clear(POLY *C1, POINT *c2, POLY_POINT *c3hat,
		  POLY_POINT *c4hat);



