modulo.o:	modulo.c modulo.h
	gcc -c -o modulo.o modulo.c

poly.o:  	poly.h poly.c
	gcc -c -o poly.o poly.c

eliptic.o:	eliptic.c eliptic.h
	gcc -c -o eliptic.o eliptic.c

poly_eliptic.o: poly_eliptic.h poly_eliptic.c
	gcc -c -o poly_eliptic.o poly_eliptic.c

test_mod:	test_mod.c modulo.o poly.o
	gcc -o test_mod test_mod.c modulo.o poly.o -lgmp

test_irrd:	test_irrd.c modulo.o poly.o eliptic.o poly_eliptic.o
	gcc -o test_irrd test_irrd.c modulo.o poly.o eliptic.o poly_eliptic.o -lgmp

pairing_gen:	pairing_gen.c modulo.o poly.o
	gcc -o pairing_gen pairing_gen.c modulo.o poly.o -lgmp -lm

pairing_sweep:	pairing_sweep.c 
	gcc -o pairing_sweep pairing_sweep.c -lgmp -lm

pairing_phi6: pairing_phi6.c
	gcc -o pairing_phi6 pairing_phi6.c -lgmp

pairing_sweep_alpha:	pairing_sweep_alpha.c 
	gcc -o pairing_sweep_alpha pairing_sweep_alpha.c -lgmp -lm

get_curve:	get_curve.c modulo.o poly.o eliptic.o
	gcc -o get_curve get_curve.c modulo.o poly.o eliptic.o -lgmp

eliptic_test:	eliptic_test.c modulo.o eliptic.o
	gcc -o eliptic_test eliptic_test.c eliptic.o modulo.o -lgmp

pairing.o:	pairing.c pairing.h
	gcc -c -o pairing.o pairing.c

weil_6_bit_pairing:	weil_6_bit_pairing.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc -o weil_6_bit_pairing weil_6_bit_pairing.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o -lgmp

tate_6_bit_pairing:	tate_6_bit_pairing.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc -o tate_6_bit_pairing tate_6_bit_pairing.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o -lgmp

quotient_group:	quotient_group.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc -o quotient_group quotient_group.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o -lgmp

signature.o: signature.c signature.h
	gcc -c -o signature.o signature.c

signatures_11: signatures_11.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o signatures_11 signatures_11.c  modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

signatures_11_keygen: signatures_11_keygen.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o signatures_11_keygen signatures_11_keygen.c  modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

base_curves_embed:  base_curves_embed.c modulo.o
	gcc -o base_curves_embed base_curves_embed.c modulo.o -lgmp

base_protocols.o: base_protocols.c base_protocols.h
	gcc -c -o base_protocols.o base_protocols.c

base_test:  base_test.c base_protocols.o modulo.o eliptic.o 
	gcc -o base_test base_test.c base_protocols.o eliptic.o modulo.o -lgmp -lk12

base_curve_gen: base_curve_gen.c modulo.o eliptic.o
	gcc -o base_curve_gen base_curve_gen.c eliptic.o modulo.o -lgmp

thrshldencrpt_setup: thrshldencrpt_setup.c thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o mpz_raw.o
	gcc -o thrshldencrpt_setup thrshldencrpt_setup.c thrshldencrpt_crs.o modulo.o poly.o eliptic.o poly_eliptic.o pairing.o mpz_raw.o -lgmp

thrshldencrpt_crsin.o: thrshldencrpt_crsin.c
	gcc -c -o thrshldencrpt_crsin.o thrshldencrpt_crsin.c

thrshldencrpt_keygen: thrshldencrpt_keygen.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o mpz_raw.o
	gcc -o thrshldencrpt_keygen thrshldencrpt_keygen.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o poly.o eliptic.o poly_eliptic.o pairing.o mpz_raw.o -lgmp -lk12

thrshldencrpt_preprss: thrshldencrpt_preprss.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o mpz_raw.o
	gcc -o thrshldencrpt_preprss thrshldencrpt_preprss.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o poly.o eliptic.o poly_eliptic.o mpz_raw.o -lgmp

mod_matinv.o:  mod_matinv.c
	gcc -c -o mod_matinv.o mod_matinv.c
thrshldencrpt_encrypt: thrshldencrpt_encrypt.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o mpz_raw.o
	gcc -o thrshldencrpt_encrypt thrshldencrpt_encrypt.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o poly.o eliptic.o poly_eliptic.o mpz_raw.o -lgmp

thrshldencrpt_prtldcrpt: thrshldencrpt_prtldcrpt.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o mpz_raw.o
	gcc -o thrshldencrpt_prtldcrpt thrshldencrpt_prtldcrpt.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o poly.o eliptic.o poly_eliptic.o mpz_raw.o -lgmp

thrshldencrpt_prtlvrfy: thrshldencrpt_prtlvrfy.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o mpz_raw.o
	gcc -o thrshldencrpt_prtlvrfy thrshldencrpt_prtlvrfy.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o poly.o eliptic.o poly_eliptic.o pairing.o mpz_raw.o -lgmp

thrshldencrpt_dcrypt: thrshldencrpt_dcrypt.c  thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o mpz_raw.o mod_matinv.o
	gcc -o thrshldencrpt_dcrypt thrshldencrpt_dcrypt.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o mpz_raw.o mod_matinv.o -lgmp

thrshldencrpt_bruteforcecheck: thrshldencrpt_bruteforcecheck.c  thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o mpz_raw.o mod_matinv.o
	gcc -o thrshldencrpt_bruteforcecheck thrshldencrpt_bruteforcecheck.c thrshldencrpt_crsin.o thrshldencrpt_crs.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o mpz_raw.o mod_matinv.o -lgmp

invtest: 	invtest.c mod_matinv.o
	gcc -o invtest invtest.c mod_matinv.o modulo.o -lgmp

thrshldencrpt_crs.o: thrshldencrpt_crs.c
	gcc -c -o thrshldencrpt_crs.o thrshldencrpt_crs.c

mpz_raw.o:	mpz_raw.c
	gcc -c -o mpz_raw.o mpz_raw.c
