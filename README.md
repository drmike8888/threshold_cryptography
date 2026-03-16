# threshold_cryptography
Implementation of threshold encryption algorithm described in NIST presentation https://csrc.nist.gov/presentations/2026/mpts2026-1b3

This is a demonstration of mathematics. I leave it to professional
programmers to convert this into a useful application.

The PDF contains a description of all the math taken from the above
presentation. It also includes a description of all the code and why I
wrote it the way I did. I am certain there are better ways to
accomplish the same results. The main point is that this is REALLY COOL!
and I'd like to see it appear in a more useful manner.

The Makefile is a list of programs and is further proof I'm not a real
programmer. The files are compiled with the following commands:

make thrshldencrpt_setup

make thrshldencrpt_keygen

make thrshldencrpt_preprss

make thrshldencrpt_encrypt

make thrshldencrpt_prtldcrpt

make thrshldencrpt_prtlvrfy

make thrshldencrpt_dcrypt

The file thrshldencrpt16.sh contains all the commands and files in
the proper order to test each program. I cut and pasted from this
file to test things were working when I got tired of infinite loops
and seg faults when I forgot which order things were supposed to be.
Again - the math works, the code is really crude.

The program thrshldencrpt_bruteforcecheck.c uses math on the torsion
value to quickly check each equation from the elliptic curve
multiplication functions. The only reason to include it is to show how
debugging the math can be done "in the exponent". 