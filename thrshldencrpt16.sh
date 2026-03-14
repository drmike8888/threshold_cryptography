#!/bin/bash
#          run all threshold encrypt routines to create files in order
./thrshldencrpt_setup threshold.curve
./thrshldencrpt_keygen threshold_setup.bin phrase_file.txt
./thrshldencrpt_preprss threshold_setup.bin public_keys.bin encrypt_quorum_16.rand
./thrshldencrpt_encrypt 5 threshold_setup.bin encryption_key.bin tag.file message.txt
./thrshldencrpt_prtldcrpt threshold_setup.bin public_keys.bin secret_keys.bin cipher_text.bin Slist16.txt
./thrshldencrpt_prtlvrfy threshold_setup.bin public_keys.bin cipher_text.bin  partial.sig
./thrshldencrpt_dcrypt threshold_setup.bin encryption_key.bin cipher_text.bin partial.sig

