This code implements the Ring-TESLA signature scheme using parameters described in Improved Parameters for the Ring-TESLA Digital Signature Scheme by Arjun Chopra (eprint.iacr.org/2016/1099). 

Example uses of key generation, signing and verifying are given in test.c.

Compile with 'make all'. This will compile test.c for n=1024,827,709. Then run test_827, test_709, test_1024