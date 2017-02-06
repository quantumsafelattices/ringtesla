CC=gcc 
CFLAGS=-O3 -Wall -Wextra -std=gnu99 -pedantic -Wfatal-errors 

all: test_827 test_1024 test_709


test_827: test.c test_vectors_827.o tesla_827.o tesla_utils_827.o tesla_rand_827.o tesla_rand_openssl_aes_827.o FFT/FFT_1024_5767169.o FFT/FFT_827_33554699.o FFT/FFT_709_10486111.o
	$(CC) $(CFLAGS) -D TESLA_N=827  -o $@ $^ -lm -lcrypto

test_vectors_827.o: test_vectors.c
	$(CC) -c $(CFLAGS) -D TESLA_N=827  -o $@ $^ 

tesla_827.o: tesla.c
	$(CC) -c $(CFLAGS) -D TESLA_N=827  -o $@ $^ 

tesla_utils_827.o: tesla_utils.c
	$(CC) -c $(CFLAGS) -D TESLA_N=827  -o $@ $^ 

tesla_rand_827.o: tesla_rand.c
	$(CC) -c $(CFLAGS) -D TESLA_N=827  -o $@ $^ 

tesla_rand_openssl_aes_827.o: tesla_rand_openssl_aes.c
	$(CC) -c $(CFLAGS) -D TESLA_N=827  -o $@ $^ 



test_1024: test.c test_vectors_1024.o tesla_1024.o tesla_utils_1024.o tesla_rand_1024.o tesla_rand_openssl_aes_1024.o FFT/FFT_1024_5767169.o FFT/FFT_827_33554699.o FFT/FFT_709_10486111.o
	$(CC) $(CFLAGS) -D TESLA_N=1024  -o $@ $^ -lm -lcrypto

test_vectors_1024.o: test_vectors.c
	$(CC) -c $(CFLAGS) -D TESLA_N=1024  -o $@ $^ 

tesla_1024.o: tesla.c
	$(CC) -c $(CFLAGS) -D TESLA_N=1024  -o $@ $^ 

tesla_utils_1024.o: tesla_utils.c
	$(CC) -c $(CFLAGS) -D TESLA_N=1024  -o $@ $^ 

tesla_rand_1024.o: tesla_rand.c
	$(CC) -c $(CFLAGS) -D TESLA_N=1024  -o $@ $^ 

tesla_rand_openssl_aes_1024.o: tesla_rand_openssl_aes.c
	$(CC) -c $(CFLAGS) -D TESLA_N=1024  -o $@ $^ 



test_709: test.c test_vectors_709.o tesla_709.o tesla_utils_709.o tesla_rand_709.o tesla_rand_openssl_aes_709.o FFT/FFT_1024_5767169.o FFT/FFT_827_33554699.o FFT/FFT_709_10486111.o
	$(CC) $(CFLAGS) -D TESLA_N=709  -o $@ $^ -lm -lcrypto

test_vectors_709.o: test_vectors.c
	$(CC) -c $(CFLAGS) -D TESLA_N=709  -o $@ $^ 

tesla_709.o: tesla.c
	$(CC) -c $(CFLAGS) -D TESLA_N=709  -o $@ $^ 

tesla_utils_709.o: tesla_utils.c
	$(CC) -c $(CFLAGS) -D TESLA_N=709  -o $@ $^ 

tesla_rand_709.o: tesla_rand.c
	$(CC) -c $(CFLAGS) -D TESLA_N=709  -o $@ $^ 

tesla_rand_openssl_aes_709.o: tesla_rand_openssl_aes.c
	$(CC) -c $(CFLAGS) -D TESLA_N=709  -o $@ $^ 

clean:
	rm -f FFT/FFT_1024_5767169.o FFT/FFT_827_33554699.o FFT/FFT_709_10486111.o tesla_*.o test_827 test_1024 test_709 test test_vectors_*.o
