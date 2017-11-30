CC=gcc
CFLAGS=-g -gdwarf-2

fq-jt : fq-jt.c
	echo "Making fq-jt..."
	$(CC) $(CFLAGS) fq-jt.c -lz -o fq-jt

tgm : test-gsl-mat.c
	echo "Making tgm..."
	$(CC) -I/usr/include -L/usr/lib/x86_64-linux-gnu test-gsl-mat.c -o tgm -lgsl -lgslcblas -lm
