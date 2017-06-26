CC=gcc
CFLAGS=-g -gdwarf-2

fq-jt : fq-jt.c
	echo "Making fq-jt..."
	$(CC) $(CFLAGS) fq-jt.c -lz -o fq-jt
