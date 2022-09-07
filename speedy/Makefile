
DEFINES= -DKISS_FFT
CC=gcc
LIBDIR=
INCDIR=../kiss_fft130
CFLAGS=-Wall -g -fPIC -pthread -I$(INCDIR) $(DEFINES)

all: libspeedy.a

libspeedy.a: speedy.o
	ar cqs libspeedy.a speedy.o

speedy.o:
	$(CC) $(CFLAGS) -c speedy.c

clean:
	rm -f speedy.o libspeedy.a
