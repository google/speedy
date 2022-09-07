
CC=gcc
CPLUSPLUS=g++
CFLAGS=-g -DKISS_FFT

DEFINES="DEFINES=-DSONIC_INTERNAL"

all: libsonic.a speedy/libspeedy.a sonic/libsonic.a speedy_wave

speedy_wave: speedy_wave.cc libsonic.a sonic/wave.o
	$(CPLUSPLUS) $(CFLAGS) speedy_wave.cc -o speedy_wave libsonic.a sonic/wave.o -lm

libsonic.a:	soniclib.o speedy/speedy.o sonic/sonic.o sonic/spectrogram.o kiss_fft130/kiss_fft.o
	ar cqs libsonic.a soniclib.o speedy/speedy.o sonic/sonic.o sonic/spectrogram.o kiss_fft130/kiss_fft.o

soniclib.o: soniclib.c
	$(CC) $(CFLAGS) -c soniclib.c

speedy/libspeedy.a:
	cd speedy; make INCDIR=../kiss_fft130

speedy/speedy.o:
	cd speedy; make INCDIR=../kiss_fft130 speedy.o

sonic/libsonic.a:
	cd sonic; make INCDIR=../kiss_fft130 LIBDIR=../kiss_fft130 $(DEFINES) libsonic.a

sonic/sonic.o:
	cd sonic; make INCDIR=../kiss_fft130 LIBDIR=../kiss_fft130 $(DEFINES) sonic.o

sonic/wave.o:
	cd sonic; make INCDIR=../kiss_fft130 LIBDIR=../kiss_fft130 $(DEFINES) wave.o

sonic/spectrogram.o:
	cd sonic; make INCDIR=../kiss_fft130 LIBDIR=../kiss_fft130 $(DEFINES) spectrogram.o

kiss_fft130: kiss_fft130/kiss_fft.a
	cd kiss_fft130; make kiss_fft.a 

clean:
	cd speedy; make clean
	cd sonic; make clean
	cd kiss_fft130; make clean
	rm -f soniclib.o libsonic.a speedy_wave

