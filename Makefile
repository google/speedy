# Simple makefile to build the speedy library.
# This depends upon Sonic, which is assumed to be in a parallel directory.

# To install fftw this might work on your system
#		sudo apt-get install -y fftw-dev

CC=gcc
CPLUSPLUS=g++
CFLAGS=-g -DFFTW -fPIC

all: libspeedy.so speedy_wave

speedy_wave: speedy_wave.cc libspeedy.so ../sonic/libsonic_internal.so 
	$(CPLUSPLUS) $(CFLAGS) speedy_wave.cc libspeedy.so ../sonic/libsonic_internal.so -lc -lfftw3 -o speedy_wave

libspeedy.so: soniclib.o speedy.o 
	$(CC) -shared soniclib.o speedy.o -o libspeedy.so

soniclib.o: sonic.h speedy.h

speedy.o: speedy.h

clean:
	rm *.o *.so

