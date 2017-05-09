CC = gcc # C compiler
CFLAGS = -fPIC -Wall -Wextra -std=c++11 -DNDEBUG -fopenmp -O3 -I../ -I/home/brg/include # C flags
LDFLAGS =  -fopenmp # linking flags
RM = rm -f  # rm command
TARGET_LIB = libIceBRG_main.a # target lib

OBJECTS = ../src/lib/globals.o ../src/lib/interpolator_derivative.o ../src/lib/interpolator.o ../src/lib/random_functions.o ../src/lib/sgsmooth.o ../src/lib/tk_spline.o

all: $(TARGET_LIB)

default: dep $(TARGET_LIB)

.cpp.o:
	g++ -c $(CFLAGS) -c $< -o $@

$(TARGET_LIB): $(OBJECTS)
	ar rcs $(TARGET_LIB) $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(TARGET_LIB) Makefile.bak
