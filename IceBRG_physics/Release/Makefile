CC = gcc # C compiler
CFLAGS = -fPIC -Wall -Wextra -std=c++11 -fopenmp -O3 -I../ -I/home/brg/include -I../../IceBRG_main # C flags
LDFLAGS =  -fopenmp # linking flags
RM = rm -f  # rm command
TARGET_LIB = libIceBRG_main.a # target lib

OBJECTS = ../src/lib/astro_caches.o  ../src/lib/astro.o  ../src/lib/correlation_function_estimator.o  ../src/lib/density_profile/density_profile.o  ../src/lib/density_profile/density_profile_functors.o  ../src/lib/density_profile/point_mass_profile.o  ../src/lib/density_profile/tNFW_profile.o  ../src/lib/lensing_correlation_function_estimator.o  ../src/lib/sky_obj/galaxy.o  ../src/lib/sky_obj/galaxy_group.o  ../src/lib/sky_obj/position_grid_cache.o  ../src/lib/sky_obj/sky_obj.o 

all: $(TARGET_LIB)

default: dep $(TARGET_LIB)

.cpp.o:
	g++ -c $(CFLAGS) -c $< -o $@

$(TARGET_LIB): $(OBJECTS)
	ar rcs $(TARGET_LIB) $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(TARGET_LIB) Makefile.bak
