CC = gcc # C compiler
CFLAGS = -fPIC -Wall -Wextra -std=c++11 -fopenmp -O3 -I../ -I/home/brg/include -I../../IceBRG_main -I../../IceBRG_physics # C flags
LDFLAGS =  -fopenmp # linking flags
RM = rm -f  # rm command
TARGET_LIB = libIceBRG_main.a # target lib

OBJECTS = ../src/lib/lensing_tNFW_caches.o ../src/lib/lensing_tNFW_profile.o ../src/lib/lens_source_pair.o ../src/lib/magnification/expected_count_cache.o ../src/lib/magnification/expected_count_derivative_cache.o ../src/lib/magnification/expected_count_fit_loader.o ../src/lib/magnification/expected_count_loader.o ../src/lib/magnification/mag_calibration_cache.o ../src/lib/magnification/mag_calibration_loader.o ../src/lib/magnification/mag_correlation_function_estimator.o ../src/lib/magnification/magnification_alpha.o ../src/lib/magnification/magnification_functors.o ../src/lib/magnification/mag_signal_integral_cache.o ../src/lib/magnification/mag_weight_integral_cache.o ../src/lib/magnification/Schechter_like_functor.o ../src/lib/pair_bin.o ../src/lib/pair_binner.o ../src/lib/pair_bins_summary.o ../src/lib/pair_bin_summary.o ../src/lib/shifting/shifting_cache.o ../src/lib/shifting/shifting_loader.o


all: $(TARGET_LIB)

default: dep $(TARGET_LIB)

.cpp.o:
	g++ -c $(CFLAGS) -c $< -o $@

$(TARGET_LIB): $(OBJECTS)
	ar rcs $(TARGET_LIB) $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(TARGET_LIB) Makefile.bak
