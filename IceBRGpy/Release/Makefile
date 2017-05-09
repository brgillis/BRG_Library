CC = gcc # C compiler
CFLAGS = -fPIC -Wall -Wextra -std=c++11 -fopenmp -O3 -DNDEBUG -I../ -I/home/brg/include -I../../IceBRG_main -I../../IceBRG_physics -I../../IceBRG_lensing -I/usr/include/python2.7/ -I/home/brg/.local/lib/python2.7/site-packages/numpy/core/include # C flags
LDFLAGS =  -fopenmp -L/home/brg/lib -L../../IceBRG_main/Release -L../../IceBRG_physics/Release -L../../IceBRG_lensing/Release -lIceBRG_lensing -lIceBRG_physics -lIceBRG_main -lpython2.7 -lboost_serialization -lboost_log # linking flags
RM = rm -f  # rm command
TARGET_LIB = _cIceBRGpy.so # target lib
OBJECT = ../pysrc/cIceBRGpy_wrap.o


all: swig $(OBJECT) $(TARGET_LIB) 

default: $(TARGET_LIB)

$(OBJECT):	
	g++ -c -o $(OBJECT) $(CFLAGS) ../pysrc/cIceBRGpy_wrap.cxx

$(TARGET_LIB): 
	g++ -shared -o $(TARGET_LIB) -Wl,--whole-archive  $(LDFLAGS) -Wl,--no-whole-archive  ../pysrc/cIceBRGpy_wrap.o

swig:
	swig -c++ -python -outdir ../Release/ ../pysrc/cIceBRGpy.i

object: $(OBJECT)

lib: $(TARGET_LIB)

swig-clean:
	rm -f ../pysrc/cIceBRGpy_wrap.cxx ../Release/cIceBRGpy.py

object-clean:
	rm -f ../pysrc/cIceBRGpy_wrap.o

lib-clean:
	rm -f $(TARGET_LIB)

clean: swig-clean object-clean lib-clean
	rm -f Makefile.bak

# g++ -std=c++0x -DNDEBUG -DEIGEN_DONT_PARALLELIZE -I"/disk2/brg/git/brg_library/IceBRG_lensing" -I"/disk2/brg/git/brg_library/IceBRG_physics" -I"/disk2/brg/git/brg_library/IceBRG_main" -I/home/brg/include -I/usr/include/python2.7 -O3 -fopenmp -funroll-loops -march=native -mfpmath=sse -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"pysrc/cIceBRGpy_wrap.d" -MT"pysrc/cIceBRGpy_wrap.d" -o "pysrc/cIceBRGpy_wrap.o" "../pysrc/cIceBRGpy_wrap.cxx"


# g++ -L"/disk2/brg/git/brg_library/IceBRG_lensing/Release" -L"/disk2/brg/git/brg_library/IceBRG_physics/Release" -L"/disk2/brg/git/brg_library/IceBRG_main/Release" -L/home/brg/lib -fopenmp -shared -o "_cIceBRGpy.so"  ./pysrc/cIceBRGpy_wrap.o   -lIceBRG_lensing -lIceBRG_physics -lIceBRG_main -lpython2.7 -lboost_serialization -lboost_log



