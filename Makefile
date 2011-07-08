CXXFLAGS = `gsl-config --cflags` -g -Wall -I /linuxfs-home/student/hejajama/lib/include #-fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm -L /linuxfs-home/student/hejajama/lib/lib 

SOURCES = src/main.cpp src/amplitude.cpp src/tools.cpp src/datafile.cpp \
	src/solver_force.cpp src/solver_chebyshev.cpp src/chebyshev_amplitude.cpp \
	src/chebyshev.cpp src/hankel.cpp  \
	src/interpolation.cpp src/spectrum.cpp src/wave_function.cpp \
	src/virtual_photon.cpp
FTSOURCES = src/fourier/fourier.c
FTOBJECTS = $(FTSOURCES:.c=.o)

FFTSOURCES = src/fft4g.c
LIBBCISOURCES = #libbci-1.1.0/bci.c libbci-1.1.0/tdspl.c libbci-1.1.0/tools.c
OBJECTS=$(SOURCES:.cpp=.o)
FFTOBJECTS=$(FFTSOURCES:.c=.o)

all: bk

bk: $(OBJECTS) $(FFTOBJECTS) $(FTOBJECTS)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(FTOBJECTS) $(FFTOBJECTS) -o bk 
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(FFTOBJECTS)
	rm -f $(FTOBJECTS)
	rm -f $(OBJECTS)
	rm -f bk	
