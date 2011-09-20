
SOURCES = src/main.cpp src/amplitude.cpp src/datafile.cpp \
	src/solver_force.cpp src/solver_chebyshev.cpp src/chebyshev_amplitude.cpp \
	src/chebyshev.cpp src/hankel.cpp  \
	src/spectrum.cpp src/wave_function.cpp \
	src/virtual_photon.cpp \
	../amplitudelib/tools/tools.cpp ../amplitudelib/tools/interpolation.cpp
FTSOURCES = src/fourier/fourier.c
FTOBJECTS = $(FTSOURCES:.c=.o)

FFTSOURCES = src/fft4g.c
LIBBCISOURCES = #libbci-1.1.0/bci.c libbci-1.1.0/tdspl.c libbci-1.1.0/tools.c
OBJECTS=$(SOURCES:.cpp=.o)
FFTOBJECTS=$(FFTSOURCES:.c=.o)


