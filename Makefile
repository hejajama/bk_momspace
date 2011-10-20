CXXFLAGS = `gsl-config --cflags` -g -Wall -I ../lib/include -I ../amplitudelib -fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm -L ../lib/lib 
all: bk

include filelist.m

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
