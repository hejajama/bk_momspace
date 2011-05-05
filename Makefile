CXXFLAGS = `gsl-config --cflags` -g -Wall # -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm

SOURCES = src/main.cpp src/amplitude.cpp src/tools.cpp src/datafile.cpp
LIBBCISOURCES = #libbci-1.1.0/bci.c libbci-1.1.0/tdspl.c libbci-1.1.0/tools.c
OBJECTS=$(SOURCES:.cpp=.o)
LIBBCIOBJECTS=$(LIBBCISOURCES:.c=.o)

all: bk

bk: $(OBJECTS) $(LIBBCIOBJECTS)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(LIBBCIOBJECTS) -o bk 
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f bk	
