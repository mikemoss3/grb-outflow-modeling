
# Flag Variables
CC=gcc
CXX=g++ 
CPPFLAGS=-g -Wall -I/opt/local/include 
CPPFLAGSAFT=-L/opt/local/lib/ -lcfitsio # Link cfistsio library 
RM=rm -f
# -Wall : this flag is used to turn on most compiler warnings
# -g : this flag adds debugging information to the executable file

# Define source file names
MAINSRCS=main.cpp Response.cpp LightCurve.cpp Spectrum.cpp cosmology.cpp
MAINOBS=$(subst .cpp,.o,$(MAINSRCS))

TSTSRCS=testbench.cpp Response.cpp LightCurve.cpp Spectrum.cpp cosmology.cpp
# TSTSRCS=testbench.cpp
TSTOBJS=$(subst .cpp,.o,$(TSTSRCS))

# Define executable file making methods
main: $(MAINOBS)
	$(CXX) $(CPPFLAGS) -o main $(MAINOBS) $(CPPFLAGSAFT)

testbench: $(TSTOBJS)
	$(CXX) $(CPPFLAGS) -o testbench $(TSTOBJS) $(CPPFLAGSAFT)


# Define dependency file making methods
cosmology.o: cosmology.cpp cosmology.hpp

Spectrum.o: Spectrum.cpp Spectrum.hpp cosmology.hpp 

LightCurve.o: LightCurve.cpp LightCurve.hpp Spectrum.hpp cosmology.hpp

Response.o: Response.cpp Response.hpp Spectrum.hpp cosmology.cpp

main.o: Main.cpp Response.hpp LightCurve.hpp Spectrum.hpp cosmology.hpp

testbench.o: testbench.cpp Response.hpp LightCurve.hpp Spectrum.hpp cosmology.hpp
# testbench.o: testbench.cpp

# Define cleaning methods
clean:
	$(RM) $(MAINOBS)
	$(RM) $(TSTOBJS)
distclean: clean
	$(RM) main
	$(RM) testbench
testclean: 
	$(RM) $(TSTOBJS)
	$(RM) testbench

