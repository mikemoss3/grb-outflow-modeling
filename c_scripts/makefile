
# Flag Variables
CC=gcc
CXX=g++ 
CPPFLAGS=-g -Wall
RM=rm -f
# -Wall : this flag is used to turn on most compiler warnings
# -g : this flag adds debugging information to the executable file

# Define source file names
SPECSRCS=spectrum.cpp spectrum_funcs.cpp cosmology.cpp
SPECOBJS=$(subst .cpp,.o,$(SPECSRCS))

LCSRCS=light_curve.cpp light_curve_funcs.cpp spectrum_funcs.cpp cosmology.cpp
LCOBJS=$(subst .cpp,.o,$(LCSRCS))

TSTSRCS=testbench.cpp light_curve_funcs.cpp spectrum_funcs.cpp cosmology.cpp
TSTOBJS=$(subst .cpp,.o,$(TSTSRCS))


all: spectrum light_curve

# Define executable file making methods
spectrum: $(SPECOBJS)
	$(CXX) $(CPPFLAGS) -o spectrum $(SPECOBJS)
	
light_curve: $(LCOBJS)
	$(CXX) $(CPPFLAGS) -o light_curve $(LCOBJS)

testbench: $(TSTOBJS)
	$(CXX) $(CPPFLAGS) -o testbench $(TSTOBJS)

# Define dependency file making methods
testbench.o: testbench.cpp light_curve_funcs.hpp spectrum_funcs.hpp cosmology.hpp 

light_curve.o: light_curve.cpp light_curve_funcs.hpp spectrum_funcs.hpp cosmology.hpp 

spectrum.o: spectrum.cpp spectrum_funcs.hpp cosmology.hpp

light_curve_funcs.o: light_curve_funcs.cpp light_curve_funcs.hpp spectrum_funcs.hpp cosmology.hpp 

spectrum_funcs.o: spectrum_funcs.cpp spectrum_funcs.hpp cosmology.hpp

cosmology.o: cosmology.cpp cosmology.hpp

# Define cleaning methods
clean:
	$(RM) $(SPECOBJS)
	$(RM) $(LCOBJS)
	$(RM) $(TSTOBJS)
distclean: clean
	$(RM) spectrum
	$(RM) light_curve
	$(RM) testbench
testclean: 
	$(RM) $(TSTOBJS)
	$(RM) testbench

