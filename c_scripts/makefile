
# Flag Variables
CC=gcc
CXX=g++ 
CPPFLAGS=-g -Wall
RM=rm -f
# -Wall : this flag is used to turn on most compiler warnings
# -g : this flag adds debugging information to the executable file

# Define source file names
SPECSRCS=spectrum_make.cpp spectrum_funcs.cpp cosmology.cpp
SPECOBJS=$(subst .cpp,.o,$(SPECSRCS))

LCSRCS=light_curve_make.cpp light_curve_funcs.cpp spectrum_funcs.cpp cosmology.cpp
LCOBJS=$(subst .cpp,.o,$(LCSRCS))

RSPSRCS=response_make.cpp response_funcs.cpp spectrum_funcs.cpp
RSPOBJS=$(subst .cpp,.o,$(RSPSRCS))

TSTSRCS=testbench.cpp light_curve_funcs.cpp spectrum_funcs.cpp cosmology.cpp
TSTOBJS=$(subst .cpp,.o,$(TSTSRCS))


all: spectrum light_curve response

# Define executable file making methods
spectrum: $(SPECOBJS)
	$(CXX) $(CPPFLAGS) -o spectrum_make $(SPECOBJS)
	
light_curve: $(LCOBJS)
	$(CXX) $(CPPFLAGS) -o light_curve_make $(LCOBJS)

response: $(RSPOBJS)
	$(CXX) $(CPPFLAGS) -o response_make $(RSPOBJS)

testbench: $(TSTOBJS)
	$(CXX) $(CPPFLAGS) -o testbench $(TSTOBJS)

# Define dependency file making methods
cosmology.o: cosmology.cpp cosmology.hpp

make_spectrum.o: spectrum_make.cpp spectrum_funcs.hpp cosmology.hpp
spectrum_funcs.o: spectrum_funcs.cpp spectrum_funcs.hpp cosmology.hpp

make_light_curve.o: light_curve_make.cpp light_curve_funcs.hpp spectrum_funcs.hpp cosmology.hpp 
light_curve_funcs.o: light_curve_funcs.cpp light_curve_funcs.hpp spectrum_funcs.hpp cosmology.hpp 

make_response.o: response_make.cpp response_funcs.hpp spectrum_funcs.hpp 
response_funcs.o: response_funcs.cpp response_funcs.hpp spectrum_funcs.hpp

testbench.o: testbench.cpp light_curve_funcs.hpp spectrum_funcs.hpp cosmology.hpp 

# Define cleaning methods
clean:
	$(RM) $(SPECOBJS)
	$(RM) $(LCOBJS)
	$(RM) $(RSPOBJS)
	$(RM) $(TSTOBJS)
distclean: clean
	$(RM) spectrum_make
	$(RM) light_curve_make
	$(RM) response_make
	$(RM) testbench
testclean: 
	$(RM) $(TSTOBJS)
	$(RM) testbench

