
# Flag Variables
CC=gcc
CXX=g++ 
# CPPFLAGS=-g -Wall -I/opt/local/include -L/opt/local/lib/ -lcfitsio -I/opt/local/include/fitsio.h
CPPFLAGS=-g -Wall -I/opt/local/include 
CPPFLAGSAFT=-L/opt/local/lib/ -lcfitsio 
RM=rm -f
# -Wall : this flag is used to turn on most compiler warnings
# -g : this flag adds debugging information to the executable file

# Define source file names
SPECSRCS=spectrum_make.cpp spectrum_funcs.cpp Spectrum.cpp cosmology.cpp
SPECOBJS=$(subst .cpp,.o,$(SPECSRCS))

LCSRCS=light_curve_make.cpp light_curve_funcs.cpp spectrum_funcs.cpp Spectrum.cpp cosmology.cpp
LCOBJS=$(subst .cpp,.o,$(LCSRCS))

RSPSRCS=response_make.cpp response_funcs.cpp Response.cpp spectrum_funcs.cpp Spectrum.cpp cosmology.cpp
RSPOBJS=$(subst .cpp,.o,$(RSPSRCS))

TSTSRCS=testbench.cpp light_curve_funcs.cpp response_funcs.cpp Response.cpp spectrum_funcs.cpp Spectrum.cpp cosmology.cpp
# TSTSRCS=testbench.cpp
TSTOBJS=$(subst .cpp,.o,$(TSTSRCS))


all: spectrum light_curve response

# Define executable file making methods
spectrum: $(SPECOBJS)
	$(CXX) $(CPPFLAGS) -o spectrum_make $(SPECOBJS)
	
light_curve: $(LCOBJS)
	$(CXX) $(CPPFLAGS) -o light_curve_make $(LCOBJS)

response: $(RSPOBJS)
	$(CXX) $(CPPFLAGS) -o response_make $(RSPOBJS) $(CPPFLAGSAFT)

testbench: $(TSTOBJS)
	$(CXX) $(CPPFLAGS) -o testbench $(TSTOBJS) $(CPPFLAGSAFT)

# Define dependency file making methods
cosmology.o: cosmology.cpp cosmology.hpp

spectrum_make.o: spectrum_make.cpp spectrum_funcs.hpp Spectrum.hpp cosmology.hpp
spectrum_funcs.o: spectrum_funcs.cpp spectrum_funcs.hpp Spectrum.cpp cosmology.hpp
Spectrum.o: Spectrum.cpp Spectrum.hpp

light_curve_make.o: light_curve_make.cpp light_curve_funcs.hpp spectrum_funcs.hpp Spectrum.hpp cosmology.hpp 
light_curve_funcs.o: light_curve_funcs.cpp light_curve_funcs.hpp spectrum_funcs.hpp Spectrum.hpp cosmology.hpp 

response_make.o: response_make.cpp response_funcs.hpp Response.hpp spectrum_funcs.hpp Spectrum.hpp cosmology.cpp
response_funcs.o: response_funcs.cpp response_funcs.hpp Response.hpp spectrum_funcs.hpp Spectrum.hpp cosmology.cpp
Response.o: Response.cpp Response.hpp Spectrum.hpp cosmology.cpp

testbench.o: testbench.cpp light_curve_funcs.hpp response_funcs.hpp Response.hpp spectrum_funcs.hpp Spectrum.hpp cosmology.hpp
# testbench.o: testbench.cpp

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

