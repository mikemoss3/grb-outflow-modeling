
# Flag Variables
CC=gcc
CXX=g++ 
CPPFLAGS=-g -Wall --std=c++17
CPPFLAGSAFT=-lcfitsio # Link cfistsio library 
RM=rm -f
# -Wall : this flag is used to turn on most compiler warnings
# -g : this flag adds debugging information to the executable file

# Define source file names
MAINSRCS=main.cpp cosmology.cpp utilfuncs.cpp Spectrum.cpp LightCurve.cpp TTEs.cpp Response.cpp FitStats.cpp ShellDist.cpp ModelParams.cpp SynthGRB.cpp SynthGRBLibrary.cpp ObsGRB.cpp DataAnalysis.cpp
MAINOBS=$(subst .cpp,.o,$(MAINSRCS))

TSTSRCS=testbench.cpp cosmology.cpp utilfuncs.cpp Spectrum.cpp LightCurve.cpp TTEs.cpp Response.cpp FitStats.cpp ShellDist.cpp ModelParams.cpp SynthGRB.cpp SynthGRBLibrary.cpp ObsGRB.cpp DataAnalysis.cpp
# TSTSRCS=testbench.cpp ShellDist.cpp cosmology.cpp
TSTOBJS=$(subst .cpp,.o,$(TSTSRCS))

# Define executable file making methods
main: $(MAINOBS)
	$(CXX) $(CPPFLAGS) -o main $(MAINOBS) $(CPPFLAGSAFT)

testbench: $(TSTOBJS)
	$(CXX) $(CPPFLAGS) -o testbench $(TSTOBJS) $(CPPFLAGSAFT)

# Define dependency file making methods
cosmology.o: cosmology.cpp cosmology.hpp
ShellDist.o: ShellDist.cpp ShellDist.hpp utilfuncs.hpp cosmology.hpp
Spectrum.o: Spectrum.cpp Spectrum.hpp
LightCurve.o: LightCurve.cpp LightCurve.hpp
TTEs.o: TTEs.cpp TTEs.hpp LightCurve.hpp Spectrum.hpp
Response.o: Response.cpp Response.hpp

utilfuncs.o: utilfuncs.cpp utilfuncs.hpp Spectrum.hpp Response.hpp cosmology.hpp

ModelParams.o: ModelParams.cpp ModelParams.hpp
SynthGRB.o: SynthGRB.cpp SynthGRB.hpp ModelParams.hpp ShellDist.hpp LightCurve.hpp Spectrum.hpp utilfuncs.hpp cosmology.hpp
SynthGRBLibrary.o: SynthGRBLibrary.cpp SynthGRBLibrary.hpp SynthGRB.hpp ModelParams.hpp
ObsGRB.o: ObsGRB.cpp ObsGRB.hpp TTEs.hpp LightCurve.hpp Spectrum.hpp

FitStats.o: FitStats.cpp FitStats.hpp
DataAnalysis.o: DataAnalysis.cpp DataAnalysis.hpp FitStats.hpp SynthGRBLibrary.hpp SynthGRB.hpp ModelParams.hpp ObsGRB.hpp utilfuncs.hpp cosmology.hpp

main.o: main.cpp DataAnalysis.hpp FitStats.hpp ObsGRB.hpp SynthGRBLibrary.hpp SynthGRBLibrary.hpp Response.hpp TTEs.hpp LightCurve.hpp Spectrum.hpp utilfuncs.hpp cosmology.hpp

testbench.o: main.cpp DataAnalysis.hpp FitStats.hpp ObsGRB.hpp SynthGRBLibrary.hpp SynthGRBLibrary.hpp Response.hpp TTEs.hpp LightCurve.hpp Spectrum.hpp utilfuncs.hpp cosmology.hpp
# testbench.o: testbench.cpp ShellDist.hpp cosmology.hpp

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

