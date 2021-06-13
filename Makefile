
# Makefile generated automatically by /Users/gsalam/scripts/mkcxx.pl '-f' '-1'
# run 'make make' to update it if you add new files

CXX = g++
CXXFLAGS = -Wall -g -O2

# also arrange for fortran support
FC = gfortran
FFLAGS = -Wall -O2
CXXFLAGS += -std=c++11
LDFLAGS += -std=c++11

FJCONFIG = fastjet-config
INCLUDE += `$(FJCONFIG) --cxxflags`
LIBRARIES  += `$(FJCONFIG) --libs --plugins`
INCLUDE += $(LCLINCLUDE)

COMMONSRC = 
F77SRC = 
COMMONOBJ = 

PROGSRC = example1.cc
PROGOBJ = example1.o

INCLUDE += 
LIBRARIES += 


all:  example1 


example1: example1.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)


make:
	/Users/gsalam/scripts/mkcxx.pl '-f' '-1'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  example1 

.cc.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.C.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.f.o:         $<
	$(FC) $(FFLAGS) -c $< -o $@
.f90.o:         $<
	$(FC) $(FFLAGS) -c $< -o $@


depend:
	makedepend  $(INCLUDE) -Y --   -- $(COMMONSRC) $(PROGSRC)
# DO NOT DELETE

example1.o: TwoBodyCuts.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/Selector.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/PseudoJet.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/config.h
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/config_auto.h
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/internal/numconsts.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/internal/base.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/internal/IsBase.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/SharedPtr.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/internal/deprecated.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/Error.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/PseudoJetStructureBase.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/RangeDefinition.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/LimitedWarning.hh
example1.o: /Users/gsalam/work/fastjet/fastjet-installation/include/fastjet/internal/thread_safety_helpers.hh
example1.o: Boson.hh RealRange.hh ArrVec.hh
