# Makedefs used
include Makedefs

# compile flags
CPPFLAGS = -I$(NETCDFINC)
FCFLAGS = -c -I$(NETCDFINC) -I$(PWD) $(FFLAGS)

# lib flags
LDFC = $(NETCDFLIB)

# program name
PROGRAM = toolsf
MOD= ./*.o
# source files and objects

SRCS = tools_fort.F

all:$(RCS)
	$(CXX) $(SRCS) tools_fort_.f
	$(F90) -c module_file.F -fPIC;
	$(F90) -c spln2d.F -fPIC;
	$(FC) $(FCFLAGS) $(MOD) -m $(PROGRAM) tools_fort_.f $(LDFC)
	mv *.so ../


clean:
	rm -f *.o *.mod *.f
	rm -f ../*.so

