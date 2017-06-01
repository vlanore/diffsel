#To compile with gsl, reactivate the following line (remove '#')
#USE_GSL=-DGSL

#To statically build the executables, uncomment the following line
#STATIC=-static

ifdef USE_GSL
#GSLLIB=-lgsl -lgslcblas
else
GSLLIB=
endif

LIBS=$(GSLLIB)
#LIBS=$(GSLLIB) -lboost_regex

CC=g++
SYSLIB=
INCLUDES=
#CPPFLAGS= $(USE_GSL) $(INCLUDES)
CPPFLAGS= -std=c++11 -g -Wall -O3 $(USE_GSL) $(INCLUDES)
LDFLAGS=$(STATIC)
INSTALL_DIR=
INSTALL_LIB=
SRCS= BranchSitePath.cpp Chrono.cpp CodonSequenceAlignment.cpp CodonStateSpace.cpp CodonSubMatrix.cpp GTRSubMatrix.cpp PhyloProcess.cpp Random.cpp SequenceAlignment.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp linalg.cpp linalg2.cpp Chain.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=..
ALL= globom siteom
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

# If we are on a windows platform, executables are .exe files
# [S.Bigaret] we assume here that we are compiling w/ MinGW, not Cygwin
EXEEXT=
uname := $(shell uname -s 2>/dev/null)
ifeq (MINGW,$(findstring MINGW,$(uname)))
	EXEEXT=.exe
	PROGS=$(addsuffix $(EXEEXT), $(ALL))
endif

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

globom$(EXEEXT): $(PROGSDIR)/globom$(EXEEXT)
$(PROGSDIR)/globom$(EXEEXT): SingleOmega.o $(OBJS)
	$(CC) SingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

siteom$(EXEEXT): $(PROGSDIR)/siteom$(EXEEXT)
$(PROGSDIR)/siteom$(EXEEXT): SiteOmega.o $(OBJS)
	$(CC) SiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

diffsel$(EXEEXT): $(PROGSDIR)/diffsel$(EXEEXT)
$(PROGSDIR)/diffsel$(EXEEXT): DiffSel.o $(OBJS)
	$(CC) DiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

