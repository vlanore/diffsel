CPPFLAGS= -std=c++11 -g -Wall -Wextra -Wpedantic -O3

SRCS= BranchSitePath.cpp Chrono.cpp CodonSequenceAlignment.cpp CodonStateSpace.cpp CodonSubMatrix.cpp GTRSubMatrix.cpp MSCodonSubMatrix.cpp PhyloProcess.cpp Random.cpp SequenceAlignment.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp linalg.cpp linalg2.cpp Chain.cpp Sample.cpp
OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_CPP=$(wildcard *.cpp)
ALL_SRCS = $(wildcard *.cpp) $(wildcard *.hpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_CPP))

PROGSDIR=..
ALL= flatdiffsel readflatdiffsel flatglobom flatgeneglobom
EXEEXT=_bin
PROGS=$(addsuffix $(EXEEXT), $(addprefix $(PROGSDIR)/, $(ALL)))

.PHONY: all clean format
all: $(PROGS)

# Rules to generate the dependencies automatically
%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$

# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CXX) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets
flatglobom$(EXEEXT): $(PROGSDIR)/flatglobom$(EXEEXT)
$(PROGSDIR)/flatglobom$(EXEEXT): SingleOmega.o $(OBJS)
	$(CXX) SingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

flatgeneglobom$(EXEEXT): $(PROGSDIR)/flatgeneglobom$(EXEEXT)
$(PROGSDIR)/flatgeneglobom$(EXEEXT): MultiGeneSingleOmega.o $(OBJS)
	$(CXX) MultiGeneSingleOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

flatdiffsel$(EXEEXT): $(PROGSDIR)/flatdiffsel$(EXEEXT)
$(PROGSDIR)/flatdiffsel$(EXEEXT): DiffSel.o $(OBJS)
	$(CXX) DiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readflatdiffsel$(EXEEXT): $(PROGSDIR)/readflatdiffsel$(EXEEXT)
$(PROGSDIR)/readflatdiffsel$(EXEEXT): ReadDiffSel.o $(OBJS)
	$(CXX) ReadDiffSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

format:
	clang-format -i $(ALL_SRCS)

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

