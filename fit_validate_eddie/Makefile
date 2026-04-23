SHELL = /bin/sh
UNAME = $(shell uname)

ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX          = g++
RM           = rm -f
AR           = ar cru

CXXFLAGS     = -O3 -g -fPIC -funroll-loops

EXENAME     = FastOpticalFitValidate
SRCEXT      = cpp
SRCDIR      = src
INCDIR      = include
OBJDIR      = build
EXEDIR      = bin
SRCS        := $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
OBJS        := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

GARBAGE     = $(OBJDIR)/*.o $(EXEDIR)/$(EXENAME)

ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS)
LINKFLAGS    = -g $(shell root-config --nonew) $(shell root-config --ldflags) -Wl,--no-as-needed
endif

ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS)
LINKFLAGS    =
endif

LIBS       += $(ROOTLIBS) -lHtml -lThread -lMathMore

all : $(EXEDIR)/$(EXENAME)

$(EXEDIR)/$(EXENAME) : $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LINKFLAGS) $(LIBS)

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE)
