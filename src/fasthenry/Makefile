# Default flags (add more specific ones based on the OS below)
CFLAGS := -O -DFOUR -m64 -fcommon -std=gnu89 -DNO_GETHOSTNAME

# Compiler path for different systems
ifeq ($(OS), Windows_NT)
    # Windows-specific settings
    CC := gcc  # Use the appropriate path to gcc if needed
    CFLAGS += -DWINDOWS

else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S), Linux)
        # Linux-specific settings
        CC := gcc  # Replace with the correct path if not default
        CFLAGS += -DLINUX
    endif
    ifeq ($(UNAME_S), Darwin)
        # macOS-specific settings
        CC := /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gcc-14  # Replace if necessary
        CFLAGS += -DMACOS
    endif
endif

# Other variables
SHELL = /bin/sh
RM = /bin/rm
BIN = ../../bin

SPLIB = sparse/sparse.a

# Object files
OBJS = induct.o gmres.o savemat_mod.o readGeom.o joelself.o writefastcap.o \
       SetupMulti.o calcp.o SetupComputePsi.o mulSetup.o BreakupSeg.o \
       Precond.o addgroundplane.o findpaths.o fillM.o cx_ludecomp.o \
       parse_command_line.o default_opts.o Prec_cost.o hole.o \
       dist_betw_fils.o mutual.o newPrecond.o deg_mutual.o \
       barnoldi.o regurgitate.o

MOBJS = mulGlobal.o mulDo.o mulMulti.o mulLocal.o mulMats.o direct.o \
        uglieralloc.o capsolve.o 

NONUNIOBJS = find_nonuni_path.o read_tree.o contact.o

HEADER = induct.h cmplx.h resusage.h
MULHEAD = mulStruct.h mulGlobal.h patran.h resusage.h
NONUNIHEAD = gp.h

# Main target
fasthenry: $(OBJS) $(MOBJS) $(NONUNIOBJS) $(SPLIB)
	$(CC) -o fasthenry $(CFLAGS) $(OBJS) $(MOBJS) $(NONUNIOBJS) $(SPLIB) -lm
	mv fasthenry $(BIN)/fasthenry

# Sparse library target
$(SPLIB):
	cd sparse; make

# Clean rule
clean:
	$(RM) -f *.o
	cd sparse; make clean

# Dependency rules
$(OBJS): $(HEADER) $(MULHEAD)
$(MOBJS): $(MULHEAD)
$(NONUNIOBJS): $(HEADER) $(MULHEAD) $(NONUNIHEAD)