CC = gcc

MACH=$(HOSTNAME)

#ICC stuff
#----------
ifeq ($(CC),icc)

 OPTS += -vec-report1
 #OPTS += -fp-model precise
 OPTS += -align

 #ruineraa
 ifeq ($(MACH),ruineraa)
   CC = /software/local/intel/compiler91/bin/icc
   LIBS += -L/software/local/intel/compiler91/lib 
   OPTS += -msse3
 endif

 #particle2
 ifeq ($(MACH),particle2)
  CC = $(PROGSDIR)/intel/bin/intel64/icc
  LIBS += -L$(PROGSDIR)/intel/lib64
  OPTS += -msse4.1
  OPTS += -axSSE4.1
  OPTS += -msse4.1
  LIBS += -L$(PROGSDIR)/intel/lib/intel64 -lintlc
 endif

  LIBS += -limf
else
 LIBS += -lm
endif

OPTS += -O3
OPTS += -g
OPTS += -Wall
#OPTS += -ansi
OPTS += -std=gnu99
#OPT += -DHAVEGSL
#OPTS += -DENABEL_MPI
#OPTS += -DFIRST_SLAB_ONLY   # if this is true, then only the equilibrium of the first surface slab is computed 

#automatically adds the CFLAG when compiled through the makefile of AMUSE
#ifdef IS_AMUSE_WORKER
#	OPTS += -DIS_AMUSE_WORKER
#endif

ifeq (-DENABEL_MPI,$(findstring -DENABEL_MPI,$(OPTS)))

    #setting MPI -L and -I paths
    ifdef AMUSE_PREFIX
    MPICH2DIR=$(AMUSE_PREFIX)
    endif

    LIBS += -L$(MPICH2DIR)/lib -lmpich -lgcc_s -lpthread -lrt
    OPTS += -I$(MPICH2DIR)/include
endif


EXEC = chemistry

ifeq (-DHAVEGSL,$(findstring -DHAVEGSL,$(OPT)))
  ifeq ($(MACH), particle2)
    LIBS += -L$(PROGSDIR)/gsl-1.12-mod-0/lib 
    OPTS += -I$(PROGSDIR)/gsl-1.12-mod-0/include  
  endif
endif

INCL = constants.h  definitions.h  prototypes.h  vars.h

ifeq (-DHAVEGSL,$(findstring -DHAVEGSL,$(OPT)))
  LIBS += -lgsl -lgslcblas
endif

CFLAGS = $(OPTS) 

OBJ = init.o vars.o read_reaction_data.o chemical_balance.o finestruct_cooling.o\
      Matrix_inversion.o metastable_cooling.o rovibrational_cooling.o self_shielding.o\
      thermal_balance.o heating.o util.o pfileIO.o read_self_shielding.o cooling.o \
      database.o

LIBOBJ = init.o vars.o read_reaction_data.o chemical_balance.o finestruct_cooling.o\
         Matrix_inversion.o metastable_cooling.o rovibrational_cooling.o self_shielding.o\
         thermal_balance.o heating.o util.o pfileIO.o read_self_shielding.o cooling.o\
         database.o

#----------------------------------------------------------------------------------------------------
.PHONY : all clean nothing
nothing	:
#----------------------------------------------------------------------------------------------------

$(EXEC) : $(OBJ) Makefile $(INCL)
	$(CC) main.c $(OBJ) $(LIBS) $(CFLAGS) -o $(EXEC)

$(OBJ): $(INCL) $(INCL)

amuse_interface: libpdr.a

libpdr.a: $(LIBOBJ)
	ar crs $@ $^ 

all: $(EXEC) $(VERSION) Makefile
	@echo "MACHINE             : " $(MACH)
ifeq (particle2,$(findstring particle2,$(MACH)))
	   @git-getBranchName.sh > VERSION
	   @echo "Dumped Version file : " $$(cat VERSION) 
endif

clean:
	rm -fr *.o *.a $(EXEC) *~ 

