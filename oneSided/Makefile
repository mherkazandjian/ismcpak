MPICXX ?= mpicxx
MPICC ?= mpicc

CFLAGS   += -Wall -g -O3 -I$(PREFIX)/include
CXXFLAGS += $(CFLAGS) -L$(PREFIX)/lib
LDFLAGS  += -lm $(MUSE_LD_FLAGS) 

OBJS = interface.o

CODELIB = src/libpdr.a

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

export AMUSE_PREFIX = $(PREFIX)
#export IS_AMUSE_WORKER = 1

all: worker __init__.py

clean:
	$(RM) -f *.so *.o *.pyc worker.cc worker.h 
	$(RM) *~ worker
	make -C src clean

$(CODELIB):
	make -C src amuse_interface

worker.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py pdrInterface -o $@

worker.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py pdrInterface -o $@

worker: worker.cc worker.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $@.cc $(OBJS) $(CODELIB) -o $@

.c.o: $<
	$(MPICC) $(CFLAGS) -c -o $@ $< 

__init__.py: 
	touch $@