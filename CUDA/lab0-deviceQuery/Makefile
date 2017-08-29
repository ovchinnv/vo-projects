SHELL=/bin/sh
CC=gcc
NVCC=nvcc -ccbin $(CC)
NVCCFLAGS=-O3 -gencode arch=compute_32,code=sm_32
CUHEADS=cutil_inline.h
LD=$(NVCC)

CUINCLUDES+=-I/usr/local/cuda/samples/common/inc

NAME=deviceQuery

CUOBJS=deviceQuery.o

CUSRC=$(CUOBJS:%.o=%.cu)

EXE=$(NAME)
OBJS=$(CUOBJS)

default: $(EXE)

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXE)

.SUFFIXES: .o .cu

%.o:%.cu $(CUHEADS)
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -c $< -o $@


clean::
	rm -f $(OBJS) $(EXE)
