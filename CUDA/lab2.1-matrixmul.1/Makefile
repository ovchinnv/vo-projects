SHELL=/bin/sh
CC=gcc
CPP=g++
CCFLAGS=-O3 -fno-math-errno -mfloat-abi=hard
NVCC=nvcc -ccbin $(CC)
NVCCFLAGS=-O3 -m32 -Xcompiler -mfloat-abi=hard -gencode arch=compute_32,code=sm_32 --use_fast_math
HEADS=assist.h matrixmul.h matrixmul_kernel.cu
LD=$(NVCC)

CUINCLUDES+=-I/usr/local/cuda/samples/common/inc

NAME=matmul

CUOBJS=matrixmul.o
COBJS=matrixmul_gold.o

CUSRC=$(CUOBJS:%.o=%.cu)
CSRC=$(COBJS:%.o=%.cu)

EXE=$(NAME)
OBJS=$(CUOBJS) $(COBJS)

default: $(EXE)

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXE)

.SUFFIXES: .o .cu .cpp

%.o:%.cu $(HEADS)
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -c $< -o $@

%.o:%.cpp $(HEADS)
	$(CPP) $(CCFLAGS) -c $< -o $@


clean::
	rm -f $(OBJS) $(EXE)
