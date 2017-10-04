SHELL=/bin/sh
FC=gfortran
FCLAGS=-O3 -fno-math-errno -fPIC
CC=gcc
CPP=g++
CCFLAGS=-O3 -fno-math-errno -mfloat-abi=hard
NVCC=nvcc -ccbin $(CC)
NVCCFLAGS=-O2 -m32 -Xcompiler -mfloat-abi=hard -gencode arch=compute_32,code=sm_32 --use_fast_math
HEADS=matrixmul_kernel.cu interface.h
LD=$(FC)
FPP=cpp
FPPFLAGS = -P -C -nostdinc
#LD=$(NVCC)
#LD=$(CC)

DEFS+=-DBLOCK_SIZE_X=8 -DBLOCK_SIZE_Y=8 -D__IMPNONE='implicit none' -D__FCINT=integer -D__FCFLOAT=real*4 -D__CUFLOAT=float -D__CFLOAT=float -D__CINT=int

LIBS+=-L/usr/local/cuda/lib -lcudart -lcuda -lgfortran -lstdc++ -lgcc -lm
#LIBS+=-L/usr/local/cuda/lib -lcudart -lcuda

CUINCLUDES+=-I/usr/local/cuda/samples/common/inc

NAME=matmul

CUOBJS=matrixmul.o CUmgutil.o
FOBJS=matrix.o
#COBJS=matmul_c.o

EXE=$(NAME)
OBJS=$(CUOBJS) $(FOBJS) $(COBJS)

default: $(EXE)

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)

.SUFFIXES: .o .cu .cpp

%.h:%.H $(LOCALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) $(LOCALMACROS) -P $< > $@

%.o:%.cu $(HEADS)
	$(NVCC) $(DEFS) $(NVCCFLAGS) $(CUINCLUDES) -c $< -o $@

%.o:%.c $(HEADS)
	$(CC) $(CCFLAGS) -c $< -o $@

%.o:%.cpp $(HEADS)
	$(CPP) $(CCFLAGS) -c $< -o $@

%.o:%.f90 $(HEADS)
	$(FC) $(FCLAGS) -c $< -o $@


clean::
	rm -f $(OBJS) $(EXE)
