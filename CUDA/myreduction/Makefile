SHELL=/bin/sh
FC=gfortran
FCLAGS=-O3 -fno-math-errno -fPIC
CC=gcc
CPP=g++
CCFLAGS=-O3 -fno-math-errno -mfloat-abi=hard
NVCC=nvcc -ccbin $(CC)
#NVCCFLAGS=-O2 -m32 -Xcompiler -mfloat-abi=hard -gencode arch=compute_32,code=sm_32 # --use_fast_math
NVCCFLAGS=-O2 -m32 -Xcompiler -mfloat-abi=hard -gencode arch=compute_32,code=sm_32 --use_fast_math
HEADS=kernel.cu
LD=$(FC)
DEFS+=-D__FTYPE=real*4 -D__CTYPE=float -D_FMT="%12.5f\n"
#DEFS+=-D__FTYPE=integer -D__CTYPE=int
#DEFS+=-D__TEX
#LD=$(NVCC)
#LD=$(CC)

LIBS+=-L/usr/local/cuda/lib -lcudart -lcuda -lgfortran -lstdc++ -lgcc -lm
#LIBS+=-L/usr/local/cuda/lib -lcudart -lcuda

CUINCLUDES+=-I/usr/local/cuda/samples/common/inc

NAME=reduction

CUOBJS=reduction_c.o
FOBJS=reduction.o

EXE=$(NAME)
OBJS=$(CUOBJS) $(FOBJS) $(COBJS)

default: $(EXE)

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)

.SUFFIXES: .o .cu .cpp

%.o:%.cu $(HEADS)
	$(NVCC) $(DEFS) $(NVCCFLAGS) $(CUINCLUDES) -c $< -o $@

%.o:%.c $(HEADS)
	$(CC) $(DEFS) $(CCFLAGS) -c $< -o $@

%.o:%.cpp $(HEADS)
	$(CPP) $(DEFS) $(CCFLAGS) -c $< -o $@

%.o:%.F90 $(HEADS)
	$(FC) $(DEFS) $(FCLAGS) -c $< -o $@


clean::
	rm -f $(OBJS) $(EXE)
