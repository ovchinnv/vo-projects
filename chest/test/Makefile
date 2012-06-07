ROOT=../
include $(ROOT)/Makefile.inc

#FC=gfortran
#FFLAGS=-O2

TESTS=\
 test2d-1\
 test2d-2\
 test2d-3\
 test2d-4\
 test3d-1\
 test3d-2

SCRATCH=\
	xbc0_test.dat\
	xbc1_test.dat\
	ybc0_test.dat\
	ybc1_test.dat\
	zbc0_test.dat\
	zbc1_test.dat\
	xg_test.dat\
	yg_test.dat\
	zg_test.dat\
	test.dat\
	xy_test.xyz\
	xyz_test.xyz\
	solution.dat\
	rhs_test.dat\
	uexact_test.dat\
	rhs_test.dat\
	test.F90\
	compare.F90

default : all
all : print $(TESTS)
$(TESTS) : force_look
	./$@ "$(FC)" "$(FFLAGS)" "$(DEFS)">& $@.out ; echo " $@ : " `grep " TEST " $@.out`

force_look:
	true

print:
	echo
	echo " =========================== "
	echo " RUNNING ALL AVAILABLE TESTS "
	echo " =========================== "
	pass=1

.SILENT : force_look print clean $(TESTS)

clean:
	rm -f $(SCRATCH)

