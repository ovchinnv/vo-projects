ROOT=../../..
include $(ROOT)/Makefile.inc

#FC=gfortran
#FFLAGS=-O2

TESTS=\
 test2d-1\
 test2d-2\
 test2d-3\
 test2d-4\
 test2d-5\
 test2d-6\
 test2d-7\
 test2d-8\
 test2d-9\
 test3d-1\
 test3d-2\
 test3d-3\
 test3d-4

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

all : print $(TESTS)
$(TESTS) : force_look
	if [ ! -f ../ches ]; then ln -s ../../driver/muphen ../ches ; fi
	./$@ "$(FC)" "$(FFLAGS)" "$(DEFS)">& $@.out ; echo " $@ : " `grep " TEST " $@.out`
#	./$@ "$(FC)" "$(FFLAGS)" "$(DEFS)">& $@.out ; echo " $@ : " `grep " TEST " $@.out`


print:
	echo
	echo " =========================== "
	echo " RUNNING ALL AVAILABLE TESTS "
	echo " =========================== "
	pass=1

.SILENT : force_look print clean $(TESTS)

clean::
	rm -f $(SCRATCH)
