for d in bestfit constants dynamol erf files lu multicom multidiag output parselist parser random string timer vectors chest; do (cd $d; make); done
make[1]: Entering directory `/home/taly/projects/bestfit'
for file in bestfit.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../bestfit.stamp
make[1]: Leaving directory `/home/taly/projects/bestfit'
make[1]: Entering directory `/home/taly/projects/constants'
for file in constants.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../constants.stamp
make[1]: Leaving directory `/home/taly/projects/constants'
make[1]: Entering directory `/home/taly/projects/dynamol'
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I ../include -I /usr/local/mpich2-path64/include stats.F90 -o stats.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for stats_
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I ../include -I /usr/local/mpich2-path64/include system.F90 -o system.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead

 function system_getind(selection) result(ind)
                                          ^    
pathf95-287 pathf90: WARNING SYSTEM_GETIND, File = system.F90, Line = 884, Column = 43 
  The result of function name "IND" in the function subprogram is not defined.

pathf95: PathScale(TM) Fortran Version 1.0.0 (f14) Fri Jul 27, 2012  21:31:44
pathf95: 931 source lines
pathf95: 0 Error(s), 1 Warning(s), 0 Other message(s), 0 ANSI(s)
pathf95: "explain pathf95-message number" gives more information about each message
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for system_
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I ../include -I /usr/local/mpich2-path64/include sysinfo.F90 -o sysinfo.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Assembly: NULL DST passed to CG for sysinfo_
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I ../include -I /usr/local/mpich2-path64/include sysmanip.F90 -o sysmanip.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for sysmanip_
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I ../include -I /usr/local/mpich2-path64/include verlet.F90 -o verlet.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for verlet_
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I ../include -I /usr/local/mpich2-path64/include molsim.F90 -o molsim.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for molsim_
for file in stats.o sysinfo.o system.o sysmanip.o verlet.o molsim.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../dynamol.stamp
make[1]: Leaving directory `/home/taly/projects/dynamol'
make[1]: Entering directory `/home/taly/projects/erf'
for file in erfappx.o erfsun_2_f90.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../erf.stamp
make[1]: Leaving directory `/home/taly/projects/erf'
make[1]: Entering directory `/home/taly/projects/files'
for file in files.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../files.stamp
make[1]: Leaving directory `/home/taly/projects/files'
make[1]: Entering directory `/home/taly/projects/lu'
cpp -P -C -I .. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../source.defs   -include ../source.msg -P lu.ftn > lu.F90
lu.ftn:8:4: warning: missing terminating " character [enabled by default]
lu.ftn:10:27: warning: missing terminating " character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I /usr/local/mpich2-path64/include lu.F90 -o lu.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for lu_
for file in lu.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../lu.stamp
make[1]: Leaving directory `/home/taly/projects/lu'
make[1]: Entering directory `/home/taly/projects/multicom'
for file in multicom.o multicom_aux.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../multicom.stamp
make[1]: Leaving directory `/home/taly/projects/multicom'
make[1]: Entering directory `/home/taly/projects/multidiag'
cpp -P -C -I .. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../source.defs   -include ../source.msg -P multidiag.ftn > multidiag.F90
multidiag.ftn:22:125: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I /usr/local/mpich2-path64/include multidiag.F90 -o multidiag.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for multidiag_
for file in multidiag.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../multidiag.stamp
make[1]: Leaving directory `/home/taly/projects/multidiag'
make[1]: Entering directory `/home/taly/projects/output'
for file in output.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../output.stamp
make[1]: Leaving directory `/home/taly/projects/output'
make[1]: Entering directory `/home/taly/projects/parselist'
for file in parselist.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../parselist.stamp
make[1]: Leaving directory `/home/taly/projects/parselist'
make[1]: Entering directory `/home/taly/projects/parser'
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../modules -I ../include -I /usr/local/mpich2-path64/include parser.F90 -o parser.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_U4) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_U4) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_U4) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for parser_
for file in parser.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../parser.stamp
make[1]: Leaving directory `/home/taly/projects/parser'
make[1]: Entering directory `/home/taly/projects/random'
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../modules -I ../include -I /usr/local/mpich2-path64/include clcg.f -o clcg.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../modules -I ../include -I /usr/local/mpich2-path64/include rng.F90 -o rng.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for rng_
for file in rng.o clcg.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../random.stamp
make[1]: Leaving directory `/home/taly/projects/random'
make[1]: Entering directory `/home/taly/projects/string'
cpp -P -C -I .. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../source.defs  -include source.defs -include ../source.msg -P sm_var.ftn > sm_var.F90
DIR=string/sm0k.stamp ; DIR=../${DIR%.*} ; if [ -d $DIR ] ; then cd $DIR; make ; fi
make[2]: Entering directory `/home/taly/projects/string/sm0k'
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include ../source.defs -include ../../source.msg -P sm0k.ftn > sm0k.F90
sm0k.ftn:524:69: warning: missing terminating ' character [enabled by default]
sm0k.ftn:629:44: warning: missing terminating ' character [enabled by default]
sm0k.ftn:924:24: warning: missing terminating ' character [enabled by default]
sm0k.ftn:1753:50: warning: missing terminating ' character [enabled by default]
sm0k.ftn:2288:41: warning: missing terminating ' character [enabled by default]
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include ../source.defs -include ../../source.msg -P ../sm_config.ftn > ../sm_config.F90
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include ../sm_config.F90 -o ../sm_config.o
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for sm_config_
file=sm_config.mod ; \
if [ -f $file ] ; then mv -f $file ../$file 2>/dev/null; touch -m ../sm_config.mod 2>/dev/null ; else \
FILE=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ; \
if [ -f $FILE ] ; then mv -f $FILE ../$FILE 2>/dev/null; \
touch -m ../$FILE 2>/dev/null ; \
rm -f ../$file ; \
ln -s $FILE ../$file ; fi ; fi
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include sm0k.F90 -o sm0k.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for sm0k_
for file in sm0k.o; do \
ln -fs `pwd`/$file ../../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../../include/$file; fi ;\
done; 
echo -n > ../../string/sm0k.stamp
rm ../sm_config.F90 ../sm_config.o
make[2]: Leaving directory `/home/taly/projects/string/sm0k'
DIR=string/smcv.stamp ; DIR=../${DIR%.*} ; if [ -d $DIR ] ; then cd $DIR; make ; fi
make[2]: Entering directory `/home/taly/projects/string/smcv'
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_types.ftn > cv_types.F90
cv_types.ftn:7:42: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_types.F90 -o cv_types.o
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for cv_types_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_common.ftn > cv_common.F90
cv_common.ftn:2555:20: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include ../sm_var.F90 -o ../sm_var.o
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for sm_var_
file=sm_var.mod ; \
if [ -f $file ] ; then mv -f $file ../$file 2>/dev/null; touch -m ../sm_var.mod 2>/dev/null ; else \
FILE=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ; \
if [ -f $FILE ] ; then mv -f $FILE ../$FILE 2>/dev/null; \
touch -m ../$FILE 2>/dev/null ; \
rm -f ../$file ; \
ln -s $FILE ../$file ; fi ; fi
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_common.F90 -o cv_common.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_common_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_frames.ftn > cv_frames.F90
cv_frames.ftn:19:113: warning: missing terminating ' character [enabled by default]
cv_frames.ftn:25:64: warning: missing terminating ' character [enabled by default]
cv_frames.ftn:242:93: warning: missing terminating ' character [enabled by default]
cv_frames.ftn:379:45: warning: missing terminating ' character [enabled by default]
cv_frames.ftn:495:24: warning: missing terminating ' character [enabled by default]
cv_frames.ftn:862:23: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_frames.F90 -o cv_frames.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Assembly: NULL DST passed to CG for cv_frames_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_cvrms.ftn > cv_cvrms.F90
cv_cvrms.ftn:4:46: warning: missing terminating ' character [enabled by default]
cv_cvrms.ftn:10:28: warning: missing terminating ' character [enabled by default]
cv_cvrms.ftn:117:43: warning: missing terminating ' character [enabled by default]
cv_cvrms.ftn:120:69: warning: missing terminating ' character [enabled by default]
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_posi_com.ftn > cv_posi_com.F90
cv_posi_com.ftn:133:43: warning: missing terminating ' character [enabled by default]
cv_posi_com.ftn:136:69: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_posi_com.F90 -o cv_posi_com.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_U4) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_posi_com_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_dist_com.ftn > cv_dist_com.F90
cv_dist_com.ftn:4:49: warning: missing terminating ' character [enabled by default]
cv_dist_com.ftn:9:28: warning: missing terminating ' character [enabled by default]
cv_dist_com.ftn:99:43: warning: missing terminating ' character [enabled by default]
cv_dist_com.ftn:102:69: warning: missing terminating ' character [enabled by default]
cv_dist_com.ftn:119:128: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_dist_com.F90 -o cv_dist_com.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_dist_com_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_angle_com.ftn > cv_angle_com.F90
cv_angle_com.ftn:5:46: warning: missing terminating ' character [enabled by default]
cv_angle_com.ftn:13:28: warning: missing terminating ' character [enabled by default]
cv_angle_com.ftn:106:43: warning: missing terminating ' character [enabled by default]
cv_angle_com.ftn:109:69: warning: missing terminating ' character [enabled by default]
cv_angle_com.ftn:213:82: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_angle_com.F90 -o cv_angle_com.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_angle_com_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_anglvec.ftn > cv_anglvec.F90
cv_anglvec.ftn:4:46: warning: missing terminating ' character [enabled by default]
cv_anglvec.ftn:17:28: warning: missing terminating ' character [enabled by default]
cv_anglvec.ftn:203:43: warning: missing terminating ' character [enabled by default]
cv_anglvec.ftn:359:82: warning: missing terminating ' character [enabled by default]
cv_anglvec.ftn:544:101: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_anglvec.F90 -o cv_anglvec.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_anglvec_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_dihe_com.ftn > cv_dihe_com.F90
cv_dihe_com.ftn:4:55: warning: missing terminating ' character [enabled by default]
cv_dihe_com.ftn:9:28: warning: missing terminating ' character [enabled by default]
cv_dihe_com.ftn:103:43: warning: missing terminating ' character [enabled by default]
cv_dihe_com.ftn:106:69: warning: missing terminating ' character [enabled by default]
cv_dihe_com.ftn:269:43: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_dihe_com.F90 -o cv_dihe_com.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_dihe_com_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_rmsd.ftn > cv_rmsd.F90
cv_rmsd.ftn:4:45: warning: missing terminating ' character [enabled by default]
cv_rmsd.ftn:9:28: warning: missing terminating ' character [enabled by default]
cv_rmsd.ftn:50:79: warning: missing terminating ' character [enabled by default]
cv_rmsd.ftn:167:43: warning: missing terminating ' character [enabled by default]
cv_rmsd.ftn:170:69: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_rmsd.F90 -o cv_rmsd.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_rmsd_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_drmsd.ftn > cv_drmsd.F90
cv_drmsd.ftn:4:46: warning: missing terminating ' character [enabled by default]
cv_drmsd.ftn:5:48: warning: missing terminating ' character [enabled by default]
cv_drmsd.ftn:10:28: warning: missing terminating ' character [enabled by default]
cv_drmsd.ftn:93:79: warning: missing terminating ' character [enabled by default]
cv_drmsd.ftn:212:43: warning: missing terminating ' character [enabled by default]
cv_drmsd.ftn:215:69: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_drmsd.F90 -o cv_drmsd.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_drmsd_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_proj.ftn > cv_proj.F90
cv_proj.ftn:4:45: warning: missing terminating ' character [enabled by default]
cv_proj.ftn:9:28: warning: missing terminating ' character [enabled by default]
cv_proj.ftn:41:43: warning: missing terminating ' character [enabled by default]
cv_proj.ftn:44:69: warning: missing terminating ' character [enabled by default]
cv_proj.ftn:241:30: warning: missing terminating ' character [enabled by default]
cv_proj.ftn:313:30: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_proj.F90 -o cv_proj.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_proj_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_quaternion.ftn > cv_quaternion.F90
cv_quaternion.ftn:15:93: warning: missing terminating ' character [enabled by default]
cv_quaternion.ftn:21:64: warning: missing terminating ' character [enabled by default]
cv_quaternion.ftn:356:26: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_quaternion.F90 -o cv_quaternion.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_quaternion_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P cv_qcomp.ftn > cv_qcomp.F90
cv_qcomp.ftn:9:28: warning: missing terminating ' character [enabled by default]
cv_qcomp.ftn:143:43: warning: missing terminating ' character [enabled by default]
cv_qcomp.ftn:146:69: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_qcomp.F90 -o cv_qcomp.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_U4) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_qcomp_
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include cv_cvrms.F90 -o cv_cvrms.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for cv_cvrms_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P smcv_master.ftn > smcv_master.F90
smcv_master.ftn:483:52: warning: missing terminating ' character [enabled by default]
smcv_master.ftn:688:107: warning: missing terminating ' character [enabled by default]
smcv_master.ftn:1066:113: warning: missing terminating ' character [enabled by default]
smcv_master.ftn:2042:113: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include smcv_master.F90 -o smcv_master.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for smcv_master_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P smcv.ftn > smcv.F90
smcv.ftn:544:12: warning: missing terminating ' character [enabled by default]
smcv.ftn:697:79: warning: missing terminating ' character [enabled by default]
smcv.ftn:1952:48: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include smcv.F90 -o smcv.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P smcv_add.ftn > smcv_add.F90
smcv_add.ftn:984:69: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../../include -I .. -I /usr/local/mpich2-path64/include smcv_add.F90 -o smcv_add.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
for file in cv_types.o cv_common.o cv_frames.o cv_cvrms.o smcv_master.o smcv.o smcv_add.o cv_posi_com.o cv_dist_com.o cv_angle_com.o cv_anglvec.o cv_dihe_com.o cv_rmsd.o cv_drmsd.o cv_proj.o cv_quaternion.o cv_qcomp.o; do \
ln -fs `pwd`/$file ../../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../../include/$file; fi ;\
done; 
echo -n > ../../string/smcv.stamp
rm ../sm_var.o
make[2]: Leaving directory `/home/taly/projects/string/smcv'
DIR=string/ftsm.stamp ; DIR=../${DIR%.*} ; if [ -d $DIR ] ; then cd $DIR; make ; fi
make[2]: Entering directory `/home/taly/projects/string/ftsm'
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P ftsm_var.ftn > ftsm_var.F90
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I .. -I ../../include -I .. -I /usr/local/mpich2-path64/include ftsm_var.F90 -o ftsm_var.o
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for ftsm_var_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P ftsm_rex.ftn > ftsm_rex.F90
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I .. -I ../../include -I .. -I /usr/local/mpich2-path64/include ftsm_rex.F90 -o ftsm_rex.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for ftsm_rex_
cpp -P -C -I ../.. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../../source.defs  -include source.defs -include ../source.defs -include ../../source.msg -P ftsm.ftn > ftsm.F90
ftsm.ftn:3306:30: warning: missing terminating ' character [enabled by default]
ftsm.ftn:3384:30: warning: missing terminating ' character [enabled by default]
ftsm.ftn:3616:128: warning: missing terminating ' character [enabled by default]
ftsm.ftn:4217:75: warning: missing terminating ' character [enabled by default]
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I .. -I ../../include -I .. -I /usr/local/mpich2-path64/include ftsm.F90 -o ftsm.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for ftsm_
for file in ftsm_var.o ftsm_rex.o ftsm.o; do \
ln -fs `pwd`/$file ../../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../../include/$file; fi ;\
done; 
echo -n > ../../string/ftsm.stamp
make[2]: Leaving directory `/home/taly/projects/string/ftsm'
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I /usr/local/mpich2-path64/include sm_var.F90 -o sm_var.o
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for sm_var_
cpp -P -C -I .. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../source.defs  -include source.defs -include ../source.msg -P sm_config.ftn > sm_config.F90
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I /usr/local/mpich2-path64/include sm_config.F90 -o sm_config.o
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for sm_config_
cpp -P -C -I .. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../source.defs  -include source.defs -include ../source.msg -P splines.ftn > splines.F90
splines.ftn:354:11: warning: missing terminating ' character [enabled by default]
splines.ftn:810:28: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I /usr/local/mpich2-path64/include splines.F90 -o splines.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
cpp -P -C -I .. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../source.defs  -include source.defs -include ../source.msg -P sm_util.ftn > sm_util.F90
sm_util.ftn:26:96: warning: missing terminating ' character [enabled by default]
sm_util.ftn:130:97: warning: missing terminating ' character [enabled by default]
sm_util.ftn:326:97: warning: missing terminating ' character [enabled by default]
sm_util.ftn:698:97: warning: missing terminating ' character [enabled by default]
sm_util.ftn:1037:97: warning: missing terminating ' character [enabled by default]
sm_util.ftn:1232:102: warning: missing terminating ' character [enabled by default]
sm_util.ftn:1332:46: warning: missing terminating ' character [enabled by default]
sm_util.ftn:1564:20: warning: missing terminating ' character [enabled by default]
sm_util.ftn:1739:20: warning: missing terminating ' character [enabled by default]
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I /usr/local/mpich2-path64/include sm_util.F90 -o sm_util.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- Create MU and CHI list: Alias analysis: f90 pointer meets non-f90 pointer
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
cpp -P -C -I .. -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -U__CHARMM -include ../source.defs  -include source.defs -include ../source.msg -P sm_main.ftn > sm_main.F90
pathf90 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../include -I /usr/local/mpich2-path64/include sm_main.F90 -o sm_main.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
for file in sm_var.o sm_config.o splines.o sm_util.o sm_main.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../string.stamp
make[1]: Leaving directory `/home/taly/projects/string'
make[1]: Entering directory `/home/taly/projects/timer'
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP  -I ../modules -I ../include -I /usr/local/mpich2-path64/include timer.F90 -o timer.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for timer_
for file in timer.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../timer.stamp
make[1]: Leaving directory `/home/taly/projects/timer'
make[1]: Entering directory `/home/taly/projects/vectors'
for file in ivector.o rvector.o ivector_list.o rvector_list.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../vectors.stamp
make[1]: Leaving directory `/home/taly/projects/vectors'
make[1]: Entering directory `/home/taly/projects/chest'
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -D SIZE=gridsize -D __TRUNCATESUMS -D __erf=erfo7 -D __TIMER  -I ../include -I . -I /usr/local/mpich2-path64/include gridsize.F90 -o gridsize.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for gridsize_
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -D SIZE=gridsize -D __TRUNCATESUMS -D __erf=erfo7 -D __TIMER  -I ../include -I . -I /usr/local/mpich2-path64/include state.F90 -o state.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for state_
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -D SIZE=gridsize -D __TRUNCATESUMS -D __erf=erfo7 -D __TIMER  -I ../include -I . -I /usr/local/mpich2-path64/include bc.F90 -o bc.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for bc_
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -D SIZE=gridsize -D __TRUNCATESUMS -D __erf=erfo7 -D __TIMER  -I ../include -I . -I /usr/local/mpich2-path64/include grid.F90 -o grid.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for grid_
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -D SIZE=gridsize -D __TRUNCATESUMS -D __erf=erfo7 -D __TIMER  -I ../include -I . -I /usr/local/mpich2-path64/include multigrid.F90 -o multigrid.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_U4) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Assembly: NULL DST passed to CG for multigrid_
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
!!! DevWarn during Global Optimization -- MainOpt emitter: PREG (.preg_I8) has mismatching MTYPE-size and TY-size; refer to bug #567932
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -D SIZE=gridsize -D __TRUNCATESUMS -D __erf=erfo7 -D __TIMER  -I ../include -I . -I /usr/local/mpich2-path64/include PBmain.F90 -o PBmain.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
!!! DevWarn during Assembly: NULL DST passed to CG for pbmain_
pathf95 -extend-source -OPT:Olimit=8753 -INLINE -c -D int=integer -D float=real*8 -D bool=logical -D int4mpi=integer*4 -D mpiint=MPI_INTEGER -D mpifloat=MPI_REAL -D mpichar=MPI_CHARACTER -D mpiint4=MPI_INTEGER4 -D mpibool=MPI_LOGICAL -D __DMOL -D __RCOMP -D SIZE=gridsize -D __TRUNCATESUMS -D __erf=erfo7 -D __TIMER  -I ../include -I . -I /usr/local/mpich2-path64/include ches.F90 -o ches.o
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during IR->WHIRL Conversion: Count limit reached on the following DevWarn:
!!! DevWarn during IR->WHIRL Conversion: Should use ST_pu_type instead
!!! DevWarn during Reading WHIRL file: TODO: implement *skip* option
!!! DevWarn during Reading WHIRL file: IPA_NODE::Skip is not yet implemented
for file in gridsize.o state.o bc.o grid.o multigrid.o PBmain.o ches.o; do \
ln -fs `pwd`/$file ../obj/$file; \
file=${file%.*}.mod ; if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
file=`echo ${file%.*} | tr '[:lower:]' '[:upper:]'`.mod ;\
if [ -f $file  ] ; then ln -fs `pwd`/$file ../include/$file; fi ;\
done; 
echo -n > ../chest.stamp
pathf95  gridsize.o state.o bc.o grid.o multigrid.o PBmain.o ches.o ../obj/output.o ../obj/files.o ../obj/parser.o ../obj/constants.o ../obj/timer.o ../obj/plot3Dio.o ../obj/chest.o ../obj/util.o ../obj/formats.o ../obj/object.o ../obj/molecule.o ../obj/PDB.o ../obj/erfsun_2_f90.o ../obj/erfappx.o ../obj/ivector_list.o ../obj/ivector.o ../obj/rng.o ../obj/clcg.o ../obj/stats.o ../obj/system.o ../obj/sysinfo.o ../obj/sysmanip.o ../obj/bestfit.o ../obj/ch_param.o ../obj/atompar.o ../obj/anglpar.o ../obj/bondpar.o ../obj/dihepar.o ../obj/psf.o ../obj/psfatom.o ../obj/tlist.o ../obj/mol_formats.o ../obj/charmmio.o ../obj/pdbio.o ../obj/freeio.o ../obj/angles.o ../obj/bonds.o ../obj/dihedrals.o -L ../lib -o ches
make[1]: Leaving directory `/home/taly/projects/chest'
