/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
! **********************************************************************!
! This source file was was generated automatically from a master source !
! code tree, which may not be distributed with this code if the !
! distributor has a proprietary compilation procedure (e.g. CHARMM) !
! If you edit this file (rather than the master source file) !
! your changes will be lost if another pull from the master tree occurs.!
! In case you are wondering why, this approach makes it possible for !
! me to have the same master source code interfaced with different !
! applications (some of which are written in a way that is quite far !
! from being object-oriented) at the source level. !
! **********************************************************************!
      module gridsize ! must be upper-case for the preprocessor to recognize
! contains basic information about the size of the calculation
       private ! keep mpif.h hidden here
!
       public size_initialize
       integer, save, public :: nx=0, ny=0, nz=0
!
       integer*4, save, public :: communicator, ncpu=-1, me=-1
       logical, save, public :: size_initialized, q2D
!
      contains
       subroutine size_initialize(mx, my, mz, comm)
       use output, only: error, message, qprint
       implicit none
       integer :: mx, my, mz
       integer*4, optional :: comm
       integer :: bug
       character(len=8), parameter :: whoami='SET_SIZE'
!
       size_initialized=.false.
       q2D=.false.
!
       me=0
       ncpu=1
       communicator=-1
! note: qprint lives in output module
       if (me.eq.0) then ; qprint=.true. ; else ; qprint=.false.; endif
!
       if (mx.lt.4) then
        call error(whoami, 'NX less than 4. Aborting.',-1)
        return
       else
        nx=mx
       endif
!
       if (my.lt.4) then
        call error(whoami, 'NY less than 4. Aborting.',-1)
        return
       else
        ny=my
       endif
!
       if (mz.lt.4) then
        call message(whoami, 'NZ less than 4. Assuming 2D configuration (NZ = 3).')
        q2D=.true.
        nz=3
       else
        nz=mz
        q2D=.false.
       endif
!
       size_initialized=.true.
!
       end subroutine size_initialize
      end module gridsize
