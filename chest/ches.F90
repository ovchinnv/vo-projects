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
  PROGRAM CHES ! (Cartesian Helmholtz Equation Solver)
! 2010-2012 Victor Ovchinnikov (ovchinnv_at_georgetown_dot_edu)
! Free for academic use
! Most of the source code is distributed
! under the GNU General Public License
! certain parts of the code are under restricted license
!
  use PBmain
  use timer
  use parser
  use output
  implicit none
!
  character(len=4) :: whoami='CHES'
  integer*4 :: comm ! default communicator
  integer :: pbtimer
  real*8 :: time(10)
!
  comm=-1
!
  call timer_init()
  pbtimer=timer_start()
!

  call PB_init()


 time(1)=timer_stamp(pbtimer)

! solve the Poisson-Boltzmann equation (if requested in the input file)
  call PB_solve()

  time(2)=timer_stamp(pbtimer)

! output solution, or other properties (as requested in the input file)
  call PB_output()

  time(3)=timer_stamp(pbtimer)

! deallocate arrays
  call PB_done()
!



!

  call timer_done()

!

 call message(whoami,'============= TIMING INFORMATION ==============')
 call message(whoami,'PB initialization time (s) : '//ftoa(time(1)))
 call message(whoami,'PB solution time (s)       : '//ftoa(time(2)))
 call message(whoami,'PB output time (s)         : '//ftoa(time(3)))
 call message(whoami,'PB finalization time (s)   : '//ftoa(time(4)))
 call message(whoami,'===============================================')

  end ! CHES
