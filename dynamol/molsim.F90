/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/



/*
#ifdef __PRINT
#undef __PRINT
#endif
#define __PRINT(__WHAT) call plainmessage(__WHAT)
*/
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
module molsim
! contains setup and details of molecular simulation
! I view this module as a "wrapper" around other modules
! The problem, of course, is that no physics problem is
! truly fully object-oriented, and this is reflected in the
! design/code
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 use parser
 use output
 use system
 use stats
 use verlet
 use rng
 use files
 implicit none
 public molsim_init
 public molsim_integrate
 integer, save :: num_iterations
 integer, save :: iteration_count ! number of MD iterations (persistent)
 integer, parameter :: nfreq=5
 integer, parameter :: printf=1, outf=2, trajf=3, resf=4, statf=5
 integer :: outfreq(nfreq) ! frequency array for writing output
!
! output options persistent in case use runs sequential simulations but wants to kep same files
 character(len=200), save :: trajectoryoutname, restartoutname, statisticsoutname
 integer, save :: trajectoryfid=-1, statisticsfid=-1
!
 logical, save :: restart, verbose=.false.
 logical, save :: molsim_initialized=.false.
!
 contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine molsim_init()
 character(len=80) :: keyword
 integer :: i
!
 character(len=10), parameter :: whoami='MOLSIM_INIT'
 character(len=100) :: parmfilename, structfilename, coorfilename, velfilename
!
 if (existtag_nocase('verbose')) then;verbose=atol(getval_nocase('verbose'));else;verbose=.false.;endif
!%%%%%%%%%%%%%%%%%%% read parameter file(s)%%%%%%%%%%%%%%%%%%%%%%%%
 parmfilename=getval('parameters') ! parameter file(s)
 call system_read_parameters(parmfilename)
 i=2
 do
  write(keyword,*) i
  call adjustleft(keyword)
  if (existtag('parameters'//keyword(1:len_trim(keyword)))) then
   parmfilename=getval('parameters'//keyword(1:len_trim(keyword)))
   call system_read_parameters(parmfilename)
   i=i+1
  else
   exit
  endif
 enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if (verbose) call system_list_parameters()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 structfilename=getval_nocase('structure') ! structure file
 coorfilename=getval_nocase('coordinates') ! coordinate file
!
 call system_read_structure(structfilename)
 call system_read_coordinates(coorfilename)
!
 call system_check() ! check that all parameters are known and whether coordinates/velocities are defined
 call random_init() ! initialize random number generator
 call verlet_init() ! initialize verlet integrator
!
! note that there is no mechanism for restarts -- to be added !
!
 if (existtag_nocase('velocities')) then
  velfilename=getval_nocase('velocities') ! velocity file
  call system_read_velocities(velfilename)
 else
  call system_init_velocities()
 endif
!
 if (verbose) then
  call system_compute() ! if langevin on, already computed; ignore for now
  call system_printe()
 endif
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% request a few parameters from parser %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
if (.not.existtag_nocase("printfreq")) then
  else
   keyword=getval("printfreq");
   call message(whoami, 'Setting '// &
& "Print frequency"//' to '//&
& keyword(1:len_trim(keyword)));
   outfreq(printf) = atoi(keyword)
   if ( outfreq(printf) .lt. 0) &
& call warning(whoami,"Print frequency"//' cannot be negative ('//&
& keyword(1:len_trim(keyword))//'). Abort.',-1)
endif
!#undef
!#undef __DEFAULT
!
if (.not.existtag_nocase("outputfreq")) then
  else
   keyword=getval("outputfreq");
   call message(whoami, 'Setting '// &
& "Output frequency"//' to '//&
& keyword(1:len_trim(keyword)));
   outfreq(outf) = atoi(keyword)
   if ( outfreq(outf) .lt. 0) &
& call warning(whoami,"Output frequency"//' cannot be negative ('//&
& keyword(1:len_trim(keyword))//'). Abort.',-1)
endif
!#undef
!#undef __DEFAULT
!
if (.not.existtag_nocase("trajectoryfreq")) then
  else
   keyword=getval("trajectoryfreq");
   call message(whoami, 'Setting '// &
& "Trajectory output frequency"//' to '//&
& keyword(1:len_trim(keyword)));
   outfreq(trajf) = atoi(keyword)
   if ( outfreq(trajf) .lt. 0) &
& call warning(whoami,"Trajectory output frequency"//' cannot be negative ('//&
& keyword(1:len_trim(keyword))//'). Abort.',-1)
endif
!#undef
!#undef __DEFAULT
!
if (.not.existtag_nocase("restartfreq")) then
  else
   keyword=getval("restartfreq");
   call message(whoami, 'Setting '// &
& "Restart file output frequency"//' to '//&
& keyword(1:len_trim(keyword)));
   outfreq(resf) = atoi(keyword)
   if ( outfreq(resf) .lt. 0) &
& call warning(whoami,"Restart file output frequency"//' cannot be negative ('//&
& keyword(1:len_trim(keyword))//'). Abort.',-1)
endif
!#undef
!#undef __DEFAULT
!
if (.not.existtag_nocase("statisticsfreq")) then
  else
   keyword=getval("statisticsfreq");
   call message(whoami, 'Setting '// &
& "Statistics output frequency"//' to '//&
& keyword(1:len_trim(keyword)));
   outfreq(statf) = atoi(keyword)
   if ( outfreq(statf) .lt. 0) &
& call warning(whoami,"Statistics output frequency"//' cannot be negative ('//&
& keyword(1:len_trim(keyword))//'). Abort.',-1)
endif
!#undef
!#undef __DEFAULT
!
 if (outfreq(trajf).gt.0) then
  if (existtag_nocase('trajectoryfile')) then
   trajectoryoutname=getval_nocase('trajectoryfile')
  else
   call warning(whoami,'Trajectory output file name not specified.',-1)
  endif
 endif
!
 if (outfreq(resf).gt.0) then
  if (existtag_nocase('restartfile')) then
   restartoutname=getval_nocase('restartfile')
  else
   call warning(whoami,'Restart file name not specified.',-1)
  endif
 endif
!
 if (outfreq(statf).gt.0) then
  if (existtag_nocase('statisticsfile')) then
   statisticsoutname=getval_nocase('statisticsfile')
  else
   call warning(whoami,'Statistics output file name not specified.',-1)
  endif
 endif
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if (fatal_warning()) then
  call system_done()
 else
!
  iteration_count=0
  molsim_initialized=.true.
!
 endif
!
 end subroutine molsim_init
!
!%%%%%%%%%%%%%%%%%%%%%%% integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine molsim_integrate(numiter)
 integer, optional :: numiter ! optional number of integraion steps
 integer :: minfreq, iterations, ncycles, i
 logical :: addheader
 character(len=6) :: action
!
 character(len=16), parameter :: whoami='MOLSIM_INTEGRATE'
 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_
 integer*4 :: me
 me=0
!
 if (present(numiter)) then ; num_iterations=numiter ;
 else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% get number of iterations from input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
if (.not.existtag_nocase("iterations")) then
   call warning(whoami, "Number of simulation steps"//' unspecified.',-1)
  else
   msg___ (1)=getval("iterations");
   call message(whoami, 'Setting '// &
& "Number of simulation steps"//' to '//&
& msg___ (1)(1:len_trim(msg___ (1))));
   num_iterations = atoi(msg___ (1))
   if ( num_iterations .lt. 0) &
& call warning(whoami,"Number of simulation steps"//' cannot be negative ('//&
& msg___ (1)(1:len_trim(msg___ (1)))//'). Abort.',-1)
endif
!#undef
!#undef __DEFAULT
!
 endif ! numiter
!%%%%%%%%%%%%%%%%%%%%%%%%%%% get number of iterations from input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
 if (any(outfreq>0)) then ; minfreq=minval(outfreq,1,outfreq>0); else; minfreq=num_iterations; endif ! minimum frequency
 if (minfreq.ne.0) then
  if (sum(mod(outfreq,minfreq)).gt.0) &
 & call warning(whoami, 'One of the frequency values must divide the others.',-1)
 else
  minfreq=1 ! in this case, niter nust be zero (see above), so ncycle will be zero below
 endif
 ncycles=num_iterations/minfreq+min(mod(num_iterations,minfreq),1)
!
 do i=1, ncycles
  iterations=min(minfreq,num_iterations-iteration_count)
  call verlet_integrate(iterations)
  iteration_count=iteration_count+iterations
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (outfreq(printf).gt.0) then;if (mod(iteration_count,outfreq(printf)).eq.0) call system_printe(); endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (outfreq(resf).gt.0) then;if (mod(iteration_count,outfreq(resf)).eq.0) call system_write_restart(restartoutname); endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (outfreq(statf).gt.0) then;
   if (mod(iteration_count,outfreq(statf)).eq.0) then
    if (me.le.0) then
     if (statisticsfid.le.0) then ; action='write' ; else ; action='append'; endif
     call files_open(statisticsfid, name_=statisticsoutname, form_='FORMATTED', action_=action)
     if (statisticsfid.le.0) call warning(whoami, 'Cannot open input file. Abort.',-1)
    endif ! me
!
! call system_statistics(statisticsfid); ! must be parallel-aware; must call fatal_warning
!
    call files_close(statisticsfid)
!
   endif ! mod
  endif ! outfreq
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (outfreq(trajf).gt.0) then;
   if (mod(iteration_count,outfreq(trajf)).eq.0) then
    if (me.le.0) then
!
! note: in the future should close and reopen traj. file with APPEND option to protect from crashes
!
     if (trajectoryfid.le.0) then ; addheader=.true.; action='write'; else ; addheader=.false.; action='append'; endif
     call files_open(trajectoryfid, name_=trajectoryoutname, form_='UNFORMATTED', action_=action)
     if (trajectoryfid.le.0) call warning(whoami, 'Cannot open input file. Abort.',-1)
    endif ! me
!
    call system_write_dcd(trajectoryfid,addheader); ! must be parallel-aware; must call fatal_warning
    call files_close(trajectoryfid)
!
   endif ! mod
  endif ! outfreq
!
  if (fatal_warning()) exit
!
 enddo ! ncycles
! final restart file
 if (outfreq(resf).gt.0) call system_write_restart(restartoutname)
!optional coordinate/and/or/velocity files
 if (existtag_nocase('finalcoordinates')) call system_write_coordinates(getval_nocase('finalcoordinates'))
 if (existtag_nocase('finalvelocities')) call system_write_velocities(getval_nocase('finalvelocities'))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end subroutine molsim_integrate
end module molsim
