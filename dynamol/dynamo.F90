program dynamo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 use parser
 use output
 use system
 use stats
 use verlet
 use rng
 
 implicit none

 int :: niter            ! number of iterations (steps)
 int :: outfreq(2)       ! frequency array for writing output
 int, parameter :: prn=1, out=2
 int :: minfreq          ! minimum frequency
!
 int :: iteration_counter     ! number of MD iterations
!
 int :: iseed(4)=(/1,2,3,4/) ! random number seeds
 int :: channel = 1           ! channel from which to draw random numbers
!
 integer*4 :: numarg 
 character*(80) :: filename, fname, keyword
 int :: fid=1, flen
!
 character*6, parameter :: whoami='DYNAMO'
!
 character*100 :: parmfilename, structfilename, coorfilename, velfilename

 bool :: restart
 int :: i, iteration, ncycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 numarg=iargc() ! number of arguments
 if (numarg.ge.1) then 
    call getarg(1, filename)
    fname=filename
    call adjustleft(fname)
    flen=len_trim(fname)
    open(fid, file=fname(1:flen), status='OLD', form='FORMATTED')
 else
    fid=5 ! if file name missing, read from standard input
 endif
!    call parser
 call parse_file(fid) ! parser will store commands internally -- we can now query parser for options using command()
!
 close(fid)
!
!%%%%%%%%%%%%%%%%%%% read parameter file(s)%%%%%%%%%%%%%%%%%%%%%%%%
 parmfilename=getval('parameters')       ! parameter file(s)
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
! call system_list_parameters()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 structfilename=getval('structure')     ! structure file
 coorfilename=getval('coordinates')     ! coordinate file
 if (existtag('velocities')) then
  restart=.true.
  velfilename=getval('velocities')      ! velocity file
 else
  restart=.false.
 endif

 call system_read_structure(structfilename)
 call system_read_coordinates(coorfilename)

 call system_check() ! check that all parameters are known and whether coordinates/velocities are defined 
! call system_compute() ! compute energy and forces

!% initialize random number generator
 call random_init(iseed)
!
 if (restart) then 
  call system_read_velocities(velfilename)
 else 
  call system_init_velocities()
 endif
!
 call system_printe()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! determine integration frequency
 niter=atoi(getval('iterations'))
 if (niter.lt.0) call error(whoami, 'NUMBER OF SIMULATION STEPS NEGATIVE. ABORT',-1)
!
 if (existtag('printfreq')) then; outfreq(1)=atoi(getval('printfreq')); else; outfreq(1)=0; endif
 if (outfreq(1).lt.0) call error(whoami, 'PRINT FREQUENCY NEGATIVE. ABORT',-1)
!
 if (existtag('outputfreq')) then; outfreq(2)=atoi(getval('outputfreq')); else; outfreq(2)=0; endif
 if (outfreq(2).lt.0) call error(whoami, 'OUTPUT FREQUENCY NEGATIVE. ABORT',-1)
!
!   
 if (any(outfreq>0)) then ; minfreq=minval(outfreq,1,outfreq>0); else; minfreq=niter; endif ! minimum frequency
 if (minfreq.ne.0) then
  if (sum(mod(outfreq,minfreq)).gt.0) call error(whoami, 'ONE OF THE FREQUENCY VALUES MUST DIVIDE ALL THE OTHERS. ABORT.',-1)
 else
  minfreq=1 ! in this case, niter nust be zero (see above), so ncycle will be zero below
 endif
 ncycle=niter/minfreq+min(mod(niter,minfreq),1)

!%%%%%%%%%%%%%%%%%%%%%%% integrate %%%%%%%%%%%%%%%%%%%%%%% 
 iteration=0
 do i=1, ncycle
  call verlet_integrate(min(minfreq,niter-iteration))
  iteration=iteration+minfreq
! 
  if (outfreq(prn).gt.0) then; if  (mod(iteration,outfreq(prn)).eq.0) call system_printe(); endif
!
 enddo ! i
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! call system_write_coordinates()
 
end program dynamo



