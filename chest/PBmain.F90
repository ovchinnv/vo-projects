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
      module PBmain
!
! include other modules
      use state ! main physical variables
      use parser ! parses input unit and stores input parameters
      use gridsize ! basic size information
      use grid ! grid module
      use bc ! boundary conditions
      use output ! runtime output
      use object ! objects that modify problem coefficients
! use plot3Dio ! output in plot3D (CFD) format; compile as a library, not as a module
      use multigrid ! multigrid solver module
! use openDXio ! output in openDX format ; not implemented yet
!**************************************************************************************
! routines:
      public PB_init
      public PB_done
      public PB_solve
      public PB_output
!
! variables:
!
      logical :: PB_initialized=.false.
!
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine PB_init( &
     &)
       use files
       implicit none
!
       integer*4 :: numarg
       character(len=80) :: filename, fname, keyword
       integer :: flen
       integer :: fid=-1
       logical :: qcomm
       integer :: bug
!
       character(len=7), parameter :: whoami='PB_INIT'
!
! if (PB_initialized) call PB_done()
!
       qcomm=.false.
! process input file
! *********************************************************************************************************
       if (me.le.0) then
        numarg=command_argument_count() ! number of arguments
        if (numarg.ge.1) then
         call getarg(1, filename)
         fname=filename
         call adjustleft(fname)
         flen=len_trim(fname)
         call files_open(handle_=fid, name_=filename(1:flen), form_='FORMATTED', action_='READ')
         if (fid.lt.0) call error(whoami,'Could not open file',-1)
        else
         fid=5 ! if file name missing, read from standard input
        endif
       endif ! me
! call parser
! only the root will receive a valid file handle (see above); all nodes still call parse_file
       call parse_file(fid) ! parser will store commands internally -- we can now query parser for options using command()
!
       if (me.le.0.and.fid.ne.5) call files_close(fid)
! *********************************************************************************************************
! * INITIALIZATION
! *
! 1) size
       if (qcomm) then







/* nothing */

       else
        if (existtag_nocase('NZ')) then
         call size_initialize(atoi(getval_nocase('NX')), atoi(getval_nocase('NY')), atoi(getval_nocase('NZ')))
        else ! 2D
         call size_initialize(atoi(getval_nocase('NX')), atoi(getval_nocase('NY')), 0)
        endif
       endif
!
! prepare Cartesian grid and metrics
       call message(whoami, 'Initializing grid variables.');
       call grid_initialize()
! initialize state variables
       call message(whoami, 'Initializing state variables.');
! write(0,*) "me:",me
       call state_initialize()
! initialize boundary conditions
       call message(whoami, 'Initializing boundary condition variables.');
       call bc_initialize()
! check whether object specified; if so, initialize
       if (existtag('OBJECT') .or. existtag('object')) then
! check if the object has not yet been initialized (as could happen via grid module above); if not, initialize
        if (.not. object_initialized) then
         call message(whoami, 'Initializing object.');
         call object_init()
        endif
       endif
!
! at present, there are no cases in which a solution does not involve an object
!
! object dependent properties
       if (object_initialized) then
        call message(whoami, 'Reading object parameters.');
        call object_read_parameters()
        if (.not.size_initialized) then
         call warning(whoami, 'Size not initialized.',-1);
        elseif (.not.grid_initialized) then
         call warning(whoami, 'Grid not initialized.',-1);
        elseif (.not.state_initialized) then
         call warning(whoami, 'State variables not initialized.',-1);
        else
         call message(whoami, 'Computing gridded parameters.');
         call object_grid_objects( xcen, ycen, zcen, dxcor, dycor, dzcor, eps, kappa, rhs, nx, ny, nz )
        endif
!
        if (qcomm) then



/* nothing */

        else
         if (fatal_warning()) call terminate(whoami)
        endif
       endif
!
       PB_initialized=.true.
!
       end subroutine PB_init
!************************************************************************************
       subroutine PB_solve()
       implicit none
!
       character(len=20) :: solver
!
! initialize multigrid solver, if requested
!
       character(len=8), parameter :: whoami='PB_SOLVE'
!
! the proper way to terminate run if error occurs on some processors (but not all)
       if (.not. PB_initialized) call warning(whoami, 'POISSON_BOLTZMANN MODULE NOT INITIALIZED.',-1)
       if (fatal_warning()) call terminate(whoami)
!
       if (existtag('SOLVER')) then ; solver=getval('SOLVER'); call toupper(solver);
       elseif (existtag('solver')) then ; solver=getval('solver'); call toupper(solver);
       else ; ! no solver defined
        call warning(whoami, 'No solver specified. Nothing done.',0);
        return
       endif
!
! special case: NONE - skip solution
!
       if (solver.eq.'NONE') return
!
! if still here, attempt to call solver
!
       select case(solver)
!*********************************************************************************!
 case("MULTIGRID");
  if (.not. MULTIGRID_INITIALIZED) call MULTIGRID_INIT();
  call message(whoami, 'CALLING '// &
    "MULTIGRID" //' SOLVER')
  call MULTIGRID_SOLVE();
!*********************************************************************************!
       case default;
        call warning(whoami, 'Unknown solver '//solver(1:len_trim(solver))//' specified. Nothing done.',-1);
       end select
!
       end subroutine PB_solve
!*********************************************************************************!
       subroutine PB_output()
       use object, only : object_surface_pointer, object_initialized
       use formats
       implicit none
!
       character(len=9), parameter :: whoami='PB_OUTPUT'
       integer :: ifmt=plot3d
       logical :: qbin ! whether binary output is used
       real*8, pointer :: surf(:,:,:)
       integer :: i
!
       character(len=80) :: filename, output_format, output_mode
!
       if (existtag_nocase('output_format')) then
        output_format=getval_nocase_upper('output_format');
        ifmt=-999
        do i=1,num_fmt
         if (output_format.eq.format_name(i)) then
          ifmt=i
          exit
         endif
        enddo
       endif
!
       if (ifmt.gt.0) then
        call message(whoami,format_name(ifmt)(1:len_trim(format_name(ifmt)))//' format will be used for output');
       else
        call warning(whoami,'Format "'//output_format(1:len_trim(output_format))//'" is not recognized',0)
        return
       endif
!
       if (existtag_nocase('output_mode')) then
        output_mode=getval_nocase('output_mode')
        select case(output_mode);
         case('BINARY', 'BIN', 'binary', 'bin');
          call message(whoami,'Using binary mode');
          qbin=.true.
         case('ASCII', 'TEXT', 'ascii', 'text');
          call message(whoami,'Using ASCII mode');
          qbin=.false.
         case default
          call warning(whoami,'Unknown output mode "'//output_mode(1:len_trim(output_mode))//'". Will use ASCII',0);
          qbin=.false.
        end select
       else
        call message(whoami,'Using ASCII mode');
        qbin=.false.
       endif
!
       if (existtag_nocase('epsoutput')) then
        filename=getval_nocase('epsoutput')
        call message(whoami,'Writing dielectric to file "'//filename(1:len_trim(filename))//'"');
        select case(ifmt)
         case(plot3d); call plot3Dwrite_scalar(filename,eps,nx,ny,nz,1,qbin)
        end select
       endif
!
       if (existtag_nocase('kappaoutput')) then
        filename=getval_nocase('kappaoutput')
        call message(whoami,'Writing ionic strength to file "'//filename(1:len_trim(filename))//'"');
        select case(ifmt)
         case(plot3d); call plot3Dwrite_scalar(filename,kappa,nx,ny,nz,1,qbin)
        end select
       endif
!
       if (existtag_nocase('chargeoutput')) then
        filename=getval_nocase('chargeoutput')
        call message(whoami,'Writing source density to file "'//filename(1:len_trim(filename))//'"');
        select case(ifmt)
         case(plot3d); call plot3Dwrite_scalar(filename,rhs,nx,ny,nz,1,qbin)
        end select
       endif
!
       if (existtag_nocase('potoutput')) then
        filename=getval_nocase('potoutput')
        call message(whoami,'Writing potential to file "'//filename(1:len_trim(filename))//'"');
        select case(ifmt)
         case(plot3d); call plot3Dwrite_scalar(filename,p,nx,ny,nz,1,qbin)
        end select
       endif
!
! check for surface output options
!
       if (existtag_nocase('surfoutput').or.existtag_nocase('alloutput')) then
        nullify(surf)
        if (object_initialized) then
         surf=>object_surface_pointer()
         if (.not.associated(surf)) then
          call warning(whoami, 'Object initialized but shape density is null. Will output zero.',0)
          allocate(surf(nx,ny,nz)); surf=0d0
         endif
        else
         call warning(whoami, 'Object not initialized. Shape density is everywhere zero.',0)
         allocate(surf(nx,ny,nz)); surf=0d0
        endif
       endif
!
       if (existtag_nocase('surfoutput')) then
        filename=getval_nocase('surfoutput')
        call message(whoami,'Writing shape density to file "'//filename(1:len_trim(filename))//'"');
        select case(ifmt)
         case(plot3d); call plot3Dwrite_scalar(filename,surf,nx,ny,nz,1,qbin)
        end select
       endif
!
! plot five (four in 2D) variables that correspond to the plot3D "solution"
       if (existtag_nocase('alloutput')) then
        filename=getval_nocase('alloutput')
        call message(whoami,'Writing data and solution [EPS(1),K(2),RHS(3),P(4),RHO(5)] to file "'&
                                                        //filename(1:len_trim(filename))//'"')
        select case(ifmt)
         case(plot3d);
          call plot3Dwrite_solution(filename,eps,kappa,rhs,p,surf,nx,ny,nz,1,qbin)
        end select
       endif
!
! destroy surf, if needed
       if (existtag_nocase('surfoutput').or.existtag_nocase('alloutput')) then
        if (.not.object_initialized) then ; deallocate(surf)
        elseif (.not.associated(object_surface_pointer())) then ; deallocate(surf)
        endif
       endif
!
       if (existtag_nocase('gridoutput')) then
        filename=getval_nocase('gridoutput')
        call message(whoami,'Writing grid to file "'//filename(1:len_trim(filename))//'"');
        select case(ifmt)
         case(plot3d); call plot3Dwrite_grid(filename,xcen,ycen,zcen,nx,ny,nz,1,qbin)
        end select
       endif
!
       end subroutine PB_output
!*********************************************************************************!
       subroutine PB_done()
       use files, only: files_done
       implicit none
!
       character(len=7), parameter :: whoami='PB_DONE'
!
       character(len=20) :: solver
!
       if (existtag('OBJECT') .or. existtag('object')) then
        if (object_initialized) then
         call object_done()
        endif
       endif
!
       call state_done()
       call grid_done()
       call bc_done()
!
       if (existtag('SOLVER')) then ; solver=getval('SOLVER'); call toupper(solver);
       elseif (existtag('solver')) then ; solver=getval('solver'); call toupper(solver);
       else
        PB_initialized=.false.
        return
       endif
! special case: NONE -- do nothing
       if (solver.eq.'NONE') then
        PB_initialized=.false.
        return
       endif
!
       select case(solver)
!*********************************************************************************!
 case("MULTIGRID");
  if (MULTIGRID_INITIALIZED) call MULTIGRID_DONE();
!*********************************************************************************!
       end select
!
       call files_done() ! close all units
       call parser_done()
!
       PB_initialized=.false.
!
       end subroutine PB_done
!*********************************************************************************!
       end module PBmain
