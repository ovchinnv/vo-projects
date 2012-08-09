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
!DEC$ ATTRIBUTES INLINE :: residual
!DEC$ ATTRIBUTES INLINE :: refine
!DEC$ ATTRIBUTES INLINE :: apply_bc
!DEC$ ATTRIBUTES INLINE :: residual2d
!DEC$ ATTRIBUTES INLINE :: fzero
!DEC$ ATTRIBUTES INLINE :: Jacobi
!DEC$ ATTRIBUTES INLINE :: JacobiTiledLoMem
!DEC$ ATTRIBUTES INLINE :: GaussSeidel
!DEC$ ATTRIBUTES INLINE :: GaussSeidelRB
!DEC$ ATTRIBUTES INLINE :: GaussSeidelRBTiled
!DEC$ ATTRIBUTES INLINE :: GaussSeidelRBTiledLoMem
!DEC$ ATTRIBUTES INLINE :: GSInner
!DEC$ ATTRIBUTES INLINE :: GSOuter
!DEC$ ATTRIBUTES INLINE :: GaussSeidel2d
!DEC$ ATTRIBUTES INLINE :: copy3d
!DEC$ ATTRIBUTES INLINE :: coarsen
!DEC$ ATTRIBUTES INLINE :: coarsen2d_vec
!DEC$ ATTRIBUTES INLINE :: compute_fd_coef
!DEC$ ATTRIBUTES INLINE :: compute_fd_coef2d
 module multigrid
!
 use constants
 use timer
 private
 logical, public :: multigrid_initialized=.false.
 integer, parameter :: nxmin=3, nymin=3, nzmin=3 ! minimum allowed number of inner points
 integer, save :: maxlev ! maximum number of allowed levels (user defined)
 integer, parameter :: default_maxlev=1000
!
 real*8, pointer, save, dimension(:) :: dxall, dyall, dzall, odxall, odyall, odzall, & ! arrays that hold all metrics
                                       allp, alleps, allkappa, allrhs,& ! arrays that hold coefficients, sources, and solution
                                       alleast, allwest, allsouth, allnorth, allfront, allback ! arrays that store boundary conditions
!
 real*8, pointer, save, dimension(:,:) :: bc_wgt
!
 integer, save :: dxlen, dylen, dzlen, len3D, len3Dbc, len2Dxy, len2Dxz, len2Dyz ! corresponding 1D array lengths
 integer, save :: numlev=-1 ! number of multigrid levels
 integer, pointer, save :: mcycle(:) ! this stores cycle definition (smooth, prolong, or restrict)
 integer, pointer, save :: vcycle(:) ! default vcycle
 character(len=200), save :: multigrid_cycle_spec ! description of multigrid cycle
!
 integer, parameter :: reduce=-1, prolong=-2, done=0; ! this is just notation; currently positive numbers ar assumed to be iterations
!
 integer, parameter :: nsmoother = 9
 character(len=24) :: smoother_names(nsmoother) = &
& (/'JACOBI                  ','GAUSS-SEIDEL            ','GAUSS-SEIDEL-RB         ',&
& 'GAUSS-SEIDEL-UNROLL     ','GAUSS-SEIDEL-RB-TILED   ','GAUSS-SEIDEL-RB-TILED-LM',&
& 'GAUSS-SEIDEL-REVERSE    ','JACOBI-TILED-LM         ','JACOBI-LM               '/)
 integer, parameter :: Jacobi_=1, GaussSeidel_=2, GaussSeidelRB_=3, GSUnrollFromMiddle_=4, GaussSeidelRBTiled_=5,&
& GaussSeidelRBTiledLoMem_=6, GaussSeidelReverse_=7, JacobiTiledLomem_=8, JacobiOnTheFly_=9
 integer, save :: smoother, default_smoother=GaussSeidel_
 integer, save :: smooth_iterations;
 integer, save :: smooth_tilesize;
 integer, parameter :: default_unroll_iter=2 ; ! default inner iterations for "unrolled" GS
 integer, save :: unroll_iter; ! inner iterations for "unrolled" GS
 integer, parameter :: default_smooth_iterations=10;
 integer, parameter :: default_smooth_tilesize=256;
 integer, save :: maxcycle ! maximum number of cycles
 integer, parameter :: default_maxcycle=1000 ! maximum number of cycles; after this is reached, issue a warning and quit
 integer, save :: convergence_skip_cycles;
 integer, parameter :: default_convergence_skip_cycles=0;
 real*8, save :: max_residual ! maximum residual allowed
 real*8, save :: default_max_residual ! initialized below
 real*8, save :: init_residual
 real*8, save :: omega ! overrelaxation parameter in GS
 real*8, parameter :: default_omega=1
 logical, save :: multigrid_loud=.true.
 logical, save :: compute_initial_residual
 logical, parameter :: default_compute_initial_residual=.true.
 integer(kind=KIND('a')), save :: i3D=1, i2D=0
!
 integer :: mgtimer, smoothtimer
 real*8 :: coarsen_time=zero, refine_time=zero, smooth_time=zero, residual_time=zero, bc_time=zero
 real*8 :: mgtime=-1d0
!
 public multigrid_init
 public multigrid_done
 public multigrid_solve
! the next three routine updates coefficients (eps, kappa, and rhs)
 public multigrid_update_coef
 public multigrid_update_bc
 private parse_multigrid_cycle
!
 contains
!*******************************************************************************!
  subroutine multigrid_init()
  use gridsize, only: nx, ny, nz, me, communicator, q2D
  use grid
  use output
  use parser
  use constants
  implicit none
!
  character(len=14), parameter :: whoami = 'MULTIGRID_INIT'
  character(len=100) :: msg
  character(len=100) :: keyword
  integer :: msglen
  integer :: nnx, nny, nnz, l, len_cycle
  integer :: dxi, dxj, dxci, dxcj, dyi, dyj, dyci, dycj, dzi, dzj, dzci, dzcj
  real*8, allocatable, dimension(:) :: x0,x1,y0,y1,z0,z1,dx0,dy0,dz0
  real*8 :: fac1d, fac2d, fac3d
  integer :: i, ii, j
!
  if (multigrid_initialized) call multigrid_done()
!
  if (.not.grid_initialized) then
    call warning(whoami, 'Main grid not initialized. Nothing done.', 0)
    return
  endif
!
  if (q2d) then ; i3d=0 ; call message(whoami, '2D configuration specified'); else ; i3d=1 ; endif ; i2d=1-i3d;
!
  default_max_residual=1.d0*ERRTOL()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% request a few parameters from parser %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
if (.not.existtag_nocase("multigrid_maximum_levels")) then
   call warning(whoami, "maximum number of grid levels"//' unspecified. Defaulting to '// &
& itoa( default_maxlev ),0)
   maxlev=default_maxlev
  else
   keyword=getval("multigrid_maximum_levels");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "maximum number of grid levels"//' to '//keyword(1:l));
   maxlev = atoi(keyword)
   if ( maxlev .lt. 0) &
& call warning(whoami,"maximum number of grid levels"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% smoother selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (.not.existtag_nocase('MULTIGRID_SMOOTHER')) then
   l=len_trim(smoother_names(default_smoother))
   call warning(whoami,' Multigrid smoother unspecified. Defaulting to "'//smoother_names(default_smoother)(1:l)//'"',0)
   smoother=default_smoother
  else
   keyword=getval_nocase_upper('MULTIGRID_SMOOTHER')
   l=len_trim(keyword)
   smoother=-999
   do i=1,nsmoother
     if (keyword.eq.smoother_names(i)) then
      smoother=i
      exit
     endif
   enddo
   if (smoother.gt.0) then
    call message(whoami, 'Setting multigrid smoother to "'//smoother_names(smoother)(1:l)//'"')
   else
    call warning(whoami, 'Unknown multigrid smoother specified. Abort.',-1)
   endif
  endif
!
  if (smoother.eq.GSUnrollFromMiddle_) then
!
if (.not.existtag_nocase("multigrid_smoother_unroll")) then
   call warning(whoami, "maximum number of unrolled smoother loops"//' unspecified. Defaulting to '// &
& itoa( default_unroll_iter ),0)
   unroll_iter=default_unroll_iter
  else
   keyword=getval("multigrid_smoother_unroll");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "maximum number of unrolled smoother loops"//' to '//keyword(1:l));
   unroll_iter = atoi(keyword)
endif
!#undef __MINUSERR
!#undef
   if (unroll_iter.lt.0) &
   & call message(whoami, 'Will perform maximal unrolling');
!
  elseif ( (smoother.eq.GaussSeidelRBTiled_).or.(smoother.eq.GaussSeidelRBTiledLoMem_).or.(smoother.eq.JacobiTiledLoMem_) ) then
if (.not.existtag_nocase("multigrid_smoother_tilesize")) then
   call warning(whoami, "smoother tile size"//' unspecified. Defaulting to '// &
& itoa( default_smooth_tilesize ),0)
   smooth_tilesize=default_smooth_tilesize
  else
   keyword=getval("multigrid_smoother_tilesize");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "smoother tile size"//' to '//keyword(1:l));
   smooth_tilesize = atoi(keyword)
   if ( smooth_tilesize .lt. 0) &
& call warning(whoami,"smoother tile size"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef
  endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
if (.not.existtag_nocase("multigrid_smooth_iterations")) then
   call warning(whoami, "number of smoother iterations"//' unspecified. Defaulting to '// &
& itoa( default_smooth_iterations ),0)
   smooth_iterations=default_smooth_iterations
  else
   keyword=getval("multigrid_smooth_iterations");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "number of smoother iterations"//' to '//keyword(1:l));
   smooth_iterations = atoi(keyword)
   if ( smooth_iterations .lt. 0) &
& call warning(whoami,"number of smoother iterations"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef
!
if (.not.existtag_nocase("multigrid_maximum_cycles")) then
   call warning(whoami, "maximum number of cycles"//' unspecified. Defaulting to '// &
& itoa( default_maxcycle ),0)
   maxcycle=default_maxcycle
  else
   keyword=getval("multigrid_maximum_cycles");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "maximum number of cycles"//' to '//keyword(1:l));
   maxcycle = atoi(keyword)
   if ( maxcycle .lt. 0) &
& call warning(whoami,"maximum number of cycles"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef
!
if (.not.existtag_nocase("multigrid_skip_convergence_test")) then
   call warning(whoami, "Number of cycles before first convergence test"//' unspecified. Defaulting to '// &
& itoa( default_convergence_skip_cycles ),0)
   convergence_skip_cycles=default_convergence_skip_cycles
  else
   keyword=getval("multigrid_skip_convergence_test");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "Number of cycles before first convergence test"//' to '//keyword(1:l));
   convergence_skip_cycles = atoi(keyword)
   if ( convergence_skip_cycles .lt. 0) &
& call warning(whoami,"Number of cycles before first convergence test"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef
!
if (.not.existtag_nocase("multigrid_compute_ini_residual")) then
   call warning(whoami, "compute-initial-residual"//' unspecified. Defaulting to '// &
& ltoa( default_compute_initial_residual ),0)
   compute_initial_residual=default_compute_initial_residual
  else
   keyword=getval("multigrid_compute_ini_residual");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "compute-initial-residual"//' to '//keyword(1:l));
   compute_initial_residual = atol(keyword)
endif
!#undef __MINUSERR
!#undef
!
if (.not.existtag_nocase("multigrid_residual")) then
   call warning(whoami, "maximum residual"//' unspecified. Defaulting to '// &
& ftoa( default_max_residual ),0)
   max_residual=default_max_residual
  else
   keyword=getval("multigrid_residual");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "maximum residual"//' to '//keyword(1:l));
   max_residual = atof(keyword)
   if ( max_residual .lt. 0) &
& call warning(whoami,"maximum residual"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef
!
if (.not.existtag_nocase("multigrid_omega")) then
   call warning(whoami, "over-relaxation parameter"//' unspecified. Defaulting to '// &
& ftoa( default_omega ),0)
   omega=default_omega
  else
   keyword=getval("multigrid_omega");
   l=len_trim(keyword)
   call message(whoami, 'Setting '// &
& "over-relaxation parameter"//' to '//keyword(1:l));
   omega = atof(keyword)
   if ( omega .lt. 0) &
& call warning(whoami,"over-relaxation parameter"//' cannot be negative ('//keyword(1:l)//'). Abort.',-1)
endif
!#undef
!#undef
!
!**********************************************
! determine the maximum number of levels possible (3D/2D)
  numlev=1;
  nnx=nx-2; nny=ny-2; nnz=nz-2; ! inner points
  do while ( (mod(nnx,2).eq.0).and.(mod(nny,2).eq.0).and.(mod(nnz,2)*i3d.eq.0).and. &
& nnx.ge.2*nxmin.and.nny.ge.2*nymin.and.nnz.ge.2*nzmin*i3d.and. &
& numlev.lt.maxlev )
   numlev=numlev+1
   nnx=nnx/2; nny=nny/2; nnz=nnz/(2-i2d);
  enddo
!
  call message(whoami, 'NUMBER OF MULTIGRID LEVELS IS '//itoa(numlev))
!
! define default v-cycle
!***********************************************
  allocate(vcycle(4*(numlev-1)+1)); ! smooth & coarsen numlev-1 times; smooth; interpolate & smooth numlev-1 times
  ii=0
  do i=1,numlev-1
   ii=ii+1 ; vcycle(ii)=smooth_iterations
   ii=ii+1 ; vcycle(ii)=reduce
  enddo
  ii=ii+1 ; vcycle(ii)=smooth_iterations
  do i=1,numlev-1
   ii=ii+1 ; vcycle(ii)=prolong
   ii=ii+1 ; vcycle(ii)=smooth_iterations
  enddo
!****************************************************************************
! parse cycle definition; default will assume a vcycle with iterations specified above
!****************************************************************************
  nullify(mcycle)
!
  if (existtag_nocase('multigrid_cycle')) then
   multigrid_cycle_spec=getval_nocase('multigrid_cycle')
   call toupper(multigrid_cycle_spec)
   select case(multigrid_cycle_spec)
    case('VCYCLE', 'V', 'V-CYCLE')
     mcycle=>vcycle
    case default ! process a custom cycle definition
     call parse_multigrid_cycle(multigrid_cycle_spec, mcycle)
     if (.not.associated(mcycle)) then
      call warning(whoami, 'INVALID CYCLE SPECIFIED. ABORT.',-1)
     endif
   end select
  else
   mcycle=>vcycle
  endif
!**************************************************************************
! check that the multigrid cycle makes sense:
!
  l=1; ! start at the first level
  len_cycle=size(mcycle)
  do i=1,len_cycle
   j=mcycle(i)
   if (j.eq.done) then
! comment out the next two lines to parse until the end
    len_cycle=i-1
    exit
   elseif (j.eq.reduce) then
    l=l+1;
    if (l.gt.numlev) then ! invalid level
     call warning(whoami, 'INVALID CYCLE SPECIFIED. ABORT.',-1);
     exit
    endif
   elseif (j.eq.prolong) then
    l=l-1;
    if (l.lt.1) then ! invalid level (allowed levels are [ 0 ... numlev-1 ]
     call warning(whoami, 'INVALID CYCLE SPECIFIED. ABORT.',-1);
     exit
    endif
   elseif (j.lt.0) then ! this means an invalid entry
    call warning(whoami, 'INVALID CYCLE SPECIFIED. ABORT.',-1);
    exit
   endif ! all positive entries are valid and correspond to smoother iterations
  enddo
!
  if ( (l.ne.1) .or. (len_cycle.eq.0) ) call warning(whoami, 'INVALID CYCLE SPECIFIED. ABORT.',-1);
  if (fatal_warning()) call terminate(whoami)
!**************************************************************************
!
! create metrics for all levels
  nnx=nx-2; nny=ny-2; nnz=nz-2-i2d; ! inner points at first level; in 2D, set nnz=0 here
!
! allocate storage space for additional grids
!
  fac1d=(1d0-0.5d0**numlev)/(1d0-0.5d0);
!
  dxlen= numlev+nint(2d0*nnx*fac1d); ! two sets of metrics; (one of two metrics has an "extra" element)
  dylen= numlev+nint(2d0*nny*fac1d);
  dzlen= numlev+nint(2d0*nnz*fac1d); ! note: dzall, odzall have to have nonzero length in 2D even though they are not used
!
  allocate(dxall(dxlen),odxall(dxlen),dyall(dylen),odyall(dylen),dzall(dzlen),odzall(dzlen))
!
! metrics at first level
  dxi=1;
  dxj=dxi+nnx; dxci=dxj+1; dxcj=dxci+nnx-1;
  dyi=1;
  dyj=dyi+nny; dyci=dyj+1; dycj=dyci+nny-1;
  dzi=1;
  dzj=dzi+nnz-i2d; dzci=dzj+1; dzcj=dzci+nnz-1;
!
  dxall(dxi:dxj) = dxcen(1:nnx+1);
  dxall(dxci:dxcj)= dxcor(2:nnx+1);
  dyall(dyi:dyj) = dycen(1:nny+1);
  dyall(dyci:dycj)= dycor(2:nny+1);
  dzall(dzi:dzj) = dzcen(1:nnz+i3d);
  dzall(dzci:dzcj)= dzcor(2:nnz+1);
!
! arrays for gridding at other levels
!
  allocate(x0(nnx),x1(nnx),y0(nny),y1(nny),z0(nnz),z1(nnz),dx0(nnx),dy0(nny), dz0(nnz))
!
  x0(1:nnx)=xcen(2:nx-1); ! take internal points only
  y0(1:nny)=ycen(2:ny-1);
  z0(1:nnz)=zcen(2:nz-1-i2d); ! OK in 3D/2D
!
  do l=2,numlev
! create coarser grid:
   x1(1:nnx/2)=0.5d0*(x0(1:nnx:2)+x0(2:nnx:2)); ! also contains only internal points (on a coarse grid)
   y1(1:nny/2)=0.5d0*(y0(1:nny:2)+y0(2:nny:2));
   z1(1:nnz/2)=0.5d0*(z0(1:nnz:2)+z0(2:nnz:2)); ! OK in 3D/2D
! compute coarse-level metrics
   nnx=nnx/2; nny=nny/2; nnz=nnz/2;
! center-to-center metrics (reuse coordinate arrays)
   x0(2:nnx)=x1(2:nnx)-x1(1:nnx-1); x0(1)=x0(2); x0(nnx+1)=x0(nnx); ! first metric would contain 1st boundary point
   y0(2:nny)=y1(2:nny)-y1(1:nny-1); y0(1)=y0(2); y0(nny+1)=y0(nny);
   z0(2:nnz)=z1(2:nnz)-z1(1:nnz-1); if (.not.q2d) then ; z0(1)=z0(2); z0(nnz+1)=z0(nnz); endif;
! corner-to-corner metrics
   dx0(1:nnx)=0.5d0*(x0(1:nnx)+x0(2:nnx+1)); ! first metric contains 1st point (no ghost pt.)
   dy0(1:nny)=0.5d0*(y0(1:nny)+y0(2:nny+1));
   dz0(1:nnz)=0.5d0*(z0(1:nnz)+z0(2:nnz+1)); ! 3D/2D
!
   dxi=dxcj+1; dxj=dxi+nnx; dxci=dxj+1; dxcj=dxci+nnx-1;
   dyi=dycj+1; dyj=dyi+nny; dyci=dyj+1; dycj=dyci+nny-1;
   dzi=dzcj+1; dzj=dzi+nnz-i2d; dzci=dzj+1; dzcj=dzci+nnz-1;
!
   dxall(dxi :dxj) = x0(1:nnx+1); ! append cor metric, then cen metric to global metric array
   dxall(dxci:dxcj)=dx0(1:nnx);
   dyall(dyi :dyj) = y0(1:nny+1);
   dyall(dyci:dycj)=dy0(1:nny);
   dzall(dzi :dzj) = z0(1:nnz+i3d); ! 3D/2D
   dzall(dzci:dzcj)=dz0(1:nnz); ! 3D/2D
!
   x0(1:nnx)=x1(1:nnx); y0(1:nny)=y1(1:nny); z0(1:nnz)=z1(1:nnz); ! recurse : redefine x1 as the fine grid and repeat
!
  enddo
!
! aa:
! write(0,*) dxall-dzall
! write(0,*) '***********'
! write(0,*) dyall-dzall
!
  odxall=1d0/dxall; odyall=1d0/dyall; odzall=1d0/dzall;
!
  deallocate(x0,y0,z0,x1,y1,z1,dx0,dy0,dz0)
!*************************************** FINISHED WITH METRICS *******************************/
  nnx=nx-2; nny=ny-2; nnz=nz-2; ! inner points
! allocate data arrays
  fac2d=(1d0-0.25d0**numlev) /(1d0-0.25d0);
  fac3d=(1d0-0.125d0**numlev)/(1d0-0.125d0);
!
  ! reduce dimensionality for storage in 2D
  if (q2d) fac3d=fac2d
! total length for 3D storage without boundary points
  len3D =nint(nnx*nny*nnz*fac3d); ! 2D/3D
! xy plane is special because in 2D we keep z-ghostpoints at all levels (even though they are not used)
  len2Dxy=nint(nnx*nny*fac2d);
! for other arrays, reduce dimensionality
  if (q2d) then ; fac2d=fac1d ; endif
  len2Dxz=nint(nnx*nnz*fac2d);
  len2Dyz=nint(nny*nnz*fac2d);
!
! coefficient for 3D storage with boundary points
  len3Dbc=len3D+ &! inner points
& 2*(len2Dxy+len2Dxz+len2Dyz)+ &! 2D boundaries
& 4*(nint(fac1d*(nnx+nny+i3d*nnz)))+ &! 1D corner lines (z constant)
& (8+i2d*4)*numlev ! 1D corner points and 1-pt corner "lines" in z-direction
!
! note that the length of coefficient arrays is smaller; this is because the solution array contains boundary points, which are also stored
! write(0,*) numlev, fac1d, fac2d, fac3d !aa
! write(0,*) len3d, nnx*nny*nnz, len3dbc, nx*ny*nz !aa
!
  allocate(allp(len3Dbc),alleps(len3Dbc)) ! 3D solution array and epsilon array (with boundary points)
!
  allocate(allkappa(len3D),allrhs(len3D)) ! 3D coefficient arrays w/o boundary points
  allocate(alleast(len2Dyz), allwest(len2Dyz), allnorth(len2Dxz), allsouth(len2Dxz), allback(len2Dxy), allfront(len2Dxy)) ! 2D bc arrays
!
  allocate(bc_wgt(6,numlev)) ! coefficients for applying BC
! to populate data arrays need to call external subroutine to do a pointer cast
  call multigrid_update_bc()
!
  call multigrid_update_coef(.true., .true., .true., .true.);
!
  if (.not.fatal_warning()) multigrid_initialized=.true.
!
  end subroutine multigrid_init
!*****************************************************************************************************
  subroutine multigrid_done()
  implicit none
  if (multigrid_initialized) then
   if (associated(mcycle).and..not.associated(mcycle, target=vcycle)) deallocate(mcycle)
   deallocate(allp,alleps,allkappa,allrhs,&
& allwest, alleast, allnorth, allsouth, allback, allfront,&
& dxall, odxall, dyall, odyall, dzall, odzall, vcycle)
   multigrid_initialized=.false.
  endif
!
  end subroutine multigrid_done
!*****************************************************************************************************
! this routine handles both 2 and 3 dimensions; however, nz has to be equal to 3 for the 2D case to work
  subroutine multigrid_update_coef(update_eps, update_kappa, update_rhs, update_p)
  use gridsize, only : nx, ny, nz, q2D
  use state
  use bc
  implicit none
  logical :: update_eps, update_kappa, update_rhs, update_p
  integer :: l, inew, iold, inewb, ioldb
  integer :: nnx, nny, nnz, ibc
  real*8, pointer :: z(:)
!
  character(len=21), parameter :: whoami='MULTIGRID_UPDATE_COEF'
! BC MACRO:
! (note that if the bc is zero neumann, the value of bc_wgt does not matter)
!
  if (.not.state_initialized) then
    call warning(whoami, 'State not initialized. Nothing.', 0)
    return
  endif
!
  nnx=nx; nny=ny; nnz=nz;
! store first level
  nullify(z)
  if (update_eps ) then
   call copy3d(eps, alleps, nnx,nny,nnz )
! compute boundary points for epsilon using apply_bc via _EPS_BC MACRO:
  l=1; inewb=1;
  allocate(z(max(nnx-2,nny-2,nnz-2)**2)); z=0d0 ; ! overdimensioned array for applying BC to epsilon at all levels
   if (bc_type(west).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,west,ibc,bc_wgt(west,l))
   if (bc_type(east).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,east,ibc,bc_wgt(east,l))
   if (bc_type(south).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,south,ibc,bc_wgt(south,l))
   if (bc_type(north).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,north,ibc,bc_wgt(north,l))
   if (.not.q2d) then
     if (bc_type(front).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,front,ibc,bc_wgt(front,l))
     if (bc_type(back).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,back,ibc,bc_wgt(back,l))
   endif
!write(0,*) 'eps updated at l',l
  endif
!
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!
  if (update_kappa) call copy3d(kappa(2:nnx+1,2:nny+1,2:nnz+1), allkappa,nnx,nny,nnz)
  if (update_rhs ) call copy3d(rhs (2:nnx+1,2:nny+1,2:nnz+1), allrhs, nnx,nny,nnz)
  if (update_p ) call copy3d(p, allp, nx,ny,nz )
! coarsen to the lower levels
!
  iold=1 ! index at previous level
  ioldb=1 ! index at previous level (arrays with boundaries)
!
  do l=2,numlev
   inew=iold + nnx * nny * nnz ! index at this level
   inewb=ioldb+(nnx+2)*(nny+2)*(nnz+2) ! index at this level (for arrays with boundary points)
! interpolate:
   if (update_kappa) call coarsen(allkappa(iold), allkappa(inew),nnx,nny,nnz, i2d, 0)
!
   if (update_eps) then
!write(0,*) 'attemp to coarsen eps at l',l, ioldb, inewb
    call coarsen(alleps(ioldb), alleps(inewb), nnx,nny,nnz, i2d, 1) ! boundary offset
!write(0,*) 'eps coarsened at l',l
! compute boundary points for epsilon using apply_bc:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnx=nnx/2+2; nny=nny/2+2; nnz=nnz/(2-i2d)+2;
    if (bc_type(west).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,west,ibc,bc_wgt(west,l))
    if (bc_type(east).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,east,ibc,bc_wgt(east,l))
    if (bc_type(south).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,south,ibc,bc_wgt(south,l))
    if (bc_type(north).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,north,ibc,bc_wgt(north,l))
    if (.not.q2d) then
     if (bc_type(front).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,front,ibc,bc_wgt(front,l))
     if (bc_type(back).eq.periodic)then;ibc=periodic;else;ibc=neumann;endif;call apply_bc(alleps(inewb),z,nnx,nny,nnz,back,ibc,bc_wgt(back,l))
    endif
!write(100+l,'(16G20.10)') alleps(1:nnx*nny*nnz)
!write(100+l,'(A)') 'xxxx'
!
    nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(0,*) 'eps bc updated at l',l
   else ! update_eps
    nnx=nnx/2; nny=nny/2; nnz=nnz/(2-i2d);
   endif ! update_eps
! no rhs update, since for coarse levels it is the residual computed after upper level iterations
   iold=inew;
   ioldb=inewb;
  enddo
  if (associated(z)) deallocate(z)
!
  end subroutine multigrid_update_coef
!*****************************************************************************************************
  subroutine multigrid_update_bc()
  use gridsize, only: nx, ny, nz
  use bc
  use output
  use constants
  implicit none
!
  character(len=19), parameter :: whoami='MULTIGRID_UPDATE_BC'
  integer :: nnx, nny, nnz
  integer :: level, i1, j1, k1, ibd
  real*8 :: dplus, dminus
!
  if (.not.bc_initialized) then
    call warning(whoami, 'Boundary conditions not initialized. Nothing done.', 0)
    return
  endif
!
  nnx=nx-2; nny=ny-2; nnz=nz-2;
!
  allwest=0.; alleast=0.; allsouth=0.; allnorth=0.; allfront=0.; allback=0.;
!
! only need to store first level -- other levels are zero; in the case of PBC, bc arrays not used at all
  call copy3d(bcs(west)%p(2:nny+1,2:nnz+1), allwest, nny, nnz, 1)
  call copy3d(bcs(east)%p(2:nny+1,2:nnz+1), alleast, nny, nnz, 1)
!
  call copy3d(bcs(south)%p(2:nnx+1,2:nnz+1), allsouth, nnx, nnz, 1)
  call copy3d(bcs(north)%p(2:nnx+1,2:nnz+1), allnorth, nnx, nnz, 1)
!
  if (ndim.gt.2) then
   call copy3d(bcs(front)%p(2:nnx+1,2:nny+1), allfront, nnx, nny, 1)
   call copy3d(bcs(back )%p(2:nnx+1,2:nny+1), allback, nnx, nny, 1)
  endif
!
! now need to set scaling coefficients for bc data at all levels
  i1=nnx+2; j1=nny+2; k1=1+nnz-i2d+i3d
  do level=1,numlev ! loop over all levels
   do ibd=1,2*ndim ! loop over all dimensions; note, only the z-dimension can have trivial metric arrays
    select case(ibd)
     case(west) ; dplus=dxall(i1); dminus=dxall(nx)
     case(east) ; dplus=dxall(i1+nnx-1); dminus=dxall(2*nx-3)
     case(south); dplus=dyall(j1); dminus=dyall(ny)
     case(north); dplus=dyall(j1+nny-1); dminus=dyall(2*ny-3)
     case(front); dplus=dzall(k1); dminus=dzall(nz)
     case(back) ; dplus=dzall(k1+nnz-1); dminus=dzall(2*nz-3)
    end select
!
    select case(bc_type(ibd))
     case(dirichlet) ; bc_wgt(ibd,level)=two ! applied at boundary
     case(dirichletg) ; bc_wgt(ibd,level)=dplus*2/(dplus+dminus) ! applied at ghost point
     case(neumann) ; bc_wgt(ibd,level)=dplus ! at boundary
     case(neumanng) ; bc_wgt(ibd,level)=dplus ! at ghost point
    end select
   enddo ! ibd
!
   nnx=nnx/2; nny=nny/2; nnz=nnz/(2-i2d);
   i1=i1+3*nnx+1; j1=j1+3*nny+1; k1=k1+3*(nnz-i2d)+i3d;
!
  enddo ! level
!
! write(0,*) bc_wgt(1,:) !aa
 end subroutine multigrid_update_bc
!*************************** HERE WE GO, KIDS *****************************************************!
 subroutine multigrid_solve()
!
  use gridsize, only: nx, ny, nz, communicator, q2d
  use output
  use parser, only: itoa, ftoa
  use bc ! boundary conditions
  use state, only: p ! for copying solution to main array
!
  implicit none
  character(len=15), parameter :: whoami = 'MULTIGRID_SOLVE'
  character(len=100) :: msg
  integer :: msglen
  integer :: level, i1, j1, k1, i3, i3b, i2, j2, k2 ! indices into global arrays
  integer :: i, j, k, m, n
  integer :: nnx, nny, nnz, nall
  integer :: di, icycle
  integer :: imax, jmax, kmax
!
! define multigrid metric arrays: will compute metrics once at all levels in the cpu version since memory tends to be cheaper than flops
!
  real*8, dimension( len3d ) :: alle, allw, alls, alln, allb, allf, allo
  real*8 :: res ( (nx-2) * (ny-2) * (nz-2) ), curr_residual, d
  logical :: qres
!
! interface
! subroutine apply_bc(u,g,nnx,nny,nnz,boundary,bctype)
! real*8 :: u(nnx,nny,nnz), g(*)
! integer :: nnx,nny,nnz,boundary,bctype
! end subroutine apply_bc
! end interface
!
! include 'multigrid_int.h' ! disable argument checking
!
  if (.not.multigrid_initialized) then
   call warning(whoami, 'SOLVER NOT INITIALIZED. ABORT.',0)
   return
  endif
!
! allocate multigrid arrays
! allocate( allw(len3D), alle(len3D), alls(len3D), alln(len3D), allb(len3D), allf(len3D), allo(len3D) )
! populate multigrid arrays
! loop over all levels
  nnx=nx-2; nny=ny-2; nnz=nz-2;
  i1=1; j1=1; k1=1;
  i3=1; i3b=1
!
  do level=1,numlev
   call compute_fd_coef (allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),allo(i3),& ! multigrid coefficients
& alleps(i3b),allkappa(i3), & ! gridded coefficients
& odxall(i1),odxall(i1+nnx+1),odyall(j1),odyall(j1+nny+1),odzall(k1),odzall(k1+nnz-i2d+i3d),& ! grid metrics
& nnx,nny,nnz,i2d) ! dimension
! note that the metrics are normalized by o; furthermore, the rhs must be normalized by o at every coarsening/prolongation
   i3=i3+nnx*nny*nnz;
   i3b=i3b+(nnx+2)*(nny+2)*(nnz+2);
   i1=i1+2*nnx+1; j1=j1+2*nny+1; k1=k1+2*(nnz-i2d)+i3d;
   nnx=nnx/2; nny=nny/2; nnz=nnz/(2-i2d);
  enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mgtimer=timer_start()
  smoothtimer=timer_start()
!
! define initial indices into global arrays
  nnx=nx-2; nny=ny-2; nnz=nz-2;
  i1=1; j1=1; k1=1;
  i2=1; j2=1; k2=1;
  i3=1; i3b=1;
!
! scale rhs by o metric
  nall=nnx*nny*nnz
  do k=i3, i3-1+nall; allrhs(k)=allrhs(k)/allo(k); enddo
  icycle=0;
! update bc
  level=1
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 bc_time=bc_time+timer_elapsed(mgtimer)
! compute initial residual
  if (compute_initial_residual) then
 mgtime=timer_stamp(mgtimer)
   call residual(res, allp, allrhs, allw, alle, alls, alln, allf, allb, allo, nx, ny, nz, q2d); qres=.false.;
 residual_time=residual_time+timer_elapsed(mgtimer)
! define residual normalization (could be made different for subsequent calls)
!
!
  imax=1 ; i=imax ; curr_residual=abs(res(i)/allo(i) );
  do i=2,nall; d=abs(res(i)/allo(i) ); if (d.gt.curr_residual) then ; curr_residual=d; imax=i ; endif ; enddo
!
  if (multigrid_loud) then
   write(msg,'('//realfmt//')') curr_residual; msglen=len_trim(msg);
   call message(whoami, 'Maximum residual after cycle '//itoa(icycle)//' : '//msg(1:msglen)//&
& ' @[ '//&
& itoa(mod(imax-1,nnx)+1 +1)//','//& ! adding one to account for ghost points
& itoa(mod((imax-1)/nnx,nny)+1 +1)//','//&
& itoa((imax-1)/(nnx*nny)+1 +1)//']')
  endif
!
  else
   curr_residual=one
  endif ! compute_initial_residual
  init_residual=curr_residual
!
!
! aa
! open(666,file='resf0.dat',form='unformatted') ; write(666) res; close(666)
! open(666,file='rhsf.dat',form='unformatted') ; write(666) allrhs(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='e.dat',form='unformatted') ; write(666) alle(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='w.dat',form='unformatted') ; write(666) allw(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='s.dat',form='unformatted') ; write(666) alls(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='n.dat',form='unformatted') ; write(666) alln(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='o.dat',form='unformatted') ; write(666) allo(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='f.dat',form='unformatted') ; write(666) allf(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='b.dat',form='unformatted') ; write(666) allb(1:nall); close(666) ! note: RHS normalized by o
! open(666,file='rhsf.dat',form='unformatted') ; write(666) allrhs(1:nall)*allo(1:nall); close(666) ! note: RHS normalized by o
! stop ! aa
! begin multigrid iterations
!
  icycle=1;
  do while (curr_residual.gt.max_residual)
   if (icycle.gt.maxcycle) then
    call warning(whoami, 'Maximum number of cycles reached.',0)
    exit
   endif
!
   do i=1,size(mcycle)
    j=mcycle(i)
!
!***********************************************************************************************
    select case(j);
!***********************************************************************************************
     case(reduce); ! 1) compute residual; 2) update indices 3) interpolate residual onto a coarser grid; 4) update bc
 mgtime=timer_stamp(mgtimer)
      if (qres) then
       call residual(res,allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),allo(i3),nnx+2,nny+2,nnz+2,q2d)
 residual_time=residual_time+timer_stamp(mgtimer)
      endif ! qres
! update indices (some of these updates are not necessary in this version of subroutine)
      i1=i1+2*nnx+1; j1=j1+2*nny+1; k1=k1+2*(nnz-i2d)+i3d;
! 2D bc offsets
      i2=i2+nny*nnz; j2=j2+nnz*nnx; k2=k2+nnx*nny;
! 3D offsets
      i3=i3+nnx*nny*nnz; i3b=i3b+(nnx+2)*(nny+2)*(nnz+2);
! inerpolate
      call coarsen(res, allrhs(i3), nnx, nny, nnz, i2d, 0)
! switch to coarse grid lengths
      nnx=nnx/2; nny=nny/2; nnz=nnz/(2-i2d); ! 3D/2D
! write(0,*) size(allrhs), i3+(nnx)*(nny)*(nnz) - 1 ! aa
! write(0,*) size(allp), i3b+(nnx+2)*(nny+2)*(nnz+2) - 1 ! aa
! scale new rhs by o metric
      do k=i3, i3-1+nnx*nny*nnz; allrhs(k)=allrhs(k)/allo(k); enddo
! initialize solution (to 0, are there other sensible options?)
      call fzero(allp(i3b),nnx+2,nny+2,nnz+2)
! update bc (note that the bc update should be unnecessary in most cases, since we begin from an initial solution of zero)
      level=level+1
!
 coarsen_time=coarsen_time+timer_elapsed(mgtimer)
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 bc_time=bc_time+timer_elapsed(mgtimer)
!***********************************************************************************************
     case(prolong); ! 1) interpolate solution onto finer mesh; 2) add to solution on the finer mesh 3) update bc
 mgtime=timer_stamp(mgtimer)
!
      nnx=2*nnx; nny=2*nny; nnz=(2-i2d)*nnz; ! switch to finer grid dimensions
      di=(nnx+2)*(nny+2)*(nnz+2)
      call refine(allp(i3b-di), allp(i3b), nnx, nny, nnz, i2d)
! update offsets
! update indices (some of these updates are not necessary in this version of subroutine)
      i1=i1-2*nnx-1; j1=j1-2*nny-1; k1=k1-2*(nnz-i2d)-i3d;
! 2D bc offsets
      i2=i2-nny*nnz; j2=j2-nnz*nnx; k2=k2-nnx*nny;
! 3D offsets
      i3=i3-nnx*nny*nnz; i3b=i3b-di;
! update bc
      level=level-1
 refine_time=refine_time+timer_elapsed(mgtimer)
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 bc_time=bc_time+timer_elapsed(mgtimer)
!***********************************************************************************************
     case default; ! this means a positive integer is found --> call smoother
 mgtime=timer_stamp(smoothtimer)
      select case(smoother);
       case(Jacobi_);
       do k=1,j
        call Jacobi(allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),nnx+2,nny+2,nnz+2,omega,q2d)
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(JacobiTiledLomem_);
       do k=1,j
        call JacobiTiledLoMem(allp(i3b),allrhs(i3),alleps(i3b),allkappa(i3),&
& odxall(i1),odxall(i1+nnx+1),odyall(j1),odyall(j1+nny+1),odzall(k1),odzall(k1+nnz-i2d+i3d),&
& nnx+2,nny+2,nnz+2,omega,smooth_tilesize,i2d)
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(JacobiOnTheFly_);
       do k=1,j
        call JacobiOnTheFly(allp(i3b),allrhs(i3),alleps(i3b),allkappa(i3),&
& odxall(i1),odxall(i1+nnx+1),odyall(j1),odyall(j1+nny+1),odzall(k1),odzall(k1+nnz-i2d+i3d),&
& nnx+2,nny+2,nnz+2,omega,i2d)
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(GaussSeidel_);
! write(0,*) 'DBG SMOOTH:'
       do k=1,j
        call GaussSeidel(allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),nnx+2,nny+2,nnz+2,omega,q2d)
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(GaussSeidelReverse_);
! write(0,*) 'DBG SMOOTH:'
       do k=1,j
        call GaussSeidelReverse(allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),nnx+2,nny+2,nnz+2,omega&
& ,q2d)
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(GaussSeidelRB_);
! write(0,*) 'DBG SMOOTH:'
       do k=1,j
        call GaussSeidelRB(allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),nnx+2,nny+2,nnz+2,omega,q2d)
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(GaussSeidelRBTiled_);
       do k=1,j
        call GaussSeidelRBTiled(allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),nnx+2,nny+2,nnz+2,omega&
& ,smooth_tilesize, q2d)
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(GaussSeidelRBTiledLoMem_);
       do k=1,j
        call GaussSeidelRBTiledLoMem(allp(i3b),allrhs(i3),alleps(i3b),allkappa(i3),&
& odxall(i1),odxall(i1+nnx+1),odyall(j1),odyall(j1+nny+1),odzall(k1),odzall(k1+nnz-i2d+i3d),&
& nnx+2,nny+2,nnz+2,omega,smooth_tilesize,i2d)
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
       enddo
!***********************************************************************************************
       case(GSUnrollFromMiddle_);
! unroll as many iterations as possible; when not possible to perform requested number of
! unrolled iterations, decrease unrolling number; repeat steps until the requested number of smoothing
! iterations is performed
       if (unroll_iter.le.-1) then ; k=j ; else ; k=unroll_iter+1 ; endif ! when unroll_iter < 0, attempt to unroll completely
! k is the total number of iterations (unroll_iter=0 corresponds to k=1)
! k cannot be greater than j (in case unroll_iter was specified to be too large)
       k=min(k,j,ishft(nnx-2,-1),ishft(nny-2,-1)); if (.not.q2d) k=min(k,ishft(nnz-2,-1))
       do while (j.gt.0)
!
! write(0,*) level, j, k!,q2d,ishft(nny-2,-1)
        do m=1,j/k ! k=1 corresponds to no unrolling; with k=2, the number of loops should be halved
         call GSInner(allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),nnx+2,nny+2,nnz+2,omega,q2d,k)
         do n=k-1,1,-1; ! at each iteration, update the outermost n points; k is the total number of iterations
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
          call GSOuter(allp(i3b),allrhs(i3),allw(i3),alle(i3),alls(i3),alln(i3),allf(i3),allb(i3),nnx+2,nny+2,nnz+2,omega,q2d,n)
         enddo ! n
! update BC
 mgtime=timer_stamp(mgtimer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nnx=nnx+2; nny=nny+2; nnz=nnz+2;
  call apply_bc(allp(i3b),allwest (i2),nnx, nny, nnz, west, bc_type(west), bc_wgt(west,level) )
  call apply_bc(allp(i3b),alleast (i2),nnx, nny, nnz, east, bc_type(east), bc_wgt(east,level) )
  call apply_bc(allp(i3b),allsouth(j2),nnx, nny, nnz, south, bc_type(south), bc_wgt(south,level))
  call apply_bc(allp(i3b),allnorth(j2),nnx, nny, nnz, north, bc_type(north), bc_wgt(north,level))
  if (.not.q2d) then
   call apply_bc(allp(i3b),allfront(k2),nnx, nny, nnz, front, bc_type(front), bc_wgt(front,level))
   call apply_bc(allp(i3b),allback (k2),nnx, nny, nnz, back, bc_type(back), bc_wgt(back,level) )
  endif
  nnx=nnx-2; nny=nny-2; nnz=nnz-2;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 smooth_time=smooth_time+bc_time ! subtract the BC time from smoother time
 bc_time=bc_time+timer_elapsed(mgtimer)
 smooth_time=smooth_time-bc_time
        enddo ! m
!
        j=mod(j,k); k=j ! remaining iterations
       enddo ! j
!
      end select
 smooth_time=smooth_time+timer_elapsed(smoothtimer)
!***********************************************************************************************
    end select
!***********************************************************************************************
    qres=.true. ! true means that the residual must be recomputed (b/c arrays were modified in the cycle)
   enddo ! end of cycle
! compute residual (update bc first)
! NOTE: the residual calculation below is strictly for the purpose of output;
! we could check less frequently and save a little time
!***********************************************
   if (icycle.gt.convergence_skip_cycles) then
 mgtime=timer_stamp(mgtimer)
   call residual(res,allp,allrhs,allw,alle,alls,alln,allf,allb,allo,nx,ny,nz,q2d); qres=.false.; ! false means residual is already known
 residual_time=residual_time+timer_stamp(mgtimer)
!
!
  imax=1 ; i=imax ; curr_residual=abs(res(i)/allo(i) );
  do i=2,nall; d=abs(res(i)/allo(i) ); if (d.gt.curr_residual) then ; curr_residual=d; imax=i ; endif ; enddo
!
  if (multigrid_loud) then
   write(msg,'('//realfmt//')') curr_residual; msglen=len_trim(msg);
   call message(whoami, 'Maximum residual after cycle '//itoa(icycle)//' : '//msg(1:msglen)//&
& ' @[ '//&
& itoa(mod(imax-1,nnx)+1 +1)//','//& ! adding one to account for ghost points
& itoa(mod((imax-1)/nnx,nny)+1 +1)//','//&
& itoa((imax-1)/(nnx*nny)+1 +1)//']')
  endif
!
   endif
!***********************************************
   icycle=icycle+1
  enddo ! while not solved
!
 mgtime=timer_elapsed_total(mgtimer)
 call message(whoami,'============= TIMING INFORMATION ==============')
 call message(whoami,'Time within multigrid cycle (s) : '//ftoa(mgtime));
!
 call message(whoami,'Coarsening time (s) : '//ftoa(coarsen_time))
 call message(whoami,'Refinement time (s) : '//ftoa(refine_time))
 call message(whoami,'Smoothing time (s)  : '//ftoa(smooth_time))
 call message(whoami,'Residual time (s)   : '//ftoa(residual_time))
 call message(whoami,'BC time (s)         : '//ftoa(bc_time))
 call message(whoami,'Other time (s)      : '//ftoa(mgtime-coarsen_time-refine_time-smooth_time-residual_time-bc_time))
 call message(whoami,'===============================================')
 call timer_erase(mgtimer)
 call timer_erase(smoothtimer)
! copy solution to main array
  call copy3d(allp,p,nx,ny,nz)
! aa
! open(666,file='resf.dat',form='unformatted') ; write(666) res; close(666)
! open(666,file='p.dat',form='unformatted') ; write(666) p; close(666)
! open(666,file='p.dat',form='formatted') ; write(666,*) p(:,:,2); close(666)
!
 end subroutine multigrid_solve
!*****************************************************************************************!
 subroutine parse_multigrid_cycle(str, icycle)
  use ivector
  use output
  use gridsize, only: me, communicator
  use parser
  implicit none
!
  character, parameter :: csmooth(2) = (/'S','s'/);
  character, parameter :: cprolong(2) = (/'P','p'/);
  character, parameter :: ccoarsen(2) = (/'C','c'/);
!
  integer, pointer :: icycle(:)
!
  character(len=21), parameter :: whoami='PARSE_MULTIGRID_CYCLE'
  character(len=*) :: str
  character :: cmd
  integer :: ind, ind2, strlen, iter
  integer :: i
  type (int_vector) :: cycle
!
  ind=1
  strlen=len_trim(str)
  call int_vector_init(cycle)
!
  do while (ind.le.strlen)
!
   cmd=str(ind:ind)
   if (any(csmooth.eq.cmd)) then
! look for number of iterations
    ind2=ind+1;
    do
     cmd=str(ind2:ind2)
     if (ind2.gt.strlen) exit
     if (any(digits2.eq.cmd) .or. (ind2.eq.ind+1.and.cmd.eq.hyphen)) then
       ind2=ind2+1;
     else
      exit
     endif
    enddo
    ind2=ind2-1
!
    if (ind2.gt.ind) then
     iter=atoi(str(ind+1:ind2))
     if (fatal_warning()) exit
     if (iter.le.0) then
      call warning(whoami, 'Number of smoothing iterations must be positive. Nothing done.',-1)
      stop
      exit
     endif
    else
     iter=smooth_iterations
    endif ! ind2
    i=int_vector_add(cycle,iter)
    ind=ind2
   elseif (any(cprolong.eq.cmd)) then
    i=int_vector_add(cycle,prolong)
   elseif (any(ccoarsen.eq.cmd)) then
    i=int_vector_add(cycle,reduce)
   elseif (any(space.eq.cmd)) then
! nothing
   else ! invalid character
    call warning(whoami, 'Invalid character encountered. Nothing done.',-1)
    exit
   endif
   ind=ind+1
  enddo
!
  if (.not.fatal_warning()) then
   allocate(icycle(cycle%last))
   icycle=cycle%i(1:cycle%last)
  endif
!
  call int_vector_done(cycle)
!
 end subroutine parse_multigrid_cycle
!
 end module multigrid
!*****************************************************************************************!
! utility subroutines; they are useful for changing array ranks
! there is no error checking in these routines
!*****************************************************************************************!
! subroutine copy2d(a,b,nx,ny)
! implicit none
! integer :: nx, ny
! real*8, dimension(nx,ny) :: a, b
! b=a
! end subroutine copy2d
!*****************************************************************************************!
 subroutine copy3d(a,b,nx,ny,nz)
 implicit none
 integer :: nx, ny, nz
 real*8, dimension(nx,ny,nz) :: a, b
 b=a
 end subroutine copy3d
!*****************************************************************************************!
 subroutine fzero(a,nx,ny,nz)
 use constants
 implicit none
 integer :: nx, ny, nz
 real*8, dimension(nx,ny,nz) :: a
 a=zero
 end subroutine fzero
!*****************************************************************************************!
 subroutine copy_inner3d(a,b,nx,ny,nz)
 implicit none
 integer :: nx, ny, nz
 real*8, dimension(*) :: a, b
! copy the inner region (i.e. excluding ghostpoints) of array a to array b
 end subroutine copy_inner3d
!*****************************************************************************************!
! Solution prolongation
 subroutine refine(f,c,nx,ny,nz,i2d)
! used to interpolate solution at coarse level
! note that both arrays are dimensioned to include boundary points
 use constants, only: three, nine, twentyseven
 implicit none
 integer :: nx, ny, nz
 integer(kind=KIND('a')) :: i2d
 real*8 :: f(nx+2,ny+2,nz+2), c(nx/2+2,ny/2+2,nz/(2-i2d)+2)
 integer :: nnx, nny, nnz, nxp, nyp, nzp, nnxp, nnyp, nnzp, nnxp2, nnyp2, nnzp2
 real*8 :: coef
 real*8, parameter :: cb=2d0/3d0
!
 integer :: i,j,k,jp,kp,ii,jj,kk
 real*8 :: cwsf, cesf, cwnf, cenf, cwsb, cesb, cwnb, cenb
!
 coef=1d0/(4**(3-i2d)) ! 1/64 in 3D; 1/16 in 2D
!
 nnx=nx/2; nny=ny/2; nnz=nz/(2-i2d);
 nxp=nx+1; nyp=ny+1; nzp=nz+1;
 nnxp=nnx+1; nnyp=nny+1; nnzp=nnz+1;
 nnxp2=nnxp+1; nnyp2=nnyp+1; nnzp2=nnzp+1;
! note: both nnx and nx correspond to the number of inner points; thus, they have different meaning from those in the coarsening routines
! assign corner points to coarse array (this is ad hoc but is needed for the interpolation below)
! *************** corner lines (12)
! k-
! k-direction is first, because the commands hold in 3D and 2D
 c(1 ,1 ,2:nnzp)=c(1 ,2 ,2:nnzp)+c(2 ,1 ,2:nnzp)-c(2 ,2 ,2:nnzp)
 c(nnxp2,1 ,2:nnzp)=c(nnxp2,2 ,2:nnzp)+c(nnxp,1 ,2:nnzp)-c(nnxp,2 ,2:nnzp)
 c(nnxp2,nnyp2,2:nnzp)=c(nnxp2,nnyp ,2:nnzp)+c(nnxp,nnyp2,2:nnzp)-c(nnxp,nnyp,2:nnzp)
 c(1 ,nnyp2,2:nnzp)=c(2 ,nnyp2,2:nnzp)+c(1 ,nnyp ,2:nnzp)-c(2 ,nnyp,2:nnzp)
!
 if (i2d.eq.0) then ! additional boundary points in 3D
! i-
 c(2:nnxp,1 ,1) =c(2:nnxp,1 ,2) +c(2:nnxp,2 ,1) -c(2:nnxp,2 ,2 )
 c(2:nnxp,nnyp2,1) =c(2:nnxp,nnyp2,2) +c(2:nnxp,nnyp,1) -c(2:nnxp,nnyp,2 )
 c(2:nnxp,nnyp2,nnzp2)=c(2:nnxp,nnyp2,nnzp) +c(2:nnxp,nnyp,nnzp2)-c(2:nnxp,nnyp,nnzp)
 c(2:nnxp,1 ,nnzp2)=c(2:nnxp,2 ,nnzp2)+c(2:nnxp,1,nnzp) -c(2:nnxp,2, nnzp)
! j-
 c(1 ,2:nnyp,1 )=c(2 ,2:nnyp,1 )+c(1 ,2:nnyp,2 )-c(2 ,2:nnyp,2 )
 c(1 ,2:nnyp,nnzp2)=c(2 ,2:nnyp,nnzp2)+c(1 ,2:nnyp,nnzp)-c(2 ,2:nnyp,nnzp)
 c(nnxp2,2:nnyp,nnzp2)=c(nnxp ,2:nnyp,nnzp2)+c(nnxp2,2:nnyp,nnzp)-c(nnxp,2:nnyp,nnzp)
 c(nnxp2,2:nnyp,1 )=c(nnxp ,2:nnyp,1 )+c(nnxp2,2:nnyp,2 )-c(nnxp,2:nnyp,2 )
! ************** corner points (8)
 c(1 ,1,1)=cb * ( c(2 ,1,1) + c(1, 2,1) + c(1 ,1,2) ) - c(2 ,2,2)
!
 c(nnxp2,1,1)=cb * ( c(nnxp,1,1) + c(nnxp2,2,1) + c(nnxp2,1,2) ) - c(nnxp,2,2)
 c(1,nnyp2,1)=cb * ( c(1,nnyp,1) + c(1,nnyp2,2) + c(2,nnyp2,1) ) - c(2,nnyp,2)
 c(1,1,nnzp2)=cb * ( c(1,1,nnzp) + c(2,1,nnzp2) + c(1,2,nnzp2) ) - c(2,2,nnzp)
!
 c(nnxp2,nnyp2,1 )=cb*(c(nnxp ,nnyp2,1) +c(nnxp2,nnyp,1) +c(nnxp2,nnyp2,2 ))-c(nnxp,nnyp,2)
 c(1 ,nnyp2,nnzp2)=cb*(c(1 ,nnyp ,nnzp2)+c(1 ,nnyp2,nnzp )+c(2 ,nnyp2,nnzp2))-c(2 ,nnyp,nnzp)
 c(nnxp2,1 ,nnzp2)=cb*(c(nnxp2,1 ,nnzp )+c(nnxp ,1 ,nnzp2)+c(nnxp2,2 ,nnzp2))-c(nnxp,2 ,nnzp)
!
 c(nnxp2,nnyp2,nnzp2)=cb*(c(nnxp ,nnyp2,nnzp2)+c(nnxp2,nnyp,nnzp2) +c(nnxp2,nnyp2,nnzp ))-c(nnxp,nnyp,nnzp)
!
!
!
! calculate internal points (note that we make use of all points, including the corner points)
 do k=1,nnz+1
  kk=k*2; kp=k+1
  do j=1,nny+1
   jj=j*2 ; jp=j+1 ; cwsf=coef*c(1,j,k) ; cwnf=coef*c(1,jp,k); cwsb=coef*c(1,j,kp) ; cwnb=coef*c(1,jp,kp);
! jj=j*2 ; jp=j+1 ; cwsf=c(1,j,k) ; cwnf=c(1,jp,k); cwsb=c(1,j,kp) ; cwnb=c(1,jp,kp);
   do i=1,nnx+1 ; ii=i*2
    cesf=coef*c(i+1,j,k) ; cenf=coef*c(i+1,jp,k) ; cesb=coef*c(i+1,j,kp) ; cenb=coef*c(i+1,jp,kp);
! cesf=c(i+1,j,k) ; cenf=c(i+1,jp,k) ; cesb=c(i+1,j,kp) ; cenb=c(i+1,jp,kp);
!
    f(ii,jj,kk) =f(ii,jj,kk) + (cwsf + (cesf + cwnf + cwsb)*three + (cenf + cesb + cwnb)*nine + cenb*twentyseven)!*coef
    f(ii-1,jj,kk) =f(ii-1,jj,kk) + (cesf + (cwsf + cenf + cesb)*three + (cwnf + cwsb + cenb)*nine + cwnb*twentyseven)!*coef
    f(ii,jj-1,kk) =f(ii,jj-1,kk) + (cwnf + (cenf + cwsf + cwnb)*three + (cesf + cenb + cwsb)*nine + cesb*twentyseven)!*coef
    f(ii-1,jj-1,kk)=f(ii-1,jj-1,kk)+ (cenf + (cwnf + cesf + cenb)*three + (cwsf + cwnb + cesb)*nine + cwsb*twentyseven)!*coef
!
    f(ii,jj,kk-1) =f(ii,jj,kk-1) + (cwsb + (cesb + cwnb + cwsf)*three + (cenb + cesf + cwnf)*nine + cenf*twentyseven)!*coef
    f(ii-1,jj,kk-1) =f(ii-1,jj,kk-1) + (cesb + (cwsb + cenb + cesf)*three + (cwnb + cwsf + cenf)*nine + cwnf*twentyseven)!*coef
    f(ii,jj-1,kk-1) =f(ii,jj-1,kk-1) + (cwnb + (cenb + cwsb + cwnf)*three + (cesb + cenf + cwsf)*nine + cesf*twentyseven)!*coef
    f(ii-1,jj-1,kk-1)=f(ii-1,jj-1,kk-1)+ (cenb + (cwnb + cesb + cenf)*three + (cwsb + cwnf + cesf)*nine + cwsf*twentyseven)!*coef
!
    cwsf=cesf ; cwnf=cenf; cwsb=cesb ; cwnb=cenb;
   enddo ! i
  enddo ! j
 enddo ! k
 else ! i2d = 1
 do j=1,nny+1
  jj=j*2 ; jp=j+1 ; cwsf=coef*c(1,j,2) ; cwnf=coef*c(1,jp,2);
  do i=1,nnx+1 ; ii=i*2
   cesf=coef*c(i+1,j,2) ; cenf=coef*c(i+1,jp,2);
!
   f(ii,jj,2) =f(ii,jj,2) + cwsf + (cesf + cwnf)*three + cenf*nine
   f(ii-1,jj-1,2)=f(ii-1,jj-1,2) + cenf + (cesf + cwnf)*three + cwsf*nine
   f(ii,jj-1,2) =f(ii,jj-1,2) + cwnf + (cenf + cwsf)*three + cesf*nine
   f(ii-1,jj,2) =f(ii-1,jj,2) + cesf + (cenf + cwsf)*three + cwnf*nine
!
   cwsf=cesf ; cwnf=cenf
  enddo
 enddo
!
 endif ! i2d
!
!
 end subroutine refine
!*****************************************************************************************!
 subroutine refine2d_vec(f,c,nx,ny)
! used to interpolate solution at coarse level
! note that both arrays are dimensioned to include boundary points
 implicit none
 integer :: nx, ny
 real*8 :: f(nx+2,ny+2), c(nx/2+2,ny/2+2)
 integer :: nnx, nny, nxp, nyp
 real*8 :: coef=0.0625d0 ! 1/16
 nnx=nx/2; nny=ny/2
 nxp=nx+1; nyp=ny+1
! note: both nnx and nx correspond to the number of inner points; thus, they have different meaning from those in the coarsening routines
! assign corner points to coarse array (ad hoc and not really needed, but see comments below)
! c(1,1) =0.5d0*(c(1,2) +c(2,1))
! c(nnx+2,1) =0.5d0*(c(nnx+1,2) +c(nnx+2,2))
! c(1,nny+2) =0.5d0*(c(2,nny+2) +c(1,nny+1))
! c(nnx+2,nny+2)=0.5d0*(c(nnx+1,nny+2) +c(nnx+2,nny+1))
!
! calculate internal points (note that we make use of all points, even the corner points, so make sure they are not undefined or inf/NaN)
 f(2:nxp:2,2:nyp:2) = f(2:nxp:2,2:nyp:2) &
                      +coef*( c(1:nnx,1:nny) &
                      + 3d0*( c(2:nnx+1,1:nny) &
                      + c(1:nnx,2:nny+1) &
                      + 3d0* c(2:nnx+1,2:nny+1)));
!
 f(3:nxp:2,2:nyp:2) = f(3:nxp:2,2:nyp:2) &
                      +coef*( c(3:nnx+2,1:nny) &
                      + 3d0*( c(2:nnx+1,1:nny) &
                      + c(3:nnx+2,2:nny+1) &
                      + 3d0* c(2:nnx+1,2:nny+1)));
!
 f(2:nxp:2,3:nyp:2) = f(2:nxp:2,3:nyp:2) &
                      +coef*( c(1:nnx,3:nny+2) &
                      + 3d0*( c(1:nnx,2:nny+1) &
                      + c(2:nnx+1,3:nny+2) &
                      + 3d0* c(2:nnx+1,2:nny+1)));
!
 f(3:nxp:2,3:nyp:2) = f(3:nxp:2,3:nyp:2) &
                      +coef*( c(3:nnx+2,3:nny+2) &
                      + 3d0*( c(3:nnx+2,2:nny+1) &
                      + c(2:nnx+1,3:nny+2) &
                      + 3d0* c(2:nnx+1,2:nny+1)));
! note: it may be more efficient (albeit less clear) to rewrite as a series of nested loops
!
 end subroutine refine2d_vec
!*****************************************************************************************!
! Residual coarsening
 subroutine coarsen(f,c,nx,ny,nz,i2d,ibc)
 implicit none
! note that nx,ny,nz correspond to the number of inner points (in 2D, nnz = 1)
 integer(kind=KIND('a')) :: i2d
 integer(kind=KIND('a')), intent(in) :: ibc ! include optional offset for arrays with boundaries (if ibc = 1)
 integer :: nx, ny, nz, nnx, nny, nnz
 real*8 :: f(nx+2*ibc,ny+2*ibc,nz+2*ibc), c(nx/2+2*ibc,ny/2+2*ibc,nz/(2-i2d)+2*ibc)
 integer :: ib,ie,je,ke
 real*8 :: coef
 integer :: i,j,k,im,jm,km
!
 ib=1+ibc; nnx=nx/2+ibc; nny=ny/2+ibc; nnz=nz/(2-i2d)+ibc
!
 if (i2d.eq.0) then
 coef=0.125d0
!############################ serial calculation ###########################
 km=ib
 do k=ib,nnz ; jm=ib ;
!
  do j=ib,nny ; im=ib ;
!
   do i=ib,nnx;
    c(i,j,k)=coef * ( f(im,jm,km) +f(im+1,jm,km) +f(im+1,jm+1,km) +f(im,jm+1,km) + &
& f(im,jm,km+1)+f(im+1,jm,km+1)+f(im+1,jm+1,km+1)+f(im,jm+1,km+1) )
    im=im+2
   enddo ! i
   jm=jm+2
  enddo ! j
  km=km+2
 enddo ! k
 else
 coef=0.25d0
!################################### serial calculation ################################
 jm=ib;
 do j=ib,nny ; im=ib ;
!
  do i=ib,nnx;
   c(i,j,ib)=coef * ( f(im,jm,ib)+f(im+1,jm,ib)+f(im+1,jm+1,ib)+f(im,jm+1,ib) )
   im=im+2
  enddo
  jm=jm+2
 enddo
!
 endif
!
 end subroutine coarsen
!*****************************************************************************************!
 subroutine coarsen2d_vec(fine,coarse,nx,ny)
! used for coarsening the residual
! note that nx,ny correspond to the number of inner points
 implicit none
 integer :: nx, ny
 real*8 :: fine(nx,ny), coarse(nx/2,ny/2)
 coarse=0.25d0*(&
   fine(1::2,1::2)+&
   fine(2::2,1::2)+&
   fine(1::2,2::2)+&
   fine(2::2,2::2)&
  )
 end subroutine coarsen2d_vec
!*****************************************************************************************!
! Residual calculation
 subroutine residual(r,p,rhs,w,e,s,n,f,b,o,nx,ny,nz,q2d)
! note that nx,ny,nz include ghost points
! note, also, that q2d is optional (since f=b=0 in 2D) but included to (hopefully) increase speed
 implicit none
 integer :: nx, ny, nz
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: r,rhs,w,e,s,n,f,b,o ! note that the metrics e--o are variable (because they include eps)
 integer :: nnx, nny, nnz
 integer :: nxp, nyp, nzp
 logical :: q2d
!
 integer :: i,j,k,im,jm,km
 real*8 :: pw,pe,po
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 if (q2d) then
! compute in a loop (a la GS.src)
!##################################################################################
  do k=2,nzp; km=k-1;
   do j=2,nyp; jm=j-1; po=p(2,j,k); pw=p(1,j,k);
    do i=2,nxp; im=i-1;
     pe=p(i+1,j,k)
     r(im,jm,km)=o(im,jm,km) * ( rhs(im,jm,km) -&
                                   w(im,jm,km)*pw - e(im,jm,km)*pe -&
                                   s(im,jm,km)*p(i,j-1,k) - n(im,jm,km)*p(i,j+1,k) -&
                                   po )
     pw=po; po=pe;
  enddo ; enddo ; enddo ! i ,j, k
 else ! 3D
!##################################################################################
  do k=2,nzp; km=k-1;
   do j=2,nyp; jm=j-1; po=p(2,j,k); pw=p(1,j,k);
    do i=2,nxp; im=i-1;
     pe=p(i+1,j,k)
     r(im,jm,km)=o(im,jm,km) * ( rhs(im,jm,km) -&
                                   w(im,jm,km)*pw - e(im,jm,km)*pe -&
                                   s(im,jm,km)*p(i,j-1,k) - n(im,jm,km)*p(i,j+1,k) -&
                                   f(im,jm,km)*p(i,j,k-1) - b(im,jm,km)*p(i,j,k+1) -&
                                   po )
     pw=po; po=pe;
  enddo ; enddo ; enddo ! i ,j, k
!
 endif
!
 end subroutine residual
!*****************************************************************************************!
 subroutine residual2d_vec(r,p,rhs,w,e,s,n,o,nx,ny)
! note that nx,ny include ghost points
 implicit none
 integer :: nx, ny
 real*8, dimension(nx,ny) :: p
 real*8, dimension(nx-2,ny-2) :: r,rhs,w,e,s,n,o ! note that the metrics e--o are variable (because they include eps)
 integer :: nnx, nny
 integer :: nxp, nyp
!
 nnx=nx-2; nny=ny-2;
 nxp=nnx+1; nyp=nny+1;
!
 r=o*(rhs-&
      w*p(1:nnx,2:nyp)-&
      e*p(3:nx ,2:nyp)-&
      s*p(2:nxp,1:nny)-&
      n*p(2:nxp,3:ny )-&
        p(2:nxp,2:nyp))
!
 end subroutine residual2d_vec
!*****************************************************************************************!
! Jacobi smoothers
 subroutine Jacobi(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,q2d)
 use bc ! boundary conditions
! note that nx,ny,nz include ghost points
! note that rhs passed into this routine is normalized by o (see below)
 implicit none
 integer :: nx, ny, nz
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p &
& ,q ! temporary array
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b ! note that the metrics e--b are variable (because they include eps & kappa; oo=1/o)
 integer :: nnx, nny, nnz
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km
 real*8 :: pw, pe, po
 logical :: q2d
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
!
!#define __SAVEP
!##################################################################################
!
  if (q2d) then
!
   do k=2,nzp; km=k-1;
    do j=2,nyp; jm=j-1;
     do i=2,nxp; im=i-1;
      q(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
                                 s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                                 p(i,j,k) - rhs(im,jm,km) )
     enddo ; enddo ; enddo ! i,j,k
!
  else ! 3D
!##################################################################################
   do k=2,nzp; km=k-1;
    do j=2,nyp; jm=j-1;
     do i=2,nxp; im=i-1;
      q(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
                                 s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                                 f(im,jm,km)*p(i,j,km) + b(im,jm,km)*p(i,j,k+1) +&
                                 p(i,j,k) - rhs(im,jm,km) )
   enddo ; enddo ; enddo ! i,j,k
!
!##################################################################################
!
  endif ! q2d
!
  do k=2,nzp; do j=2,nyp; do i=2,nxp ; p(i,j,k)=q(i,j,k) ; enddo ; enddo ; enddo
!
 end subroutine Jacobi
!**********************************************************************************************************************!
! EXPERIMENTAL (now, very slow) ----------- !!!!!!!!!!!!!!!!!!!!
 subroutine JacobiOnTheFly(p,rhs,eps,kappa,odxcen,odxcor,odycen,odycor,odzcen,odzcor,nx,ny,nz,dt,i2d)
! compute metrics 'on the fly'
 use constants, only : one
 implicit none
 integer :: nx, ny, nz
 integer(kind=KIND('a')) :: i2d
 real*8 :: dt
 real*8 :: odxcen(nx-1), odxcor(nx-2)
 real*8 :: odycen(ny-1), odycor(ny-2)
 real*8 :: odzcen(nz-1-i2d-i2d), odzcor(nz-2-i2d) ! 3D/2D
 real*8, dimension(nx,ny,nz) :: p &
& ,q ! temporary array
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs, eps, kappa
 real*8 :: w, e, s, n, f, b, o
 real*8 :: pw, pe, po, d
 integer :: nnx, nny, nnz, ntot
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
! if (i2d.eq.1) then
! else ! 3D
!#ifdef __SAVEP
! prefactor quarter instead of half if metrics are doubled (see below)
 d = 0.25d0 * dt / (maxval(odxcor)**2+maxval(odycor)**2+maxval(odzcor)**2)
!write(0,*)"------------ ",d
!##################################################################################
  do k=2,nzp; km=k-1;
   do j=2,nyp; jm=j-1; po=p(2,j,k); pw=p(1,j,k);
    do i=2,nxp; im=i-1;
     pe=p(i+1,j,k)
! compute metrics
     w = odxcor(im)*odxcen(im)*(eps(im,jm,km)+eps(max(im-1,1),jm, km )) ! beware of the 0.5 factor
     e = odxcor(im)*odxcen(i )*(eps(im,jm,km)+eps(min(i,nnx), jm, km ))
     s = odycor(jm)*odycen(jm)*(eps(im,jm,km)+eps(im, max(jm-1,1),km ))
     n = odycor(jm)*odycen(j )*(eps(im,jm,km)+eps(im, min(j,nny), km ))
     f = odzcor(km)*odzcen(km)*(eps(im,jm,km)+eps(im, jm, max(km-1,1)))
     b = odzcor(km)*odzcen(k )*(eps(im,jm,km)+eps(im, jm, min(k,nnz) ))
     o = -(w+e+s+n+f+b) + 2d0 * kappa(im,jm,km) ;! o=one/o; ! add twice kappa if missing 0.5 above
!write(0,*) d*o
!stop
     q(i,j,k)=po + d* (o*(po-rhs(im,jm,km))+&
& w*pw + &
& e*pe + &
& s*p(i,jm,k) + &
& n*p(i,j+1,k) + &
& f*p(i,j,km) + &
& b*p(i,j,k+1) &
& )
! q(i,j,k)=po - dt * ( po + &
!& o*( w*pw + &
!& e*pe + &
!& s*p(i,jm,k) + &
!& n*p(i,j+1,k) + &
!& f*p(i,j,km) + &
!& b*p(i,j,k+1) &
!& ) &
!& -rhs(im,jm,km) )
      pw=po; po=pe;
     enddo ; enddo ; enddo ! i,j,k
!#else
! endif ! i2d
!##################################################################################
  do k=2,nzp; do j=2,nyp; do i=2,nxp ; p(i,j,k)=q(i,j,k) ; enddo ; enddo ; enddo
!
 end subroutine JacobiOnTheFly
!**********************************************************************************************************************!
 subroutine JacobiTiledLoMem(p,rhs,eps,kappa,odxcen,odxcor,odycen,odycor,odzcen,odzcor,nx,ny,nz,dt,ts,i2d)
! conceptual GPU template
! compute metrics 'on the fly'
 use constants, only : one, third
 implicit none
 integer :: nx, ny, nz
 integer(kind=KIND('a')) :: i2d
 real*8 :: dt
 real*8 :: odxcen(nx-1), odxcor(nx-2)
 real*8 :: odycen(ny-1), odycor(ny-2)
 real*8 :: odzcen(nz-1-i2d-i2d), odzcor(nz-2-i2d) ! 3D/2D
 real*8, dimension(nx,ny,nz) :: p, eps
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs, kappa
 real*8 :: w, e, s, n, f, b, o
! real*8 :: d1, d2, d3 ! dummies
 integer :: nnx, nny, nnz, ntot
 integer :: nxp, nyp, nzp
 integer :: i, j, k, im, jm, km, ii, jj, kk, jeps, m
 integer :: iter ! inner iterations
 integer :: offset
!
 integer :: tilesize, ts ! conceptually, might be related to CPU cache size
 integer :: tiley, tilez, tx, ty, tz, ntx, nty, ntz, remainder ! sizes do not include ghost cells
 real*8, allocatable, dimension(:,:,:) :: rhs1, eps1, kappa1
 real*8, pointer, dimension(:,:,:) :: p1, q1
 real*8, allocatable, dimension(:) :: odxcen1, odxcor1, odycen1, odycor1, odzcen1, odzcor1 ! local metrics
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 ntot=nnx*nny*nnz
!
 if (i2d.eq.1) then
  tilesize=min(max(ts,4),ntot) ! tile should not be larger than number of points
! determine tile dimensions (reuse some variables)
  tilez=INT(sqrt(one*tilesize));
  offset=NOT(ISHFTC(1,-1)) ! upper bound -- largest representable integer : all ones except 1st bit
! vary tile sizes a little to get a better cover (heuristic)
ji_loop:&
  do j=0,1; jm = tilez + j; nty = (nny/jm) + min(1,mod(nny,jm)); ! number of tiles in y-direction
   do i=0,2; im = tilesize / jm + i ; ntx = (nnx/im) + min(1,mod(nnx,im)); ! number of tiles in x-direction
    remainder = ntx*nty*im*jm - ntot ! number of elements by which the size of the cover exceeds grid size
!aa
! write(0,*) im, jm, ntx, nty, remainder, ntot
    if (offset.gt.remainder) then
     tx=im; ty=jm;
     offset=remainder
    endif
    if (remainder.eq.0) exit ji_loop
  enddo; enddo ji_loop
!
! allocate data, now that tile size is known (pad epsilon for easy corner interpolation)
  allocate(rhs1(tx,ty,1),eps1(tx+2,ty+2,1),kappa1(tx,ty,1),p1(tx+2,ty+2,1),q1(tx+2,ty+2,1))
  allocate(odxcen1(tx+1),odxcor1(tx),odycen1(ty+1),odycor1(ty))
! slide tile over all points and perform RB relaxation within each tile
  k=2 ! 2D
  km=k-1
  jj=1;
write(0,*) "TILES:",tx, ty ! aa
  do while (jj.le.nny) ; nty=min(ty,nny-jj+1); ii=1
! local y-metrics
   odycen1(1:nty+1)=odycen(jj:jj+nty); odycor1(1:nty)=odycor(jj:jj+nty-1); ! center-to-center metric has extra element
   do while (ii.le.nnx) ; ntx=min(tx,nnx-ii+1); ! number of elements to read
! local x-metrics
    odxcen1(1:ntx+1)=odxcen(ii:ii+ntx); odxcor1(1:ntx)=odxcor(ii:ii+ntx-1);
!
! populate tile
!
    jm=jj
    do j=1,nty; im=ii;
     do i=1,ntx
      rhs1(i,j,1)=rhs(im,jm,km); eps1(i,j,1)=eps(im,jm,km); kappa1(i,j,1)=kappa(im,jm,km); p1(i,j,1)=p(im,jm,k)
      im=im+1;
     enddo ; ! i
! : last inner point and ghost point
     do i=ntx+1, ntx+2;
      p1(i,j,1)=p(im,jm,k); eps1(i,j,1)=eps(im,jm,k); ! remaining boundary points
      im=im+1
     enddo ; ! i
     jm=jm+1;
    enddo ! j
! : bottom inner points
! epsilon: bottom ghost points
    do j=nty+1, nty+2 ; im=ii;
     do i=1,ntx+2
      p1(i,j,1)=p(im,jm,k); eps1(i,j,1)=eps(im,jm,k);
      im=im+1;
     enddo ;
     jm=jm+1;
! p: bottom ghost points
    enddo ! j
!************************************************************************************************
! smooth within tile (iteration involves local data only)
!************************************************************************************************
    do j=2,nty+1; jm=j-1;
     do i=2,ntx+1; im=i-1;
! compute metrics
! d1=eps1(im,jm,1)
! d2=odxcor1(im)
! d3=odycor1(jm)
      w = odxcor1(im)*odxcen1(im)*(eps1(i,j,1)+eps1(im, j, 1)) ! beware of 0.5 (included using o normalization)
      e = odxcor1(im)*odxcen1(i )*(eps1(i,j,1)+eps1(i+1,j, 1))
      s = odycor1(jm)*odycen1(jm)*(eps1(i,j,1)+eps1(i, jm, 1))
      n = odycor1(jm)*odycen1(j )*(eps1(i,j,1)+eps1(i, j+1,1))
! w = d2*odxcen1(im)*(d1+eps1(im-1,jm,1)) ! 0.5 included using o normalization
! e = d2*odxcen1(i )*(d1+eps1(i ,jm,1))
! s = d3*odycen1(jm)*(d1+eps1(im,jm-1,1))
! n = d3*odycen1(j )*(d1+eps1(im ,j ,1))
      o = -(w+e+s+n) + kappa(im,jm,1) * 2 ; o=one/o ; ! add twice kappa if missing 0.5 above
!
      q1(i,j,1)=p1(i,j,1) - dt * ( p1(i,j,1) &
& + o*( w*p1(im,j,1) + e*p1(i+1,j,1) + s*p1(i,jm,1) + n*p1(i,j+1,1) ) &
& - rhs1(im,jm,1) )
! NOTE: rhs above has been divided by o, which means that o is already known;
! this routine is _just a test of computing metrics on-the-fly vs. memory access_
     enddo ; enddo ; ! i,j
!************************************************************************************************
! put solution back into global array (inner points only)
    jm=jj+1
    do j=2,nty+1; im=ii+1; do i=2,ntx+1
     p(im,jm,k)=q1(i,j,1); im=im+1; enddo ; jm=jm+1 ;
    enddo
! shift tile
    ii=ii+tx
   enddo ! ii
   jj=jj+ty
  enddo ! jj
! free tile memory
  deallocate(eps1, kappa1, rhs1, p1, q1)
  deallocate(odxcen1, odxcor1, odycen1, odycor1)
!
 else ! i2d
!**********************************************************************************************************
! 3D
  tilesize=min(max(ts,8),ntot) ! tile should not be larger than number of points
  tilez=INT( (one*tilesize)**third ) ;
  offset=1 ;
  offset=NOT(ISHFTC(offset,-1)) ! upper bound -- largest representable integer : all ones except 1st bit
! vary tile sizes a little to get a better cover (heuristic)
 kji_loop:&
  do k=0,1; km = tilez + k; ntz = (nnz/km) + min(1,mod(nnz,km)); ! number of tiles in z-direction
   tiley = INT(sqrt(one*tilesize/km))
   do j=0,1; jm = tiley + j ; nty = (nny/jm) + min(1,mod(nny,jm)); ! number of tiles in y-direction
    do i=0,2; im = tilesize/km/jm + i; ntx = (nnx/im) + min(1,mod(nnx,im)); ! number of tiles in x-direction
     remainder = ntx*nty*ntz*im*jm*km - ntot ! number of elements by which the size of the cover exceeds grid size
!aa
! write(0,*) im, jm, km, ntx, nty, ntz, remainder, ntot
!
     if (offset.gt.remainder) then
      tx=im; ty=jm; tz=km
      offset=remainder
     endif
     if (remainder.eq.0) exit kji_loop ! break out of three nested loops since we cannot do better than zero !
  enddo; enddo; enddo kji_loop ! i, j, k
!
! allocate data, now that tile size is known (pad epsilon for easy corner interpolation)
  allocate(rhs1(tx,ty,tz),eps1(tx+2,ty+2,tz+2),kappa1(tx,ty,tz),p1(tx+2,ty+2,tz+2),q1(tx+2,ty+2,tz+2))
  allocate(odxcen1(tx+1),odxcor1(tx),odycen1(ty+1),odycor1(ty),odzcen1(tz+1),odzcor1(tz))
! slide tile over all points and perform RB relaxation within each tile
! NOTE: REUSING NTX, NTY, NYZ: they now correspond to the number of points within tile
  kk=1
  do while (kk.le.nnz) ; ntz=min(tz,nnz-kk+1); jj=1
! local z-metrics
   odzcen1(1:ntz+1)=odzcen(kk:kk+ntz); odzcor1(1:ntz)=odzcor(kk:kk+ntz-1); ! center-to-center metric has extra element
   do while (jj.le.nny) ; nty=min(ty,nny-jj+1); ii=1
! local y-metrics
    odycen1(1:nty+1)=odycen(jj:jj+nty); odycor1(1:nty)=odycor(jj:jj+nty-1);
    do while (ii.le.nnx) ; ntx=min(tx,nnx-ii+1); ! number of elements to read
! local x-metrics
     odxcen1(1:ntx+1)=odxcen(ii:ii+ntx); odxcor1(1:ntx)=odxcor(ii:ii+ntx-1);
! populate tile
!*********************************************************************************************************************
     km=kk
     do k=1,ntz; jm=jj;
      do j=1,nty; im=ii;
       do i=1,ntx ! all inner points below
        rhs1(i,j,k)=rhs(im,jm,km); eps1(i,j,k)=eps(im,jm,km); kappa1(i,j,k)=kappa(im,jm,km); p1(i,j,k)=p(im,jm,km)
        im=im+1;
       enddo ; ! i
! the remaining eps- and p-points are ghost points
! : last inner point and right ghost point
       do i=ntx+1,ntx+2 ; eps1(i,j,k)=eps(im,jm,km) ; p1(i,j,k)=p(im,jm,km) ; im=im+1; enddo ; ! remaining inner and boundary points
       jm=jm+1;
      enddo ! j
!
      do j=nty+1, nty+2; im=ii;
       do i=1,ntx+2 ; eps1(i,j,k)=eps(im,jm,km) ; p1(i,j,k)=p(im,jm,km) ; im=im+1; enddo ; ! remaining inner and boundary points
       jm=jm+1
      enddo ! j
      km=km+1
     enddo ! k
!
     do k=ntz+1, ntz+2; jm=jj;
      do j=1,nty+2; im=ii;
       do i=1,ntx+2 ; eps1(i,j,k)=eps(im,jm,km) ; p1 (i,j,k)=p(im,jm,km) ; im=im+1; enddo ;
       jm=jm+1;
      enddo ! j
      km=km+1
     enddo ! k
!
!****************************************************************************************************************
! smooth within tile (iteration involves local data only)
!****************************************************************************************************************
!#define __INNER
     do k=2,ntz+1; km=k-1;
      do j=2,nty+1; jm=j-1;
       do i=2,ntx+1; im=i-1;
! compute metrics
        w = odxcor1(im)*odxcen1(im)*(eps1(i,j,k)+eps1(im, j, k)) ! be aware of 0.5 factor
        e = odxcor1(im)*odxcen1(i )*(eps1(i,j,k)+eps1(i+1,j, k))
        s = odycor1(jm)*odycen1(jm)*(eps1(i,j,k)+eps1(i, jm, k))
        n = odycor1(jm)*odycen1(j )*(eps1(i,j,k)+eps1(i, j+1,k))
        f = odzcor1(km)*odzcen1(km)*(eps1(i,j,k)+eps1(i, j, km))
        b = odzcor1(km)*odzcen1(k )*(eps1(i,j,k)+eps1(i, j, k+1))
        o = -(w+e+s+n+f+b) + kappa(im,jm,km) * 2 ; o=one/o ! add twice kappa if missing 0.5 above
!
        q1(i,j,k)=p1(i,j,k) - dt * ( p1(i,j,k) &
& + o*( w*p1(im,j,k) + e*p1(i+1,j,k) + s*p1(i,jm,k) + n*p1(i,j+1,k) + f*p1(i,j,km) + b*p1(i,j,k+1))&
& - rhs1(im,jm,km) )
     enddo ; enddo ; enddo ! i,j,k
!************************************************************************************************
! put solution back into global array (inner points only)
     km=kk+1
     do k=2,ntz+1; jm=jj+1; do j=2,nty+1; im=ii+1; do i=2,ntx+1
                                                     p(im,jm,km)=q1(i,j,k);
                                                     im=im+1;
                                                    enddo ;
                              jm=jm+1 ;
                             enddo
      km=km+1
     enddo
!
! shift tile
     ii=ii+tx
    enddo ! ii
    jj=jj+ty
   enddo ! jj
   kk=kk+tz
  enddo ! kk
! free tile memory
  deallocate(eps1, kappa1, rhs1, p1, q1)
  deallocate(odxcen1, odxcor1, odycen1, odycor1, odzcen1, odzcor1)
!
 endif
!
 end subroutine JacobiTiledLoMem
!*****************************************************************************************!
! Gauss-Seidel smoothers
 subroutine GaussSeidel(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,q2d)
! note that nx,ny,nz include ghost points
! note that rhs passed into this routine is normalized by o (see below)
 implicit none
 integer :: nx, ny, nz
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b ! note that the metrics e--b are variable (because they include eps & kappa; oo=1/o)
 integer :: nnx, nny, nnz
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km
 logical :: q2d
 real*8 :: pw, pe, po
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 if (q2d) then
!
!##################################################################################
  k=2; km=k-1;
  do j=2,nyp; jm=j-1; po=p(2,j,k); pw=p(1,j,k);
   do i=2,nxp; im=i-1;
    pe=p(i+1,j,k)
    po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,jm,k) + n(im,jm,km)*p(i,j+1,k) +&
                               po - rhs(im,jm,km) )
    p(i,j,k)=po;
    pw=po; po=pe;
   enddo ; enddo ;
!
 else ! q2d
!
!##################################################################################
 do k=2,nzp; km=k-1;
  do j=2,nyp; jm=j-1; po=p(2,j,k); pw=p(1,j,k);
   do i=2,nxp; im=i-1;
    pe=p(i+1,j,k)
    po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,jm,k) + n(im,jm,km)*p(i,j+1,k) +&
                               f(im,jm,km)*p(i,j,km) + b(im,jm,km)*p(i,j,k+1) +&
                               po - rhs(im,jm,km) )
    p(i,j,k)=po;
    pw=po; po=pe;
   enddo
  enddo
 enddo ! k
!
 endif
! may be better off unrolling, or reusing solution values already fetched from memory (e.g. pcen=p(i,j,k); pright=p(i+1,j,k) ... pleft=pcen, pcen=pright)
 end subroutine GaussSeidel
!*****************************************************************************************!
 subroutine GaussSeidelReverse(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,q2d)
! note that nx,ny,nz include ghost points
! note that rhs passed into this routine is normalized by o (see below)
 implicit none
 integer :: nx, ny, nz
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b ! note that the metrics e--b are variable (because they include eps & kappa; oo=1/o)
 integer :: nnx, nny, nnz
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km
 logical :: q2d
 real*8 :: pw, pe, po
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 if (q2d) then
!
!##################################################################################
  k=2; km=k-1;
  do j=nyp,2,-1; jm=j-1; po=p(nxp,j,k); pe=p(nxp+1,j,k);
   do i=nxp,2,-1; im=i-1;
    pw=p(im,j,k)
    po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,jm,k) + n(im,jm,km)*p(i,j+1,k) +&
                               po - rhs(im,jm,km) )
    p(i,j,k)=po;
    pe=po; po=pw;
   enddo ; enddo ;
!
 else ! q2d
!
!##################################################################################
 do k=nzp,2,-1; km=k-1;
  do j=nyp,2,-1; jm=j-1; po=p(nxp,j,k); pe=p(nxp+1,j,k);
   do i=nxp,2,-1; im=i-1;
    pw=p(im,j,k)
    po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,jm,k) + n(im,jm,km)*p(i,j+1,k) +&
                               f(im,jm,km)*p(i,j,km) + b(im,jm,km)*p(i,j,k+1) +&
                               po - rhs(im,jm,km) )
    p(i,j,k)=po;
    pe=po; po=pw;
   enddo
  enddo
 enddo ! k
!
 endif
! may be better off unrolling, or reusing solution values already fetched from memory (e.g. pcen=p(i,j,k); pright=p(i+1,j,k) ... pleft=pcen, pcen=pright)
 end subroutine GaussSeidelReverse
!*****************************************************************************************!
 subroutine GaussSeidelRB(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,q2d)
! note that nx,ny,nz include ghost points
! note that rhs passed into this routine is normalized by o (see below)
 implicit none
 integer :: nx, ny, nz
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b ! note that the metrics e--b are variable
! ! (because they include eps & kappa; oo=1/o)
 integer :: nnx, nny, nnz
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km
 logical :: q2d
 real*8 :: pw, pe
 integer :: offset
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 if (q2d) then
!
! red: sum of incides is odd; black: sum of indices is even
! simplest code: perform two sequential loops
! black:
  k=2; km=k-1;
!
  do j=2,nyp; jm=j-1; offset=mod(j,2); pw=p(1+offset,j,k);
   do i=2+offset,nxp,2; im=i-1;
    pe=p(i+1,j,k);
    p(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,jm,k) + n(im,jm,km)*p(i,j+1,k) +&
                               p(i,j,k) - rhs(im,jm,km) )
    pw=pe
  enddo ; enddo ;
!
! red:
  do j=2,nyp; jm=j-1; offset=mod(j,2); pw=p(2-offset,j,k);
   do i=3-offset,nxp,2; im=i-1;
    pe=p(i+1,j,k);
    p(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,jm,k) + n(im,jm,km)*p(i,j+1,k) +&
                               p(i,j,k) - rhs(im,jm,km) )
    pw=pe
  enddo ; enddo ;
!
!
!
 else ! q2d
!
 do k=2,nzp; km=k-1;
  do j=2,nyp; jm=j-1;
   do i=2,nxp; im=i-1;
    p(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
                               s(im,jm,km)*p(i,jm,k) + n(im,jm,km)*p(i,j+1,k) +&
                               f(im,jm,km)*p(i,j,km) + b(im,jm,km)*p(i,j,k+1) +&
                               p(i,j,k) - rhs(im,jm,km) )
   enddo ; enddo ; enddo ! k
!
 endif
!
 end subroutine GaussSeidelRB
!*****************************************************************************************!
 subroutine GaussSeidelRBTiled(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,tilesize,q2d)
! first conceptual GPU template
 implicit none
 integer :: nx, ny, nz
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b ! note that the metrics e--b are variable
! ! (because they include eps & kappa; oo=1/o)
 integer :: nnx, nny, nnz, ntot
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km, ii, jj, kk
 logical :: q2d
 integer :: offset
!
 integer :: tilesize ! conceptually, might be related to CPU cache size
 integer :: tilex, tx, ty, tz, ntx, nty, ntz, remainder ! sizes do not include ghost cells
 real*8, allocatable, dimension(:,:,:) :: e1, w1, s1, n1, f1, b1, p1, rhs1
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 ntot=nnx*nny*nnz
!
 tilesize=min(max(tilesize,4),ntot) ! tile should not be larger than number of points
!
 if (q2d) then
! determine tile dimensions (reuse some variables)
  tilex=INT(sqrt(1.0*tilesize));
  offset=NOT(ISHFTC(1,-1)) ! upper bound -- largest representable integer : all ones except 1st bit
! vary tile sizes a little to get a better cover (heuristic)
  do i=0,1; im = tilex + i; ntx = (nnx/im) + min(1,mod(nnx,im)); ! number of tiles in x-direction
   do j=0,2; jm = tilesize / im + j ; nty = (nny/jm) + min(1,mod(nny,jm)); ! number of tiles in y-direction
    remainder = ntx*nty*im*jm - ntot ! number of elements by which the size of the cover exceeds grid size
!aa
! write(0,*) im, jm, ntx, nty, remainder, ntot
    if (offset.gt.remainder) then
     tx=im; ty=jm;
     offset=remainder
    endif
  enddo; enddo ! i, j
! aa
! write(0,*) tx, ty, tx*ty, (nnx/im) + min(1,mod(nnx,im)), (nny/jm) + min(1,mod(nny,jm))
! allocate data, now that tile size is known
  allocate(w1(tx,ty,1),e1(tx,ty,1),s1(tx,ty,1),n1(tx,ty,1),rhs1(tx,ty,1),p1(tx+2,ty+2,1))
! slide tile over all points and perform RB relaxation within each tile
  k=2 ! 2D
  km=k-1
  jj=1;
  do while (jj.le.nny) ; nty=min(ty,nny-jj+1); ii=1
   do while (ii.le.nnx) ; ntx=min(tx,nnx-ii+1); ! number of elements to read
!aa
! write(0,*) ii, jj
! populate tile
    jm=jj;
    do j=1,nty; im=ii; do i=1,ntx
     w1(i,j,1)=w(im,jm,km); e1(i,j,1)=e(im,jm,km); s1(i,j,1)=s(im,jm,km); n1(i,j,1)=n(im,jm,km); rhs1(i,j,1)=rhs(im,jm,km);
     p1(i,j,1)=p(im,jm,k)
     im=im+1; enddo ; ! i
!
     p1(ntx+1,j,1)=p(im,jm,k); im=im+1 ! remaining boundary points
     p1(ntx+2,j,1)=p(im,jm,k)
     jm=jm+1;
    enddo ! j
! remaining boundary points
    do j=nty+1,nty+2; im=ii; do i=1,ntx+2
     p1(i,j,1)=p(im,jm,k); im=im+1; enddo ; jm=jm+1;
    enddo
!************************************************************************************************
! smooth within tile
!************************************************************************************************
! red: sum of indices is odd; black: sum of indices is even
! black:
    do j=2,nty+1; jm=j-1; offset=mod(j,2)
     do i=2+offset,ntx+1,2; im=i-1;
      p1(i,j,1)=p1(i,j,1) - dt * (w1(im,jm,1)*p1(im,j,1) + e1(im,jm,1)*p1(i+1,j,1) +&
                                  s1(im,jm,1)*p1(i,jm,1) + n1(im,jm,1)*p1(i,j+1,1) +&
                                  p1(i,j,1) - rhs1(im,jm,1) )
    enddo ; enddo ;
! red:
    do j=2,nty+1; jm=j-1; offset=mod(j,2)
     do i=3-offset,ntx+1,2; im=i-1;
      p1(i,j,1)=p1(i,j,1) - dt * (w1(im,jm,1)*p1(im,j,1) + e1(im,jm,1)*p1(i+1,j,1) +&
                                  s1(im,jm,1)*p1(i,jm,1) + n1(im,jm,1)*p1(i,j+1,1) +&
                                  p1(i,j,1) - rhs1(im,jm,1) )
    enddo ; enddo ;
!
!************************************************************************************************
! put solution back into global array (inner points only)
    jm=jj+1
    do j=2,nty+1; im=ii+1; do i=2,ntx+1
     p(im,jm,k)=p1(i,j,1); im=im+1; enddo ; jm=jm+1;
    enddo
! shift tile
    ii=ii+tx
   enddo ! ii
   jj=jj+ty
  enddo ! jj
! free tile memory
  deallocate(w1, e1, s1, n1, rhs1, p1)
 else ! q2d
! 3D case not done yet
 endif
!
 end subroutine GaussSeidelRBTiled
!**********************************************************************************************************************!
 subroutine GaussSeidelRBTiledLoMem(p,rhs,eps,kappa,odxcen,odxcor,odycen,odycor,odzcen,odzcor,nx,ny,nz,dt,tilesize,i2d)
! conceptual GPU template
! compute metrics 'on the fly'
 use constants, only : one
 implicit none
 integer :: nx, ny, nz
 integer(kind=KIND('a')) :: i2d
 real*8 :: dt
 real*8 :: odxcen(nx-1), odxcor(nx-2)
 real*8 :: odycen(ny-1), odycor(ny-2)
 real*8 :: odzcen(nz-1-i2d-i2d), odzcor(nz-2-i2d) ! 3D/2D
 real*8, dimension(nx,ny,nz) :: p, eps
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs, kappa
 real*8 :: w, e, s, n, f, b, o
 integer :: nnx, nny, nnz, ntot
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km, ii, jj, kk, jeps
 integer :: offset
!
 integer :: tilesize ! conceptually, might be related to CPU cache size
 integer :: tilex, tx, ty, tz, ntx, nty, ntz, remainder ! sizes do not include ghost cells
 real*8, allocatable, dimension(:,:,:) :: p1, rhs1, eps1, kappa1
 real*8, allocatable, dimension(:) :: odxcen1, odxcor1, odycen1, odycor1, odzcen1, odzcor1 ! local metrics
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 ntot=nnx*nny*nnz
!
 tilesize=min(max(tilesize,4),ntot) ! tile should not be larger than number of points
!
 if (i2d.eq.1) then
! determine tile dimensions (reuse some variables)
  tilex=INT(sqrt(1.0*tilesize));
  offset=NOT(ISHFTC(1,-1)) ! upper bound -- largest representable integer : all ones except 1st bit
! vary tile sizes a little to get a better cover (heuristic)
  do i=0,1; im = tilex + i; ntx = (nnx/im) + min(1,mod(nnx,im)); ! number of tiles in x-direction
   do j=0,2; jm = tilesize / im + j ; nty = (nny/jm) + min(1,mod(nny,jm)); ! number of tiles in y-direction
    remainder = ntx*nty*im*jm - ntot ! number of elements by which the size of the cover exceeds grid size
!aa
! write(0,*) im, jm, ntx, nty, remainder, ntot
    if (offset.gt.remainder) then
     tx=im; ty=jm;
     offset=remainder
    endif
  enddo; enddo ! i, j
!
! allocate data, now that tile size is known (pad epsilon for easy corner interpolation)
  allocate(rhs1(tx,ty,1),eps1(0:tx+1,0:ty+1,1),kappa1(tx,ty,1),p1(tx+2,ty+2,1))
  allocate(odxcen1(tx+1),odxcor1(tx),odycen1(ty+1),odycor1(ty))
! slide tile over all points and perform RB relaxation within each tile
  k=2 ! 2D
  km=k-1
  jj=1;
  do while (jj.le.nny) ; nty=min(ty,nny-jj+1); ii=1
! local y-metrics
   odycen1(1:nty+1)=odycen(jj:jj+nty); odycor1(1:nty)=odycor(jj:jj+nty-1); ! center-to-center metric has extra element
   do while (ii.le.nnx) ; ntx=min(tx,nnx-ii+1); ! number of elements to read
! local x-metrics
    odxcen1(1:ntx+1)=odxcen(ii:ii+ntx); odxcor1(1:ntx)=odxcor(ii:ii+ntx-1);
!
! populate tile
!
! epsilon: top ghost row (j=0)
    jm=jj-1; jeps=max(jm,1)
    j=0; im=ii;
    do i=1,ntx+1; eps1(i,j,1)=eps(min(im,nnx),jeps,km); im=im+1; enddo ; ! ignore left corner ghost cell
    jm=jm+1;
! inner point j-loop
    do j=1,nty; im=ii-1;
! epsilon: left ghost point
     eps1(0,j,1)=eps(max(1,im),jm,km); im=im+1;
! inner point i-loop
     do i=1,ntx
      rhs1(i,j,1)=rhs(im,jm,km); eps1(i,j,1)=eps(im,jm,km); kappa1(i,j,1)=kappa(im,jm,km); p1(i,j,1)=p(im,jm,k)
      im=im+1;
     enddo ; ! i
! epsilon: right ghost point
     eps1(ntx+1,j,1)=eps(min(im,nnx),jm,km) !
! p: last inner point and right ghost point
     p1(ntx+1,j,1)=p(im,jm,k); im=im+1 ! remaining boundary points
     p1(ntx+2,j,1)=p(im,jm,k) !
     jm=jm+1;
    enddo ! j
! p: bottom inner points
! epsilon: bottom ghost points
    j=nty+1; im=ii;
    jeps=min(jm,nny);
    do i=1,ntx+1
     eps1(i,j,1)=eps(min(im,nnx),jeps,km); ! epsilon : neumann BC
     p1(i,j,1)=p(im,jm,k); im=im+1;
    enddo ;
    p1(ntx+2,j,1)=p(im,jm,k);
    jm=jm+1;
! p: bottom ghost points
    j=nty+2; im=ii;
    do i=1,ntx+1 ! ignore corner ghost point
     p1(i,j,1)=p(im,jm,k); im=im+1;
    enddo ;
!************************************************************************************************
! smooth within tile (iteration involves local data only)
!************************************************************************************************
! red: sum of indices is odd; black: sum of indices is even
! black:
    do j=2,nty+1; jm=j-1; offset=mod(j,2)
     do i=2+offset,ntx+1,2; im=i-1;
! compute metrics
      w = odxcor1(im)*odxcen1(im)*(eps1(im,jm,1)+eps1(im-1,jm,1)) ! 0.5 included using o normalization
      e = odxcor1(im)*odxcen1(i )*(eps1(im,jm,1)+eps1(i ,jm,1))
      s = odycor1(jm)*odycen1(jm)*(eps1(im,jm,1)+eps1(im,jm-1,1))
      n = odycor1(jm)*odycen1(j )*(eps1(im,jm,1)+eps1(im ,j ,1))
      o = -(w+e+s+n) + kappa(im,jm,1) ; o=one/o
!
      p1(i,j,1)=p1(i,j,1) - dt * ( p1(i,j,1) &
& + o*( w*p1(im,j,1) + e*p1(i+1,j,1) + s*p1(i,jm,1) + n*p1(i,j+1,1) ) &
& - rhs1(im,jm,1) )
! NOTE: rhs above has been divided by o, which means that o is already known;
! this routine is _just a test of computing metrics on-the-fly vs. memory access_
    enddo ; enddo ;
! red:
    do j=2,nty+1; jm=j-1; offset=mod(j,2)
     do i=3-offset,ntx+1,2; im=i-1;
! compute metrics
      w = odxcor1(im)*odxcen1(im)*(eps1(im,jm,1)+eps1(im-1,jm,1)) ! 0.5 included using o normalization
      e = odxcor1(im)*odxcen1(i )*(eps1(im,jm,1)+eps1(i ,jm,1))
      s = odycor1(jm)*odycen1(jm)*(eps1(im,jm,1)+eps1(im,jm-1,1))
      n = odycor1(jm)*odycen1(j )*(eps1(im,jm,1)+eps1(im ,j ,1))
      o = -(w+e+s+n) + kappa(im,jm,1) ; o=one/o
!
      p1(i,j,1)=p1(i,j,1) - dt * ( p1(i,j,1) &
& + o*( w*p1(im,j,1) + e*p1(i+1,j,1) + s*p1(i,jm,1) + n*p1(i,j+1,1) ) &
& - rhs1(im,jm,1) )
    enddo ; enddo ;
!
!************************************************************************************************
! put solution back into global array (inner points only)
    jm=jj+1
    do j=2,nty+1; im=ii+1; do i=2,ntx+1
     p(im,jm,k)=p1(i,j,1); im=im+1; enddo ; jm=jm+1;
    enddo
! shift tile
    ii=ii+tx
   enddo ! ii
   jj=jj+ty
  enddo ! jj
! free tile memory
  deallocate(eps1, kappa1, rhs1, p1)
  deallocate(odxcen1, odxcor1, odycen1, odycor1)
!
 else ! i2d
! 3D case not done yet
 endif
!
 end subroutine GaussSeidelRBTiledLoMem
!*****************************************************************************************!
 subroutine GSInner(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,q2d,iter)
! note that nx,ny,nz include ghost points
! note that rhs passed into this routine is normalized by o (see below)
! combine more than one GS iteration in one subroutine call
! loop over grid points in blocks to maximize cache reuse
! default GS behavior recovered using 1 block and no unrolling
 implicit none
 integer :: nx, ny, nz, iter
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b ! note that the metrics e--b are variable
! ! (because they include eps & kappa; oo=1/o)
 integer :: nnx, nny, nnz
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km, idir, jdir, kdir, imid
 logical :: q2d
 real*8 :: pw, pe, po
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 if (q2d) then
!
!###################################################################################
 do m=0,iter-1 ! inner iteration counter
! note that, if m is larger than half the grid, no iterations will take place
!
! start in the middle of the domain and proceed toward the edges
! imid=ishft(nx,-1) ! nx/2
  k=2; km=1;
! do kdir=0,1; do k = ishft(nz,-1) + kdir, (2+m)*(1-kdir) + (nzp-m) * kdir , ishft(kdir,1)-1 ; km=k-1
  do jdir=0,1;
   do j = ny/2 + jdir, (2+m)*(1-jdir) + (nyp-m)*jdir , 2*jdir-1 ; jm=j-1
! do idir=0,1;do i = ishft(nx,-1) + idir, (2+m)*(1-idir) + (nxp-m) * idir , ishft(idir,1)-1 ; im=i-1
!
! p(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
! s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
! p(i,j,k) - rhs(im,jm,km) )
!
! write out i-loops explicitly
!
! imid=ishft(nx,-1)
   imid=nx/2; po=p(imid,j,k); pe=p(imid+1,j,k);
   do i = imid, 2+m,-1; im=i-1
    pw=p(im,j,k)
    po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                               po - rhs(im,jm,km) )
    p(i,j,k)=po;
    pe=po;
    po=pw;
   enddo
!
                           pw=p(imid,j,k); po=p(imid+1,j,k); im=imid
   do i = imid+1, nxp-m, 1
    pe=p(i+1,j,k)
    po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                               s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                               po - rhs(im,jm,km) )
    p(i,j,k)=po;
    pw=po;
    po=pe;
    im=i
! enddo; enddo;
 enddo; enddo ; enddo ; ! enddo
!
 enddo ! m-loop
!###################################################################################
!
 else ! q2d
!
 do k=2,nzp; km=k-1;
  do j=2,nyp; jm=j-1;
   do i=2,nxp; im=i-1;
    p(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
                               s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                               f(im,jm,km)*p(i,j,km) + b(im,jm,km)*p(i,j,k+1) +&
                               p(i,j,k) - rhs(im,jm,km) )
   enddo ; enddo ; enddo ! k
!
 endif
!
 end subroutine GSInner
!*****************************************************************************************!
 subroutine GSOuter(p,rhs,w,e,s,n,f,b,nx,ny,nz,dt,q2d,iter)
! complement to GSInner: update points near the boundary
! NOTE: the unrolling "acrobatics" make only a little difference
 implicit none
 integer :: nx, ny, nz, iter
 real*8 :: dt
 real*8, dimension(nx,ny,nz) :: p
 real*8, dimension(nx-2,ny-2,nz-2) :: rhs,w,e,s,n,f,b ! note that the metrics e--b are variable (because they include eps & kappa; oo=1/o)
 integer :: nnx, nny, nnz
 integer :: nxp, nyp, nzp
 integer :: i, j, k, m, im, jm, km, idir, jdir, kdir, imid, jmid, il, ir, jstep
 logical :: q2d
 real*8 :: pw, pe, po, por
!
 nnx=nx-2; nny=ny-2; nnz=nz-2;
 nxp=nnx+1; nyp=nny+1; nzp=nnz+1;
!
 if (q2d) then
!
!##########################################################################################################
! start at the edges and proceed toward the middle, filling in missing iterations at different levels
  k=2 ; km=k-1
!
  do m=0, iter-1
   il=2+m; ir=nxp-m;
   do jdir=0,1; jmid=ny/2 + jdir; jstep=2*jdir-1 ; po=p(il,jmid,k); por=p(ir,jmid,k)
!
    do j = jmid, (3+m)*(1-jdir) + (nyp-m-1)*jdir, jstep ; jm=j-1 ! i-constant sides
!
     i=il; im=i-1;
     po=po - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
                    s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                    po - rhs(im,jm,km) )
     p(i,j,k)=po
     po=p(i,j+jstep,k)
!
     i=ir; im=i-1;
     por=por - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
                      s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                      por - rhs(im,jm,km) )
     p(i,j,k)=por
     por=p(i,j+jstep,k)
    enddo ! i-constant sides
!
    j=(2+m)*(1-jdir) + (nyp-m)*jdir; jm=j-1; ! j-constant sides
!
! write out i-loops explicitly
!
    imid=nx/2; po=p(imid,j,k); pe=p(imid+1,j,k);
    do i = imid, 2+m,-1; im=i-1
     pw=p(im,j,k)
     po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                                s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                                po - rhs(im,jm,km) )
     p(i,j,k)=po;
     pe=po;
     po=pw;
    enddo
!
                            pw=p(imid,j,k); po=p(imid+1,j,k); im=imid
    do i = imid+1, nxp-m, 1
     pe=p(i+1,j,k)
     po=po - dt * ( w(im,jm,km)*pw + e(im,jm,km)*pe +&
                                s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
                                po - rhs(im,jm,km) )
     p(i,j,k)=po;
     pw=po;
     po=pe;
     im=i
    enddo;
! i-loop over two directions
! do idir=0,1;
! do i = nx/2 + idir, (2+m)*(1-idir) + (nxp-m) * idir, 2*idir-1 ; im=i-1
!
! p(i,j,k)=p(i,j,k) - dt * ( w(im,jm,km)*p(im,j,k) + e(im,jm,km)*p(i+1,j,k) +&
! s(im,jm,km)*p(i,j-1,k) + n(im,jm,km)*p(i,j+1,k) +&
! p(i,j,k) - rhs(im,jm,km) )
! enddo ! i
! enddo ! i-dir
!
   enddo ! j-dir
  enddo ! m
!###################################################################################
!
 else ! q2d
! 3D case not treated yet
 endif
!
 end subroutine GSOuter
!*****************************************************************************************!
 subroutine GaussSeidel2d(p,rhs,w,e,s,n,nx,ny,dt,iter)
! note that nx,ny,nz include ghost points
! note that rhs passed into this routine is normalized by o (see below)
 implicit none
 integer :: nx, ny, iter
 real*8 :: dt
 real*8, dimension(nx,ny) :: p
 real*8, dimension(nx-2,ny-2) :: rhs,w,e,s,n ! note that the metrics e--o are variable (because they include eps & kappa;
 integer :: nnx, nny ! the metrics and rhs are normalized by o
 integer :: nxp, nyp
 integer :: i, j, m, im, jm
!
 nnx=nx-2; nny=ny-2;
 nxp=nnx+1; nyp=nny+1;
!
 do m=1,iter
!
  do j=2,nyp;
   jm=j-1;
   do i=2,nxp
    im=i-1;
    p(i,j)=p(i,j)-dt * ( w(im,jm)*p(im,j) + e(im,jm)*p(i+1,j) +&
                         s(im,jm)*p(i,j-1) + n(im,jm)*p(i,j+1) +&
                                  p(i,j) - rhs(im,jm) )
   enddo ! i
  enddo ! j
 enddo ! m
! may be better off unrolling, or reusing solution values already fetched from memory (e.g. pcen=p(i,j,k); pright=p(i+1,j,k) ... pleft=pcen, pcen=pright)
 end subroutine GaussSeidel2d
!
!*****************************************************************************************!
 subroutine compute_fd_coef(w,e,s,n,f,b,o,eps,kappa,odxcen,odxcor,odycen,odycor,odzcen,odzcor,nx,ny,nz,i2d)
 implicit none
! nx corresponds to the number of inner points
 integer :: nx, ny, nz
 integer(kind=KIND('a')) :: i2d
 real*8, dimension(nx,ny,nz) :: w,e,s,n,f,b,o,kappa
 real*8 :: eps(nx+2, ny+2, nz+2)
 real*8 :: odxcen(nx+1), odxcor(nx)
 real*8 :: odycen(ny+1), odycor(ny)
 real*8 :: odzcen(nz+1-i2d-i2d), odzcor(nz-i2d) ! 3D/2D
! it may be too costly to have all of the metrics in memory
 integer :: i ,j ,k
 integer :: ip, jp, kp
 real*8 :: epscor(nx+1,ny+1,nz+1-i2d)
 real*8 :: oo(nx,ny,nz)
!
! x-metrics:
! compute epsilon at x cell boundary
! note that I am effectively applying 0 neumann conditions to epsilon; it might be required in the future to customize bc (e.g. for periodicity)
! aa
! write(0,*) '*******'
! write(0,*) eps
! write(0,*) odxcor
! write(0,*) odxcen
! write(0,*) '*******'
! write(0,*) '*******'
! write(0,*) kappa
! write(0,*) odycor
! write(0,*) odycen
! write(0,*) odzcor
! write(0,*) odzcen
! write(0,*) '*******'
!
 epscor(:,:ny,:nz)=0.5d0*(eps(2:,2:ny+1,2:nz+1)+eps(:nx+1,2:ny+1,2:nz+1));
!
 do i=1,nx
  ip=i+1
  w(i,:,:)=odxcor(i)*odxcen(i) *epscor(i, :ny,:nz)
  e(i,:,:)=odxcor(i)*odxcen(ip)*epscor(ip,:ny,:nz)
 enddo
!
! y-metrics:
! compute epsilon at y cell boundary
 epscor(:nx,:,:nz)=0.5d0*(eps(2:nx+1,2:,2:nz+1)+eps(2:nx+1,:ny+1,2:nz+1));
!
 do j=1,ny
  jp=j+1
  s(:,j,:)=odycor(j)*odycen(j) *epscor(:nx,j, :nz)
  n(:,j,:)=odycor(j)*odycen(jp)*epscor(:nx,jp,:nz)
 enddo
!
! z-metrics:
! compute epsilon at z cell boundary
 if (i2d.eq.0) then
  epscor(:nx,:ny,:)=0.5d0*(eps(2:nx+1,2:ny+1,2:)+eps(2:nx+1,2:ny+1,:nz+1));
 else
! nothing to do for 2D
 endif
! a single statement:
! epscor(:nx,:ny,:)=0.5d0*(eps(2:nx+1,2:ny+1,2:nz+2-i2d)+eps(2:nx+1,2:ny+1,1+i2d:nz+1));
!
 f=0d0; b=0d0; ! initialize in case they are never used (2D)
!
 do k=1,nz-i2d ! will not be executed in 2D
  kp=k+1
  f(:,:,k)=odzcor(k)*odzcen(k) *epscor(:nx,:ny,k )
  b(:,:,k)=odzcor(k)*odzcen(kp)*epscor(:nx,:ny,kp)
 enddo
!
 o=-(w+e+s+n+b+f)+kappa ;
! invert metrics
 oo=1d0/o;
 w=w*oo; e=e*oo; s=s*oo; n=n*oo; f=f*oo; b=b*oo
!
! aa
! write(7,*) 'xxxxxxxxxx'
! write(7,*) w
! write(7,*) e
! write(7,*) s
! write(7,*) n
! write(7,*) f
! write(7,*) b
! write(7,*) 'xxxxxxxxxx'
! close(7)
!
 end subroutine compute_fd_coef
!*****************************************************************************************!
