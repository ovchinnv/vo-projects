      module eddie
      private
!
!    include other modules
      use parser    ! parses input unit and stores input parameters
      use size      ! basic size information
      use grid      ! grid module
      use fout      ! runtime output
      use IO        ! input/output module
#ifdef AB
      use ABintegrator    ! Adams-Bashforth integrator
#endif
#ifdef RK3
      use RK3integrator    ! Runge-Kutta 3 integrator (not implemented)
#endif
#ifdef central_fd2
       use central_fd2    ! central finite difference, second order derivatives
#endif
!
!
#ifdef immbnd 
      use immbnd      ! use immersed boundary approach to model internal surfaces
#endif
#ifdef two_pt_corr
      use two_pt_corr ! compute two-point correlations on the fly
#endif
#ifdef rogallo
      use rogallo ! rogallo-type turbulent inflow
#endif
#ifdef save_planes
      use save_planes ! save planes of velocity and pressure
#endif
#ifdef precursor_inflow
      use precursor_inflow ! save planes of velocity and pressure
#endif
#ifdef precursor_output
      use precursor_output ! save planes of velocity and pressure
#endif
#ifdef LES
      use LES ! routines/variables for Large-Eddy simulation
#endif
#ifdef stat
      use stat ! routines/variables for computing statistics (on-the-fly)
#endif
#ifdef probe
      use probe ! routines/variables for saving time-traces of velocity and pressure
#endif
!
!     routines:
      public eddie_init
      public eddie_done
      public eddie
!
!     variables declaration section:
!
      eddie_initialized=.false.
!
!
!.....uo, vo, wo:		inst. velocity field
!.....ustar,vstar,wstar:	predicted velocities
!.....p, dp, divergence:        pressure, correction, divergence
!.....uab1, vab1, wab1:		RHS predictor step for Adams-Bashforth advancement
!.....uab2, vab2, wab2:		RHS corrector step (A/B)
!.....tv:                       viscosity
!
      float, allocatable, save, dimension(:,:,:) :: 
     &  uo, vo, wo, ustar, vstar, wstar, p, dp, divergence, tv
!
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine eddie_init(comm)
       implicit none
!
#ifdef PARALLEL
       include 'mpif.h'
#endif
       integer*4, optional :: comm
       integer*4 :: numarg 
       character*(80) :: filename, fname, keyword
       int :: flen
       int :: fid=1 
       logical :: qcomm
       int :: bug
!
       character*(10) :: whoami
!
       data whoami /'EDDIE_INIT'/
!
       if (eddie_initialized) call eddie_done()
!
       qcomm=present(comm)
!
#ifdef PARALLEL
       if (qcomm) call MPI_RANK(comm, me, bug)
#endif
!
       if (me.le.0) then 
        numarg=iargc() ! number of arguments
        if (numarg.ge.1) then 
         call getarg(1, filename)
         fname=filename
         call adjustleft(fname)
         flen=len_trim(fname)
         open(fid, file=filename(1:flen), status='OLD', form='FORMATTED')
        else
         fid=5 ! if file name missing, read from standard input
        endif
       endif ! me
!    call parser
       call parse_file(fid) ! parser will store commands internally -- we can now query parser for options using command()
!
       if (me.le.0.and.fid.ne.5) close(fid)
! *********************************************************************************************************
! * INITIALIZATION
! *
!    1) size
       if (qcomm) then 
        call size_initialize( 
     &   atoi(getval('NX')), atoi(getval('NY')), atoi(getval('NZ')), comm)
       else
        call size_initialize( 
     &   atoi(getval('NX')), atoi(getval('NY')), atoi(getval('NZ')))
       endif
!
       if (.not.size_initialized) then 
        call error(whoami, ' SIZE NOT INITIALIZED. NOTHING DONE.',-1)
        return
       endif
!    set boundary conditions
       call boundary_initialize()
!    prepare Cartesian grid and metrics
       call grid_initialize()
!    allocate basic arrays
       allocate(
     &  uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     &  ustar(nx,ny,nz),vstar(nx,ny,nz),wstar(nx,ny,nz),
     &  p(nx,ny,nz),dp(nx,ny,nz),divergence(nx,ny,nz),
     &  tv(nx,ny,nz)
#ifdef AdamsBashforth
     &  uab1(nx,ny,nz),vab1(nx,ny,nz),wab1(nx,ny,nz),
     &  uab2(nx,ny,nz),vab2(nx,ny,nz),wab2(nx,ny,nz)
#endif
!    initialize LES
#ifdef LES
       call LES_initialize()
#endif LES
!    initialize immersed boundary
#ifdef immbnd 
       call immbnd_initialize()
#endif
!    set initial condition
       call initial()
!------------------------------------------------------------------------
!     Check Divergence of Initial condition
!-----------------------------------------------------------------------
       call CALCDIV(UO,VO,WO)
!
!   other initialization
!    random number generator (might change to a better one!)
       keyword=getval('ISEED')
       if (.not.isempty(keyword)) then 
        iseed=min(-atoi(keyword),atoi(keyword))
       else
        iseed=-1
       endif
       ran1(iseed) ! initialize RNG
!
      
!...  set extra source terms (e.g. pressure gradient, or gravity)
       if (existtag('DPDX')) then 
        dpdx=atof(getval('DPDX'))
       else
        dpdx=0d0
       endif
!
       if (existtag('DPDY')) then 
        dpdy=atof(getval('DPDY'))
       else
        dpdy=0d0
       endif
!
       if (existtag('DPDZ')) then 
        dpdz=atof(getval('DPDZ'))
       else
        dpdz=0d0
       endif
!***************************************************************
       if (dpdx.ne.0) then 
        if (me.eq.0) then 
          if (bctype(west).ne.500.or.bctype(east).ne.500)
     & call warning(whoami,' DPDX =/= 0 FOR NONPERIODIC X-BC',0)
        endif
       elseif (dpdx.ne.0) then
        if (me.eq.0) then 
          if (bctype(front).ne.500.or.bctype(back).ne.500)
     & call warning(whoami,' DPDY =/= 0 FOR NONPERIODIC Y-BC',0)
        endif
       elseif (dpdz.ne.0) then
        if (me.eq.0) then 
          if (bctype(bottom).ne.500)
     & call warning(whoami,' DPDZ =/= 0 FOR NONPERIODIC Z0-BC',0)
        elseif (me.eq.ncpu-1) then
          if (bctype(top).ne.500)
     & call warning(whoami,' DPDZ =/= 0 FOR NONPERIODIC Z1-BC',0)
        endif
       endif
!*****************************************************************
!
#ifdef rogallo
       call rogallo_initialize() ! synthetic turbulence
#endif
#ifdef precursor_output
       call precursor_output_initialize()
#endif
!#ifdef precursor_inflow
!       call precursor_inflow_initialize() ! should move to bc?
!#endif
#ifdef stat
       call stat_initialize() 
#endif
#ifdef probe
       call probe_initialize() ! set velocity/pressure 'probes'
#endif
#ifdef two_pt_corr
       call two_pt_corr_initialize() 
#endif
! other modules to be added here
!
       eddie_initialized=.true.
!
       end subroutine eddie_init
!******************************************************************8
!
      
       subroutine eddie()
       implicit none
!
#ifdef PARALLEL
       include 'mpif.h'
#endif
!

!     headers will disappear
      include 'headers/common.h'
!
!
!    
!.....arrays for the solver
!    
!.....dummy:                    dummy array for Poisson solver

!.....uoc,voc,woc:		velocities at cell centers
!.....upr, vpr, wpr:		fluctuating components
!.....StrkPos,StrfIndx:         streaks
!    
!     

       character*(5) :: whoami
       data whoami /'EDDIE'/
!
!    check for initialization
!
       if (.not.eddie_initialized) then 
        call error(whoami, ' NOT INITIALIZED. NOTHING DONE.',-1)
        return
       endif
!


      real upr(nx,ny,nz),vpr(nx,ny,nz),wpr(nx,ny,nz)
      real ucn(nx,ny,nz),vcn(nx,ny,nz),wcn(nx,ny,nz)
!      real dummy(nx,ny,nz) !artifice for the doubly periodic poisson solver

      real wfrg(nx,ny,nz) ! probably for the "fringe region" boundary layer sim 
!     
      real u2d(nt,ny,nz),v2d(nt,ny,nz),w2d(nt,ny,nz)
      real u2dnd(ntnd,ny,nz),v2dnd(ntnd,ny,nz),w2dnd(ntnd,ny,nz)
!     
      real urun(nx,nz),vrun(nx,nz),wrun(nx,nz)
      real uinf(nx),duinf(nx)
      integer iuinf
!     


      integer statnum,samples,psample


!     start of inflow        
      CHARACTER*(30) InputFile
!     
      REAL TPRESSURE(3)
      REAL    qo,qold,dpdx,tke,divmax,dpdz
      INTEGER i, j, k, ii, jj, kk, iii, jjj, kkk
!     
      LOGICAL Ltstep            ! .f.-predictor .t.-corrector
      REAL    a1,a2,tmeth
!     
      REAL    delt, yplusmin, erm
      INTEGER icycle,istart,icount1,icount2,icinf,Niter,ibnd
      REAL    dtm1,dtold,adiv,zz
!     
      INTEGER ipart
!...  slice data
      INTEGER ic
      REAL    tlevel,t_prev,t_old
      INTEGER iSliceSize,ksub
!...  series of restart files


      logical  LcreateR

!     
      real   a, b, lamda, pi, u0, w0
!     
      INTEGER IERR
      real*8 cpue,cpub
      real*8 cputurvis
      real*8 cpucorr
      real*8 cpupres,cpupredict,cpurhs
!     
!
      integer type !double or single (for MPI)
!     
!     added for boundary layer simulation VO
      integer iistart,imid,imidf,ifin,nload,n,nn,status(MPI_STATUS_SIZE),counter,cpusave
      integer np,ksave
      real xloc, uscale, puscale, pvscale, pwscale, plevel, ux,vx,scx
      integer pind
      real dummy1,dummy2,dummy3
!--------------statistics variables-------------------------------------
      real uactg
!-----------------------------------------------------------------------
!cccc VARIABLES FOR STREAMWISE SPECTRA
      include 'headers/spec.h'
      integer spec_cpu(nspec), spec_k(nspec)
!cccccccccccccccccccccccccccccccccccccccccccccccccc      
!
!     misc variables
      real topbc(3,ny,nz), botbc(3,ny,nz) !these variables are needed when the code is run implicitly
      integer real1d, real2d, stat2d, cross2d !new MPI types
      real dstar ! boundary layer variable
      real oldtime, znewtime ! variables for keeping track of time per timestep

!      temporary debugging variables
      real divmax_,divmin_,divavg_
      integer idmax, jdmax, kdmax, idmin, jdmin, kdmin
      real rsu, rsv, rsw ! currently exist for backward compatibility with "correct"        
      real umx,vmx,wmx
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     input data and setup
!-----------------------------------------------------------------------
!     
!...  fill index array for various output purposes (now in strings)

!     
!...  input parameters (namelist-statements)
         WRITE(6,*) '*..data read'
         write(6,*) 'ru1=', ru1
!     

!
!     initiaize statistics variables (now in stats)
!
!     set up custom MPI types : will be distributed among modules that use these
      call MPI_TYPE_CONTIGUOUS(nz-2,type,real1d,ierr)
      call MPI_TYPE_CONTIGUOUS((nz-2)*nx,type,stat2d,ierr)
      call MPI_TYPE_COMMIT(real1d,ierr)
      call MPI_TYPE_COMMIT(stat2d,ierr)
      call MPI_TYPE_VECTOR(nz,nx,ny*nx,type,real2d,ierr)
      call MPI_TYPE_VECTOR(nz*ny,1,nx,type,cross2d,ierr)
      call MPI_TYPE_COMMIT(cross2d,ierr)
      
      call MPI_TYPE_COMMIT(real2d,ierr)
!     
!

         write(*,*) 'starting time=',tlevel
         write(*,*) 'Initial field'

       call eddie_info(uo,vo,wo,p,nx,ny,nz) ! initial info
      
!     CALL CALCTKE(tke,uo,vo,wo,urun,vrun,wrun,zw,xu,nx,ny,nz)
!-----------------------------------------------------------------------
!     setup immersed boundary -- now part of IB setup
!-----------------------------------------------------------------------

!................initialize counters...................
!.....istart:	contains number of previous iterations
!.....icycle:	iteration counter
!.....itmax:	number of iterations
!.....itini:	start space-time sampling
!.....infnum:	information counter
!.....ipart:	particle tracking counter

      ICOUNT1=0
      ICOUNT2=0
      ICINF  =0
      ICYCLE =ITINI
      ipart  =0
!-----------------------------------------------------------------------
!     compute coefficients
!-----------------------------------------------------------------------
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '*..compute coefficients'
      ENDIF
!
#ifdef LES
!
      dtm1=1. ! should be ineffectual in this first call , but included to avoid an undefined variable error
      CALL turvis(uo,vo,wo,dxdydz,xu,yv,zw,urun,
     &     sxx,syy,szz,sxy,sxz,syz,g,
     &     uoc,voc,woc,upr,vpr,wpr,gbar,
     &     lm,mm,Ilm,Imm,xcf,ycf,zcf,
     &     wrk1,wrk2,wrk3,wrk4,
     &     tv,dtm1,icycle-itini,nx,ny,nz       )
!
#ifdef immbnd
      if(ibswitch.and.Isgsmodel.ne.0) then
!
             CALL IB_forcing(tv,tvim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
!...  numerator for largangian dynamic model
           if(Isgsmodel.eq.20) then
!
             CALL Forcing(Ilm,Ilmim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
             CALL REFRESHBC(Ilm,nx*ny,(nz-2),MYSIZE,MYRANK,communicator)
           endif
         CALL REFRESHBC(tv, nx*ny,(nz-2),MYSIZE,MYRANK,communicator)
      endif
#endif

#endif
!
#ifndef LES
      tv=0.
#endif      	            

!-----------------------------------------------------------------------
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '*..start flow'
      ENDIF
!-----------------------------------------------------------------------
!     begin main loop
!-----------------------------------------------------------------------
!     
      oldtime=0.
 1000 icycle = icycle + 1
!
      znewtime=mpi_wtime()
      write(666,*) icycle, znewtime-oldtime
      oldtime=znewtime
      close(666)
!
      cpurhs=0
      cpupredict=0
      cpupres=0
      cpucorr=0
      cputurvis=0
!.....time step
      if(Lcfl) then
         CALL TSTEPCEN(XU,YV,ZW,UO,VO,WO,TV,nx,ny,nz,DTM1)
      else
         dtm1=tstep
      endif
!     for phase reshuffling -- new routine
#ifdef rogallo
      if (jdose.or.rogal) then 
       if (myrank.eq.0) then 
        dummy1=mod(uinfty*tlevel+zw(1),zlen)
        dummy2=mod(uinfty*tlevel+dtm1+zw(1),zlen)
!	write(0,'(5 G20.10)') dummy1,dummy2,tlevel,dtm1,zlen
!	if (dummy2.lt.dummy1) call reshuffle(jcounter,uh,vh,wh,ky,kz,seed)
	if (dummy2.lt.dummy1) call reshuffle_old(jcounter,uh,vh,wh,count,is,seed)
       endif !myrank

       call MPI_BCAST(uh(1),2*nmodes,type,0,communicator,ierr)		     	
       call MPI_BCAST(vh(1),2*nmodes,type,0,communicator,ierr)		     	
       call MPI_BCAST(wh(1),2*nmodes,type,0,communicator,ierr)		     	
       call MPI_BCAST(ky(1),nmodes,type,0,communicator,ierr)		     	
       call MPI_BCAST(kz(1),nmodes,type,0,communicator,ierr)		     	
      endif !jdose/rogal
#endif
!
#ifdef precursor_inflow
      if (box2) call box2dt(dtm1,icycle,index,type)
#endif
!      
      TLEVEL=TLEVEL+DTM1
!-----------------------------------------------------------------------
!     predictor-corrector step
!-----------------------------------------------------------------------
      ITSTAMP = ITMAX
      ICOUNT1 = MOD(ICYCLE,100000)
      
      IF ((ICOUNT1.EQ.0).OR.((IFIELD.NE.2).AND.(ICYCLE-ITINI.EQ.1))) THEN
!...  Euler step
         a1   =0.
         a2   =1.
         tmeth=1.
         cpub=mpi_wtime()
         CALL Rhs(uo,vo,wo,ucn,vcn,wcn,tv,wfrg,uab2,vab2,wab2,nx,ny,nz )        
         cpurhs=mpi_wtime()-cpub
      ELSE
!...  Adams-Bashforth step
         a1   =0.5
         a2   =1.5
         tmeth=1.
         cpub=mpi_wtime()
         CALL Rhs(uo,vo,wo,ucn,vcn,wcn,tv,wfrg,uab2,vab2,wab2,nx,ny,nz )   
         cpurhs=mpi_wtime()-cpub
      ENDIF

      
!.....Predictor step + boundary conditions for Ustar
!     
      if(impfla.eq.1.and.icycle-itini.eq.1) then ! : to be moved to setup since this is an INIT call
!...  boundary condition setup for implicit part

         CALL BOUN_IMP(botbc,topbc,nx,ny,nz)
      endif
!
      cpub=mpi_wtime()
      CALL predictor(UO,VO,WO,P,UAB1,VAB1,WAB1,
     &     UAB2,VAB2,WAB2,ucn,vcn,wcn,tv,
     &     USTAR,VSTAR,WSTAR,xu,
     &     dpdx,dpdy,dpdz,dtm1,a1,a2,botbc,topbc,
     &     nx,ny,nz                          )
      Ltstep = .FALSE.
!     don't think this is needed anymore since new IB does _NOT_ use ghost points
!      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbc'
!      CALL REFRESHBC(USTAR,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
!      CALL REFRESHBC(VSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
!      CALL REFRESHBC(WSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
!      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'    
!cccccIMMERSED BOUNDARY
#ifdef immbnd
!....modify viscosity at the immersed boundary points
!    DEAL WITH THIS LATER
!      if(ibswitch.and.ftype.eq.1) then !remove pressure from boundary points
!         do i=1,nbodyu
!            ustar(ibodyu(i),j,kbodyu(i)) = ustar(ibodyu(i),j,kbodyu(i))
!     &           +dpdx*dtm1
!         enddo 
!         do i=1,nbodyv
!            vstar(ibodyv(i),j,kbodyv(i)) = vstar(ibodyv(i),j,kbodyv(i))
!     &           +dpdy*dtm1
!         enddo 
!      endif
!
!....compute forcing in x momentum
      if(ibswitch) then
              CALL IB_forcing(ustar,uim,
     &                    uniml,uiiml,ujiml,untop,utk,ustk1,ustk2,unbot,ubk,usti,utib_ind,utib_c1,utib_c2,ubib_ind,
     &                    ubib_c1,ubib_c2,
     &                    uc4intl) 
!.... compute forcing in y momentum
              CALL IB_forcing(vstar,vim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
!.... compute forcing in z momentum
              CALL IB_forcing(wstar,wim,
     &                    wniml,wiiml,wjiml,wntop,wtk,wstk1,wstk2,wnbot,wbk,wsti,wtib_ind,wtib_c1,wtib_c2,wbib_ind,
     &                    wbib_c1,wbib_c2,
     &                    wc4intl) 
      endif
!
!      me=myrank
!      do j=jy1,jy2+1
!       do i=1,vnbod
!        p (vibod(i),j,vjbod(i))=0.
!        vstar(vibod(i),j,vjbod(i))=0.
!       enddo
!       do i=1,unbod
!        ustar(uibod(i),j,ujbod(i))=0.
!       enddo
!       do i=1,wnbod
!        wstar(wibod(i),j,wjbod(i))=0.
!       enddo
!      enddo
!
#endif  
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbc after predictor (or IB)'
      CALL REFRESHBC(USTAR,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
      CALL REFRESHBC(VSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
      CALL REFRESHBC(WSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'    
!cccccIMMERSED BOUNDARY
!     need to call boundary because refreshbc changes wstar at 1st point!!!
      CALL Boundary(USTAR,VSTAR,WSTAR,UO,VO,WO,
     &     u2d,v2d,w2d,tlevel,t_old,ksub,nt,
     &     urun,wrun,
     &     XU,YV,ZW,nx,ny,nz,Ltstep,uinf,duinf,dtm1,x_ref)
!
#ifdef rogallo
      if (jdose) then
       call jdorrsom(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,eta,detacen,detacor,U_,d2U,
     &                 a_w,b_w,a_u,b_u,ps_w,ps_u)
       if (myrank.eq.0) then 
!     add the perturbation
        do i=2,nx-1 
         do j=2,ny-1
          ustar(i,j,1)=ustar(i,j,1)+wose(j,i)*uscale*puscale    
          vstar(i,j,1)=vstar(i,j,1)+vose(j,i)*uscale*pvscale    
          wstar(i,j,1)=wstar(i,j,1)+uose(j,i)*uscale*pwscale
         enddo
        enddo    
       endif ! myrank
      endif !jdose
!          
      if (rogal) then 
!     generate inflow with fourier wall-normal basis cccccc
!       call rogallo(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,type)
       call rogallo_old(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,type)
       if (myrank.eq.0) then
!     add the perturbation
        do i=2,nx-1 
         do j=2,ny-1
          ustar(i,j,1)=ustar(i,j,1)+uose(j,i)*uinfty*uscale*cutoffu(i)*puscale     
          vstar(i,j,1)=vstar(i,j,1)+vose(j,i)*uinfty*uscale*cutoffw(i)*pvscale   
          wstar(i,j,1)=wstar(i,j,1)+wose(j,i)*uinfty*uscale*cutoffw(i)*pwscale
	 enddo
        enddo    
!	   write(0,*) 'max(u,v,w):',umx,vmx,wmx
       endif ! myrank
      endif !rogal
#endif rogallo
!
#ifdef precursor_inflow
      if (box2) call box2read(ustar,vstar,wstar,a_u,a_w,b_u,b_w,ps_u,ps_w)
#endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
      cpupredict=mpi_wtime()-cpub
!     
!      call correct_mass(USTAR,VSTAR,WSTAR,xu,yv,zw,icycle-itini,nx,ny,nz) ! MasBlas routine
      call correct_mass(USTAR,VSTAR,WSTAR,xu,yv,zw,nx,ny,nz) ! mass-balance routine
!     
      cpub=mpi_wtime()
      call divergence(USTAR,VSTAR,WSTAR,nx,ny,nz,DTM1,dp)
!
      if(IPresMethod .eq. 1) then
         CALL PRESDP2P(dp,dp,P,ZWG,YV,XU,DUMMY,
     &        MYRANK,communicator,TPRESSURE)
      elseif(IPresMethod.eq.11) then 
         call pfft_yz(p,dp,dp,dummy,dummy,au,ap,bp(1),cwg(1),myrank,communicator,type)
      elseif(IPresMethod .eq. 22) then

         CALL presdirect(dp,dp,p,xu,yv,zwg,au,cwg,icycle-itini,
     &        myrank,communicator,tpressure)
      else
         call error(whoami,' UNKNOWN PRESSURE SOLVER. ABORTING.',-1)
         return
      endif
!
      cpupres=mpi_wtime()-cpub
!     
!.....Update velocity values , compute residuals and store RHS
!     
      cpub=mpi_wtime()
      call Correct(USTAR,VSTAR,WSTAR,DP,DTM1,nx,ny,nz,UO,VO,WO)
!cccc pressure recovery
!      me=myrank
!      do i=1,ptrashn
!       do j=jy1,jy2+1
!        p(ptrashi(me,i),j,ptrashj(me,i))=0.
!       enddo
!      enddo

#ifdef immbnd
!      me=myrank
!      do j=jy1,jy2+1
!       do i=1,vnbod
!        p (vibod(i),j,vjbod(i))=4.
!        vo(vibod(i),j,vjbod(i))=0.
!       enddo
!       do i=1,unbod
!        uo(uibod(i),j,ujbod(i))=0.
!       enddo
!       do i=1,wnbod
!        wo(wibod(i),j,wjbod(i))=0.
!       enddo
!      enddo
#endif
!-----------------------------------------------------------------------
!     boundary conditions to the updated values
!-----------------------------------------------------------------------
!     
      Ltstep = .TRUE. 
      CALL Boundary(UO,VO,WO,USTAR,VSTAR,WSTAR,
     &     u2d,v2d,w2d,tlevel,t_old,ksub,nt,
     &     urun,wrun,
     &     XU,YV,ZW,nx,ny,nz,Ltstep,uinf,duinf,dtm1,x_ref)
!     
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbs'
      CALL REFRESHBC(UO,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
      CALL REFRESHBC(VO,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
      CALL REFRESHBC(WO,NX*NY,NZ-2,MYSIZE,MYRANK,communicator)
!
      if (box2) call box2bc(ustar,vstar,wstar,uo,vo,wo)
!     .....divergence
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
      CALL Storen(UAB1,VAB1,WAB1,UAB2,VAB2,WAB2,tmeth,nx,ny,nz)
!     
      cpucorr=mpi_wtime()-cpub
!     ccccccccccccccccccccccccccccccccccc
#if precursor_output
       if (save) call save_slice(uo,vo,wo,cpusave,ksave,xu,yv,tlevel,t_prev,icycle,index)
#endif
!
!-----------------------------------------------------------------------
!     eddy viscosity  and forcing
!-----------------------------------------------------------------------
!.... compute runing average of velocities for inflow bc
!     CALL RunAve(uo,vo,wo,urun,vrun,wrun,xu,yv,zw,
!     &            tlevel,dtm1,icycle,nx,ny,nz       )

!.....turbulent viscosity
      cpub=mpi_wtime()
#ifdef LES
!      write(0,*) myrank,' calling turvis'
      CALL turvis(uo,vo,wo,dxdydz,xu,yv,zw,urun,
     &     sxx,syy,szz,sxy,sxz,syz,g,
     &     uoc,voc,woc,upr,vpr,wpr,gbar,
     &     lm,mm,Ilm,Imm,xcf,ycf,zcf,
     &     wrk1,wrk2,wrk3,wrk4,
     &     tv,dtm1,icycle-itini,nx,ny,nz       )
!      write(0,*) myrank,' back from turvis'
#endif
      cputurvis=mpi_wtime()-cpub
!      
#ifdef immbnd
#ifdef LES
!
!....modify viscosity at the immersed boundary points
      if(ibswitch.and.Isgsmodel.ne.0) then

             CALL Forcing(tv,tvim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
!
!...  numerator and denominator for largangian dynamic model
           if(Isgsmodel.eq.20) then
!
             CALL Forcing(Ilm,Ilmim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
             CALL REFRESHBC(Ilm,nx*ny,(nz-2),MYSIZE,MYRANK,communicator)
           endif
         CALL REFRESHBC(tv, nx*ny,(nz-2),MYSIZE,MYRANK,communicator)
      endif
!      write(0,*) myrank,' back from les-immbnd loop'
#endif 
#endif     
!     output for streamwise spectra
#ifdef probe
      if (spec) then 
        if (mod(icycle,specfreq).eq.0) then 
          j=ny/2  
          call spec_write(uo,vo,wo,p,spec_cpu,spec_k,tlevel,j)
        endif
      endif	
#endif ! probe
!
#ifdef immbnd      
!     output for surface pressure
      if (pres_out) then
       if (mod(icycle, presfreq).eq.0) then
!        write(0,*) 'pressure output' 
        call surface_pressure(p, tlevel, npres, npresl, ipresl, kpresl, npres_array, npres_disp, myrank, communicator)
       endif
      endif
#endif
!     
!---------------------------------------------------------------------
 555  continue
      if (.not.(CalcStat.and.planes)) then
!ccccccccccccccccccccccccccc OUTPUT one xz PLANE cccccccccccccccccccccccccccccccccccccccccccc
         if(mod(ICYCLE-itini-1,infnum).eq.0) then
 1234            j=ny/2
                 write(66+mod(icycle,5),*) myrank,': matlab output'
                 call plane_xz(uo,vo,wo,p,xcf,ycf,zcg,STR5,j,real2d)
                 write(66+mod(icycle,5),*) myrank,': done'
         endif !mod
      endif ! .not.CalcStat
!     
!-----------------------------------------------------------------------
!     create a series of restart files
!-----------------------------------------------------------------------
      
      LcreateR=.true.
      if(MOD(ICYCLE,ievery).EQ.0) then
!.... write restart file
         idir  = 1
	 
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': restart file'
!         call IOVECTOR('res.int.vel_f',uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOSCALAR('res.int.p_f',p,nx,ny,nz,nzg,idir,tlevel)
!         call IOVECTOR('res.int.ab_f',uab1,vab1,wab1,nx,ny,nz,nzg,idir,tlevel)
!cccccccccc for large files
         call ioscalar('res.int.u_f',uo,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.int.v_f',vo,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.int.w_f',wo,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.int.uab_f',uab1,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.int.vab_f',vab1,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.int.wab_f',wab1,nx,ny,nz,nzg,idir,tlevel)
!cccccccccc temporary solution for large files
#ifdef LES
         if (Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then  
           call ioscalar('savt.res',Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar('satv.res',tv,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
	 
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
!     
!     write latest inflow phase information(note: this file overwrites the original)
         if (jdose.or.rogal) then 
!          call waves_write(jcounter,kx,ky,kz,uh,vh,wh)
          call waves_write_old(jcounter,kx,ky,kz,count,is,uh,vh,wh)
!	  
          if (myrank.eq.0) then        
            open(101,file='pscale.dat',form='formatted',status='old')
              write(101,'(3 G25.15)') puscale, pvscale, pwscale
	      write(101,'(I6)') psample
	    close(101)  
	  endif     
 999  continue
         endif
      endif
!     
!-----------------------------------------------------------------------
!     write intermediate restart files
!-----------------------------------------------------------------------
!     
      ICOUNT2 = MOD(ICYCLE,INFNUM)
      IF((ICOUNT2.EQ.0)) THEN
         IF(MYRANK.EQ.0) THEN
!-----------------------------------------------------------------------
!     screen and file information
!-----------------------------------------------------------------------
            write(6,*) '*********************************************'
            write(6,*) '*...Iteration number= ',icycle
            write(6,*) '*...Time Step       = ',dtm1
            write(6,*) '*...Total Time      = ',tlevel
            write(6,*) '*...CPU RHS         = ',cpurhs
            write(6,*) '*...CPU PRED        = ',cpupredict
            write(6,*) '*...CPU PRESS       = ',cpupres
            write(6,*) '*...CPU CORR        = ',cpucorr
            write(6,*) '*...CPU TURVIS      = ',cputurvis
         ENDIF
         CALL CALCDIV(uo,vo,wo,nx,ny,nz)
         CALL CALCME(uo,vo,wo,xu,yv,zw,dtm1,nx,ny,nz)
         CALL SCRUVWP(uo,vo,wo,p,nx,ny,nz)
      ENDIF
!     
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
!     
!.....check counters
!
  456 continue
      if (CalcStat.and.(mod(icycle,ISTATS).eq.1)) then 
         samples =samples+1      
         if (samples.gt.statsample) then ! write new stats files
           statnum=statnum+1
           samples=1
	 endif  	  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 
!     compute uinfinity 
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': compute uinft,uact'
      call uinfinity(wo,xu,zwg,uinft,uinftg,uactg,samples,real1d)      
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'	
!     compute free-stream statistics
      call stats_old(uo,vo,wo,zwg,uav,vav,wav,uuav,vvav,wwav,guav,gvav,gwav,guuav,gvvav,gwwav,samples,statnum,imid,imidf,
     &               index,real1d)
!ccccccccc DYNAMIC NOISE ADJUSTMENT ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        if (.false.) then!(jdose.or.rogal) then 
          if (myrank.eq.mysize-1) then 
!     recomupute scaling parameter
            uosm=sqrt(guuav(pind)-guav(pind)**2)/uinfty
            vosm=sqrt(gvvav(pind)-gvav(pind)**2)/uinfty
            wosm=sqrt(gwwav(pind)-gwav(pind)**2)/uinfty
	       	    
      	    puscale=(puscale*psample+puscale*plevel/(uosm))
      	    pvscale=(pvscale*psample+pvscale*plevel/(vosm))
      	    pwscale=(pwscale*psample+pwscale*plevel/(wosm))
      
            psample=psample+1
            puscale=puscale/psample
            pvscale=pvscale/psample
            pwscale=pwscale/psample
            write(0,'("inlet noise rms: u,v,w & scales ",6 G13.5,I5)') uosm,vosm,wosm,puscale,pvscale,pwscale,psample
          endif
!
          call MPI_BCAST(puscale,1,type,mysize-1,communicator,ierr)
          call MPI_BCAST(pvscale,1,type,mysize-1,communicator,ierr)
          call MPI_BCAST(pwscale,1,type,mysize-1,communicator,ierr)
          call MPI_BCAST(psample,1,MPI_int,mysize-1,communicator,ierr)
        endif   

 222  continue
!     compute turbulent boundary layer statistics
      call stats_new(uo,vo,wo,p,tv,xcf,xu,zcg,zwg,uavg,vavg,wavg,pavg,uuavg,vvavg,wwavg,uwavg,ppavg,
     &                     uuuavg,wwuavg,utaug,tauwg,
     &                     uinft,uinftg,dstg,thetag,uactg,winstg,cfg,
     &                     kpg,ubulkg,samples,statnum,real1d,stat2d,index,icycle)
!
        if (.false.) then!(myrank.eq.mysize-1) then 
            open(22,file='shape'//index(statnum)//'.gnu',status='unknown',form='formatted')
            open(27,file='cf'//index2(icycle)//'.gnu',status='unknown',form='formatted')
            do k=kz1,nzg-1 
               write(22,'(6 G20.12)') (zwg(k)-x_ref)*uactg/ru1,cfg(k),dstg(k),thetag(k),kpg(k),zwg(k)-x_ref
               write(27,'(3 G20.12)') (zwg(k)-x_ref)*uactg/ru1,winstg(k),zwg(k)-x_ref
            enddo !k
            close(22)
      	    close(27)
        endif	!myrank
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccc INSTANTANEOUS PLANES cccccccccccccccccccccccccccc
         if (planes) call write_planes(uo,vo,wo,p,
#ifdef LES
     &                        tv,
#endif
     &  xu,xcf,yv,ycf,zcf,zcg,cpuxy,kxy,type,real2d,cross2d,index,index2,icycle)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccREYNOLDS STRESS BUDGETS -- SEE 'RSBUD1.F cccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (rstat) call rstress(rsblk_1,rsblkg_1,rsblk_2,
#ifdef LES      
     &                           rsblk_3,rsblk_4,rsblkg_3,
#endif     
     &                           uo,vo,wo,p,tv,utaug,tauwg,samples,statnum,rsnames,xcf,zcg,stat2d,index)
!cccccccccccccccccccccccccccccccREYNOLDS STRESS BUDGETS -- SEE 'RSBUD1.F cccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccc two-point spanwise correlations cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (two_pt_corr) then 
	  call span_corr(uo,vo,wo,u2pt,v2pt,w2pt,um2pt,vm2pt,wm2pt,ks(myrank),ke(myrank),samples)
!      normalize & write output file       
          call write_corr(u2pt,v2pt,w2pt,um2pt,vm2pt,wm2pt,nx_act,nz_act,xn,ycf,zn,npt,type)
          if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,' two_pt_corr: done'
         endif ! two-pt-corr	   
!ccccccccccccccccccccccccccccccc two-point spanwise correlations cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 333     continue 
      endif ! stat loop

      close(66+mod(icycle,5))      
      IF (icycle .LT. itmax) GO TO 1000
!     
      Niter = Itmax
      IF(MYRANK.eq.0) THEN
         WRITE(6,*) '*.........................'
         WRITE(6,*) '*..EDDY finished.'
      ENDIF
!     
!     write final restart file
!     
      ifield = 2
      idir   = 1
       if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': final restart'
!      call IOVECTOR('res.fin.vel_f',uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
      call IOSCALAR('res.fin.p_f',p,nx,ny,nz,nzg,idir,tlevel)
!      call IOVECTOR('res.fin.ab_f',uab1,vab1,wab1,nx,ny,nz,nzg,idir,
!     &     tlevel)
!cccccccccc for large files
         call ioscalar('res.fin.u_f',uo,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.fin.v_f',vo,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.fin.w_f',wo,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.fin.uab_f',uab1,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.fin.vab_f',vab1,nx,ny,nz,nzg,idir,tlevel)
         call ioscalar('res.fin.wab_f',wab1,nx,ny,nz,nzg,idir,tlevel)
!cccccccccc temporary solution for large files

#ifdef LES
         if (Isgsmodel.eq.20) then  
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then  
           call ioscalar('savt.res',Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar('satv.res',tv,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
!     write latest inflow phase information(note: this file overwrites the original)
         if (jdose.or.rogal) then 
!          call waves_write(jcounter,kx,ky,kz,uh,vh,wh)
          call waves_write_old(jcounter,kx,ky,kz,count,is,uh,vh,wh)
!
          if (myrank.eq.0) then        
            open(101,file='pscale.dat',form='formatted',status='old')
              write(101,'(3 G25.15)') puscale, pvscale, pwscale
	      write(101,'(I6)') psample
	    close(101)  
	  endif
	 endif     
!
       if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,' done'
      call MPI_TYPE_FREE(real1d,ierr)
      call MPI_TYPE_FREE(real2d,ierr)
      call MPI_TYPE_FREE(stat2d,ierr)
      call MPI_TYPE_FREE(cross2d,ierr)
!     
      if(myrank .eq. 0) then 
         close(67)
	 close(66+mod(icycle,5))
      endif
      if (save) call save_close(cpusave)
      if (spec) call spec_close(spec_cpu)
      if (box2) call box2close()
#ifdef immbnd
      if (pres_out) call pressure_close(myrank)
#endif
!
      return
      end
!-----------------------------------------------------------------------
