! ********************************************************************************!
! Preprocessed by GNU CPP !
! ********************************************************************************!
! Plastic Network Model.
!
! Author: Jingzhi Pu & Paul Maragakis
!
! Date: Fri Oct 31 15:02:20 EST 2003
! Date: Thu Nov 6 15:23:49 EST 2003
! Date: Wed Dec 3, 2003 (increased maxcon to 150)
! Date: Fri Apr 2, 2004 (consistency sweep: fortran indexing)
! Date: June 2007 (interfaced to CHARMM)
! overhauled by VO in Aug. 2013 :! bug fixes, parallelization, generalization to an arbitrary number of networks,
! multiple simultaneous plastic networks (e.g. corrsponding to different proteins/domains); NOTE : parallelization may not help in some cases
!---------------------------------------------------------------
module pnm
!
!**CHARMM_ONLY**!##IF PNM (pnm_main)
 
  use constants
  use constants
!
!







  use ivector
  use rvector



!
  implicit none
  private
!







!
  type enet ! elastic network type
   type(int_vector) :: nodes ! pnm node list (no connectivity)
   type(int_vector) :: bonds ! pnm bond list (interacting pairs in tandem)
   type(real_vector) :: x ! x-coordinate
   type(real_vector) :: y ! y-coordinate
   type(real_vector) :: z ! z-coordinate
   type(real_vector) :: fx ! gradient of energy w.r.t x
   type(real_vector) :: fy ! y
   type(real_vector) :: fz ! z
   type(real_vector) :: econt ! energy decomposition array
   type(real_vector) :: r0 ! equilibrium bond lengths (single list corresponding to bonds)
   real*8 :: k ! force constant
   real*8 :: kb ! bonded force constant
   real*8 :: knb ! nonbonded force constant
   real*8 :: cutoff ! cutoff for network
   real*8 :: emin ! energy value at minimum
   real*8 :: eps ! mixing constant
   logical :: initialized
  end type enet
!
  logical :: qpnm ! pnm active
!
  type(enet), private, save, pointer :: networks(:) ! array of elastic network type
  integer :: num_enm ! number of networks
  integer :: max_num_enm ! max number of networks : allocated at initialization
  logical, private :: initialized=.false. ! whether the module has been initialized (possibly with no networks)
  logical, public :: calc_para=.true. ! whether to calculate forces in parallel
  real*8, private, save, pointer :: pnmene(:,:)=>null() ! matrix of pairwise mixing coefficeints
!
  integer, private, save, pointer :: models(:)=>null() ! index of separate pnm models into the networks array
  logical, private, save, pointer :: qexp(:)=>null() ! array of flags indicating whether the exponential PNM version is used
  real*8, private, save, pointer :: beta(:)=>null() ! array of PNM temperatures (only for the exponential version)
  integer, private :: num_models ! number of pnm models
!

  real*8, pointer,save :: wlapack(:)=>NULL()
  logical :: qdouble, qsingle


  public pnm_main
  public pnm_add
  public pnm_ene

  contains

!===================================================================
   subroutine pnm_main(comlyn,comlen) ! parser
   use cmd, only: maxlinelen; use prm, only : vartaglen; use parser
   use output
!
   character(len=*) :: comlyn
   integer :: comlen
! local variables
   character(len=8) :: keyword
   character(len=len("PNM_MAIN") ),parameter::whoami="PNM_MAIN" ! for maximum compatibility
   integer :: i, j
   real*8 :: t
!
   character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_=1
!
   keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
!====================================================================
   if (( keyword(1:4).eq.'DONE'(1:4) )) then
     call pnm_done()
!====================================================================
   elseif (( keyword(1:4).eq.'INIT'(1:4) )) then
     keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
     i=-1
     j=len(keyword)
     j=min(max(0,j),len(keyword));keyword(j+1:)='';call adjustleft(keyword,(/' ',tab/));j=len_trim(keyword)
     if (j.gt.0) i=atoi(keyword(1:j))
     if (i.le.0) then
      i=2
     write(msg___,600) whoami,i ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
 600 format(/A,' VALID NUMBER OF PNM NETWORKS NOT SPECIFIED. USING ',I6)
     endif ! i<=2
     call pnm_init(i)
!====================================================================
   elseif ( ( keyword(1:3).eq.'ADD'(1:3) )) then
     call pnm_add(comlyn, comlen) ! add new rtmd restraint
!====================================================================
! specify parallel calculation options (lifted from string code)
   elseif ( ( keyword(1:4).eq.'PARA'(1:4) )) then
     keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
     select case(keyword)
       case('YES','ON','TRUE','T','yes','on','true','t')
        keyword='ENABLED '; calc_para=.true.
        write(msg___,7009) whoami, keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
       case('NO','OFF','FALSE','F','no','off','false','f')
        keyword='DISABLED' ; calc_para=.false.
        write(msg___,7009) whoami, keyword ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
       case default
        call warning(whoami, 'UNKNOWN OPTION SPECIFIED', -1)
     end select
 7009 format(' ',A,':', ' PARALLEL COMPUTATION OF FORCES ',A)
!====================================================================
! specify exponential version
   elseif (( keyword(1:3).eq.'EXP'(1:3) )) then
     if (.not.initialized) call pnm_init() ! default initialization
     keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
     select case(keyword)
       case('YES','ON','TRUE','T','yes','on','true','t')
        qexp(num_models)=.true.
        t=atof(get_remove_parameter(COMLYN, 'TEMP', COMLEN), 300d0)
        if (t.le.errtol()) then
         call warning(whoami, ' PNM TEMPERATURE MUST BE GREATER THAN ZERO. ABORT.', -1)
         return
        else
         beta(num_models)=one/(kboltzmann*t)
        endif
        write(msg___,7010) whoami, itoa(num_models), ftoa(t); do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
       case('NO','OFF','FALSE','F','no','off','false','f')
        qexp(num_models)=.false.
        write(msg___,7011) whoami, itoa(num_models) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
       case default
        call warning(whoami, 'NO VALID "EXP" OPTION SPECIFIED', -1)
     end select
 7010 format(' ',A,': MODEL #',A,' WILL USE EXPONENTIAL PNM WITH T=',A,'K.')
 7011 format(' ',A,': MODEL #',A,' WILL USE STANDARD PNM.')
!====================================================================
   elseif ( ( keyword(1:4).eq.'NEWM'(1:4) )) then ! new model
    if (.not.initialized) call pnm_init() ! default initialization
    ! can only add a new model if the present model has at least one network !
    if ( num_enm .ge. models(num_models) ) then
! make sure there will be enough room to add network(s) to the model
     if (num_enm.ge.max_num_enm) then
       call warning(whoami, ' MAXIMUM NUMBER OF ALLOWED NETWORKS EXCEEDED (REINITIALIZE).', -1)
     else
       num_models=num_models+1
       models(num_models)=num_enm+1
       write(msg___,601) whoami, itoa(num_models) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
 601 format(' ',A,': WILL CREATE A NEW PNM MODEL (#',A,').')
     endif ! num_enm > max
    endif ! initialized
!====================================================================
   else
    write(msg___(1),*)'UNRECOGNIZED SUBCOMMAND: ',keyword;call warning(whoami, msg___(1), -1)
   endif
!
   end subroutine pnm_main
!====================================================================
   subroutine pnm_done()
   use output
!====================================================================
   character(len=len("PNM_DONE") ),parameter::whoami="PNM_DONE" ! for maximum compatibility
   integer :: i
   type(enet), pointer :: enm
!
   character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_=1
!
   if (.not.initialized) return
!
   if(associated(pnmene))deallocate(pnmene)
   if(associated(models))deallocate(models)
   if(associated(qexp))deallocate(qexp)
   if(associated(beta))deallocate(beta)

   if(associated(wlapack))deallocate(wlapack)

!
   do i=1, num_enm ! over all networks
    enm=>networks(i)
    call int_vector_done(enm%nodes)
    call int_vector_done(enm%bonds)
    call real_vector_done(enm%x)
    call real_vector_done(enm%y)
    call real_vector_done(enm%z)
    call real_vector_done(enm%fx)
    call real_vector_done(enm%fy)
    call real_vector_done(enm%fz)
    call real_vector_done(enm%r0)
    call real_vector_done(enm%econt)
   enddo
   if(associated(networks))deallocate(networks)
!
   num_enm=0
   num_models=0
   initialized=.false.
!
   write(msg___,99) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
 99 format(' ',A,': PNM IS OFF.')
!
   end subroutine pnm_done
!
!===========================================================================
   subroutine pnm_init(n_)
   use cmd, only: maxlinelen; use prm, only : vartaglen; use parser
   use constants
   use output
!===========================================================================
   character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_=1
   character(len=len("PNM_INIT") ),parameter::whoami="PNM_INIT" ! for maximum compatibility
   integer, optional :: n_
   integer :: i, n
   type(enet), pointer :: enm
!
   if (initialized) call pnm_done()
   if (present(n_)) then ; n=n_ ; else ; n=2; endif
!
   if (n.gt.0) then
     allocate(pnmene(n,n), ) ! allocate mixing coefficients
     allocate(networks(n), ) ! allocate network data storage
     allocate(models(n+1), ) ! allocate models array; this is `overdimensioned` for technical simplicity in pnm_ene
     allocate(qexp(n), ) ! exponential PNM flags
     allocate(beta(n), ) ! temperatures for exponential PNM
     do i=1, n
      enm=>networks(i)
!associate(nodes=>networks(i)%nodes, r0=>networks(i)%r0, xn=>networks(i)%x, yn=>networks(i)%y, zn=>networks(i)%z)
! initialize components
      call int_vector_init(enm%nodes)
      call int_vector_init(enm%bonds)
      call real_vector_init(enm%x)
      call real_vector_init(enm%y)
      call real_vector_init(enm%z)
      call real_vector_init(enm%fx)
      call real_vector_init(enm%fy)
      call real_vector_init(enm%fz)
      call real_vector_init(enm%r0)
      call real_vector_init(enm%econt)
      enm%k=zero
      enm%kb=zero
      enm%knb=zero
      enm%cutoff=zero
      enm%emin=zero
      enm%eps=zero
      enm%initialized=.false.
!end associate
     enddo
     num_enm=0
     max_num_enm=n
     pnmene=zero ! initialize mixing matrix
     models=0 ! initialize to zero
     qexp=.false.! off by default
     num_models=1
     models(num_models)=1 ! index of the first pnm model in the networks array
!
     write(msg___,101) whoami, itoa(n) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
 101 format(' ',A,': INITIALIZED PNM MODULE WITH ',A,' NETWORKS (MAX.)')
     initialized=.true.
    else
     write(msg___(1),*)' PNM MODULE INVOKED WITH AN INVALID NUMBER OF NETWORKS(',itoa(n),').';call warning(whoami, msg___(1), -1)
     return
    endif
!
   end subroutine pnm_init
!===========================================================================
  subroutine pnm_add(comlyn, comlen)
!
! add a single elastic network
!
  use psf
  use constants
  use system, only : r, rcomp, m, bfactor, occupancy
  use constants
  use cmd, only: maxlinelen; use prm, only : vartaglen; use parser
  use output
  use system, only : system_getind
  use psf

  character(len=len("PNM_ADD") ),parameter::whoami="PNM_ADD" ! for maximum compatibility
  integer, pointer, dimension(:) :: islct, jslct, kslct

  integer :: isele, i__, iend; integer, pointer::dmolselect(:);character(LEN=20)::word__
  integer :: natom

  character(len=*) :: comlyn
  integer :: comlen
!
  real*8 :: dist2, dref2, dx, dy, dz
  integer :: i, j, k, ii, jj, kk, l, inet
  logical :: qrem, qcomp
!
  type(enet), pointer :: enm
!
  character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_=1

  natom=psf_natom()

!
  if (.not.initialized) then
   call warning(whoami, ' PNM MODULE NOT INITIALIZED. INITIALIZING WITH DEFAULTS.', -1)
   call pnm_init()
  endif
!
  if (num_enm.ge.max_num_enm) then
    call warning(whoami, ' MAXIMUM NUMBER OF ALLOWED NETWORKS EXCEEDED (REINITIALIZE).', -1)
    return
  endif
!
! add network
  inet=num_enm+1 ! network index
! set network parameters
  enm=>networks(inet);
  enm%k =atof(get_remove_parameter(COMLYN, 'FORC', COMLEN), ONE) ! generic force constant
  enm%kb =atof(get_remove_parameter(COMLYN, 'FCBD', COMLEN), enm%k) ! bonded atoms
  enm%knb =atof(get_remove_parameter(COMLYN, 'FCNB', COMLEN), enm%k) ! non-bonded atoms
  enm%cutoff=atof(get_remove_parameter(COMLYN, 'CUT', COMLEN), TEN) ! cutoff for this
  enm%emin =atof(get_remove_parameter(COMLYN, 'ZERO', COMLEN), ZERO) ! equilibrium energy
  enm%eps =atof(get_remove_parameter(COMLYN, 'PMIX', COMLEN), HALF) ! enm mixing constant
!
  write(msg___,*) 'PNM_ADD :   FORCE   V_0   Cutoff   Epsilon'
  write(msg___(2),*) '------------------------------------------'
  write(msg___(3),10)'           ',enm%k, enm%emin, enm%cutoff, enm%eps
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),2);enddo;msg___=''
!
10 FORMAT(A,4F7.2)
! PNM node selection
  allocate(islct(natom), jslct(natom), kslct(natom))
  islct=0; jslct=0; kslct=0;
          if(associated(dmolselect))deallocate(dmolselect)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___ (1)=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___ (1), 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___ (1)=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___ (1)) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          if(associated(dmolselect))deallocate(dmolselect)
          dmolselect=>system_getind(msg___ (1))
! remove selection string from command line:
          msg___ (1) =comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___ (1) ) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ (1) ! selection has been removed from command
          comlen=len_trim(comlyn)
!#undef dmolselect
  islct(dmolselect)=1
! auxiliary selections to generate exclusions between two sub-selections
  qrem = ( remove_tag(COMLYN,'REMO',COMLEN) .gt. 0)
  if (qrem) then
          if(associated(dmolselect))deallocate(dmolselect)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___ (1)=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___ (1), 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___ (1)=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___ (1)) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          if(associated(dmolselect))deallocate(dmolselect)
          dmolselect=>system_getind(msg___ (1))
! remove selection string from command line:
          msg___ (1) =comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___ (1) ) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ (1) ! selection has been removed from command
          comlen=len_trim(comlyn)
!#undef dmolselect
   jslct(dmolselect)=1
          if(associated(dmolselect))deallocate(dmolselect)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___ (1)=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___ (1), 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___ (1)=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___ (1)) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          if(associated(dmolselect))deallocate(dmolselect)
          dmolselect=>system_getind(msg___ (1))
! remove selection string from command line:
          msg___ (1) =comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___ (1) ) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ (1) ! selection has been removed from command
          comlen=len_trim(comlyn)
!#undef dmolselect
   kslct(dmolselect)=1
  endif
!
! use selection array to add node indices and coordinates
  qcomp = (remove_tag(COMLYN,'COMP',COMLEN) .gt. 0) ! whether to take reference coordinates from comparison set
  if (qcomp) then
   write(msg___,104) whoami, 'COMP'
   do i=1,natom
    if (islct(i) .eq. 1) then
!=============================================
     l=int_vector_add(enm%nodes,i)
     l=real_vector_add(enm%x,rcomp(1,i))
     l=real_vector_add(enm%y,rcomp(2,i))
     l=real_vector_add(enm%z,rcomp(3,i))
     l=real_vector_add(enm%fx,zero); l=real_vector_add(enm%fy,zero); l=real_vector_add(enm%fz,zero); l=real_vector_add(enm%econt,zero);
    endif
   enddo
  else
   write(msg___,104) whoami, 'MAIN'
   do i=1,natom
    if (islct(i) .eq. 1) then
     l=int_vector_add(enm%nodes,i)
     l=real_vector_add(enm%x,r(1,i))
     l=real_vector_add(enm%y,r(2,i))
     l=real_vector_add(enm%z,r(3,i))
     l=real_vector_add(enm%fx,zero); l=real_vector_add(enm%fy,zero); l=real_vector_add(enm%fz,zero); l=real_vector_add(enm%econt,zero)
!=============================================
    endif
   enddo
  endif ! qcomp
 104 format(' ',A,': EQUILIBRIUM GEOMETRY TAKEN FROM ',A,' SET.')
!
  write(msg___(2),105) whoami, itoa(enm%nodes%last)
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___=''
 105 format(' ',A,': ',A,' ENM ATOMS FOUND.')
!
! compute connectivity
  dref2=enm%cutoff**2 ! for faster comparison
! O(N^2) loop to find connectivity
!
  do j=1, enm%nodes%last
    do k=j+1, enm%nodes%last
! check for node connection exclusion;
! fetch atom indices
     jj=enm%nodes%i(j) ! atom index of first node
     kk=enm%nodes%i(k) ! atom index of second node
     if ( .not. ( (jslct(jj).eq.1 .and. kslct(kk).eq.1) .or. (jslct(kk).eq.1 .and. kslct(jj).eq.1) ) ) then ! no exclusion
! distance between points i and j
        dx=enm%x%r(j) - enm%x%r(k); dy=enm%y%r(j) - enm%y%r(k); dz=enm%z%r(j) - enm%z%r(k);
        dist2 = dx**2 + dy**2 + dz**2
        if (dist2.le.dref2) then
! add to the network connectivity list
! index negative if nodes are covalently bonded
! not implemented
! update bond list (note that adding node and _NOT_ atom indices
         l=int_vector_add(enm%bonds,j); l=int_vector_add(enm%bonds, sign(k,kk) );
! equilibrium distance value
         l=real_vector_add(enm%r0,sqrt(dist2));
        endif ! dist <= dref
     endif ! no exclusion
    enddo ! k-nodes
  enddo ! jnodes
!
! set additional interaction terms
  do i=1, num_enm
   pnmene(i,inet) = half * (enm%eps + networks(i)%eps )
   pnmene(inet,i) = pnmene(i, inet)
  enddo
!
  enm%initialized=.true.
  num_enm=inet ! increment counter
! free arrays
  if(associated(ISLCT))deallocate(ISLCT)
  if(associated(JSLCT))deallocate(JSLCT)
  if(associated(KSLCT))deallocate(KSLCT)
  if(associated(dmolselect))deallocate(dmolselect)
!
  end subroutine pnm_add
!----------------------------------------------------------------
  subroutine enm_ene(DEU,X,Y,Z,INET)
  use constants
  use constants
!
! Subroutine that returns the elastic network forces of the network (INET)
  real*8, dimension(*) :: X, Y, Z
  integer :: inet ! network index
  integer :: ibeg, iend
  integer :: i, j, ii, jj, iii, jjj
  real*8 :: deu, dist, dref, ddx, ddy, ddz
  real*8 :: d, kf
  type(enet), pointer :: enm
  real*8, pointer, dimension(:) :: dx, dy, dz, econt ! pointers to gradient arrays
!
  if (inet.gt.size(networks)) return
!
  enm=>networks(inet)
!
  dx=>enm%fx%r ; dy=>enm%fy%r ; dz=>enm%fz%r
  econt=>enm%econt%r ! energy decomposition array
  deu=zero;
  do i=1, enm%nodes%last ; dx(i)=zero; dy(i)=zero ; dz(i)=zero; econt(i)=zero ; enddo
!





  ibeg=1; iend=enm%r0%last
!**CHARMM_ONLY**!##ENDIF

!
  do i=ibeg, iend ! over all bonds on this CPU
   jj=2*i; ii=jj-1 ; ! indices into list of pairs
   ii=enm%bonds%i( ii ) ! index of first node
   jj=enm%bonds%i( jj ) ! index of second node
   if (jj.lt.0) then ; kf=enm%kb ; else ; kf=enm%knb ; endif ! ( negative j-index corresponds to a bonded pair )
   jj=abs(jj);
   iii=enm%nodes%i(ii) ! atom index of first node
   jjj=enm%nodes%i(jj) ! atom index of second node
   ddx=(X(iii)-X(jjj)); ddy=(Y(iii)-Y(jjj)); ddz=(Z(iii)-Z(jjj));
   dist=sqrt( ddx**2 + ddy**2 + ddz**2 ) ! distance between nodes i and j
!
   dref=dist-enm%r0%r(i) ! difference from equilibrium distance
!
  
!
!======= Update total energy ====================
   d=kf * dref**2
   deu=deu+d
!======= Update the energy decomposition array ==
   d=fourth*d; ! split equally between two nodes
   econt(ii)=econt(ii) + d ; econt(jj)=econt(jj) + d ;
!==============================================
! Update forces (DX = Grad_x = -f_x in CHARMM convention)
   d=kf*dref/( max ( abs(dist),errtol() ) ); ! protect from overflow
   dx(ii)=dx(ii) + d * ddx ; dy(ii)=dy(ii) + d * ddy; dz(ii)=dz(ii) + d * ddz
   dx(jj)=dx(jj) - d * ddx ; dy(jj)=dy(jj) - d * ddy; dz(jj)=dz(jj) - d * ddz
!======= forces update =======
  enddo ! i: bondlist
!
!end associate
!
  deu = half * deu + enm%emin ! elastic energy
!
 
!
  end subroutine enm_ene
!-------------------------------------------------------------
  subroutine pnm_ene(EU,X,Y,Z,DX,DY,DZ,QECONT,ECONT,NATOM)
! The plastic network model (PNM) energy routine.
! Original Author: Paul Maragakis
! overhaul by VO 2013
! NOTE : parallelization assumes that all coodinates are known
  use constants
  use constants
  use cmd, only: maxlinelen; use prm, only : vartaglen; use parser
!
!
  character(len=len("PNM_ENE") ),parameter::whoami="PNM_ENE" ! for maximum compatibility
  real*8 :: eu, deu ! total energy, energy from a particular PNM
  integer :: natom
  real*8, dimension(natom) :: X, Y, Z, DX, DY, DZ ! coordinates and forces
  logical :: qecont ! decomposition flag
  real*8 :: econt(*) ! decomposition array
!
  integer :: i, j, k, ii
  integer ::imodel, emodel, num_enm_this ! beginning and ending indices of ENMs in the pnm model, number of ENMs in PNM
  real*8 :: dref
!
  integer :: ierr
!
  real*8, pointer, dimension(:) :: fx, fy, fz, edecomp ! short-hand pointers
! variables for diagonalization
  real*8, pointer, dimension(:,:) :: M, evec ! copy of interaction matrix, eigenvectors
  real*8, pointer, dimension(:) :: eval ! eigenvalues
  real*8 :: deval(num_enm) ! force prefactor for individual enms
!
  character(len=200) :: msg___(21)=(/'','','','','','','','','','','','','','','','','','','','',''/); integer :: i_=1
  type(enet), pointer :: enm
!
  if (num_enm.le.0) return
  if (.not.initialized) then
   call warning(whoami, ' PNM MODULE NOT INITIALIZED. NOTHING DONE.', -1)
   return
  endif
!
  eu=zero
!
! note : num_models is less than or equal to num_enm (equality with one ENM per model)
  do j=1, num_enm ! over all networks (in all models)
!
   call enm_ene(pnmene(j,j),X,Y,Z,j)
!**CHARMM_ONLY**!##ENDIF
!
  enddo ! j: over all networks
! reduce energies in parallel
! compute energy of the plastic network
! diagonalize interaction matrix and take the lowest eigenvalue
! using the default diagonalized in CHARMM, which may only work for symmetric
! matrices; this is OK as long as the interaction matrix is kept symmetric
! for the exponential version of the model, diagonalization is not needed (see below)
!
! now loop over all models and diagonalize
  do k=1, num_models ! over all pnm models
   imodel=models(k) ! network index of first model
   emodel=models(k+1)-1 ; if (emodel.le.0) emodel=num_enm ! network index of second model
   num_enm_this=emodel-imodel+1; ! number of ENMs in this PNM
!
   if (qexp(k)) then ! code for exponential version
!=========== exponential PNM
    deu=zero ; dref=zero
!=========== compute in a numerically stable way for large energies
! 1) find minimum energy
    ii=imodel ! location of minimum
    deu=pnmene(ii,ii)
    do j=imodel+1, emodel
     if (deu.gt.pnmene(j,j)) then ; ii=j ; deu=pnmene(j,j) ; endif
    enddo
! 2) compute correction to minimum energy
    do j=imodel, emodel
      deval(j)=exp( - beta(k) * (pnmene(j,j)-deu) ) ! relative to minimum
      dref = dref + deval(j)
    enddo
    deu = deu - log ( dref ) / beta(k)
    deval(imodel:emodel) = deval(imodel:emodel) / dref ! gradient contribution coefficients for all except ii

   else
!============ standard PNM
    allocate(M(num_enm_this, num_enm_this), evec(num_enm_this, num_enm_this), eval(num_enm_this))
    M=pnmene(imodel:emodel,imodel:emodel); ! copy part of interaction matrix corresponding to this PNM
! diagonalize matrix
! adopted from pca.ftn
    if (.not. associated(wlapack)) then
     qdouble=(kind(eval).eq.kind(1d0));
     qsingle=(kind(eval).eq.kind(1.0));
     if (qdouble) then ! double precision
      call dsyev('V','L', num_enm, eval, num_enm, eval, eval, iminusone, ierr)
     elseif (qsingle) then ! single precision
      call ssyev('V','L', num_enm, eval, num_enm, eval, eval, iminusone, ierr)
     else
      write(msg___(1),*)'Cannot find compatible LAPACK diagonalization routine for kind "',itoa(kind(eval)),'". Trying to abort';call warning(whoami, msg___(1), -1);
      call pnm_done()
      return
     endif
!
     if (ierr.ne.0) then
      write(msg___(1),*)'Error calculating work array size for LAPACK diagonalizer. Trying to abort';call warning(whoami, msg___(1), -1);
      call pnm_done()
      return
     endif
!
     i=nint(eval(1)) ; ! dimension of work array
     allocate(wlapack(i), ) ;
!
    endif ! associated
!
    if (qdouble) then
     call dsyev('V','L', num_enm_this, M, num_enm_this, eval, wlapack, i, ierr)
    elseif (qsingle) then
     call ssyev('V','L', num_enm_this, M, num_enm_this, eval, wlapack, i, ierr)
    else
! should be impossible to get here
    endif
!
    ii=1; deu=eval(ii); do i=2, num_enm_this ; if ( eval(i) .lt. deu ) then ; ii=i ; deu=eval(ii) ; endif ; enddo ! scan all evals to find lowest
! compute corresponding eigenvalue (energy) derivatives w.r.t individual ENM energies (diagonal matrix components);
    deval(imodel:emodel) = evec (:, ii)**2 ;
! make sure eigenvectors are normalized to unity:
    dref = sum(deval(imodel:emodel)) ; if (dref .gt. errtol()) deval(imodel:emodel)=deval(imodel:emodel) / dref ;
! deallocate arrays
    deallocate(M, evec, eval)
   endif ! qexp
   eu=eu+deu ! add this network`s energy contribution to total PNM energy
  enddo ! over all models
!



! apply forces
! all CPUs do this (even in parallel, because in that case only partial forces are computed by each CPU)
!
  if (qecont) then
   do j=1,num_enm ! over all networks
    enm=>networks(j)
    fx=>enm%fx%r; fy=>enm%fy%r; fz=>enm%fz%r; edecomp=>enm%econt%r;
    do i=1, enm%nodes%last
     ii=enm%nodes%i(i) ! index
!
!**CHARMM_ONLY**!##ELSE
     dx(ii) = dx(ii) + deval(j) * fx(i);
     dy(ii) = dy(ii) + deval(j) * fy(i);
     dz(ii) = dz(ii) + deval(j) * fz(i);
     econt(ii) = econt(ii) + edecomp(i);
!**CHARMM_ONLY**!##ENDIF
!
    enddo ! over nodes
   enddo ! over networks
  else
   do j=1,num_enm ! over all networks
    enm=>networks(j)
    fx=>enm%fx%r; fy=>enm%fy%r; fz=>enm%fz%r;
    do i=1, enm%nodes%last
     ii=enm%nodes%i(i)
!
!**CHARMM_ONLY**!##ELSE
     dx(ii) = dx(ii) + deval(j) * fx(i)
     dy(ii) = dy(ii) + deval(j) * fy(i)
     dz(ii) = dz(ii) + deval(j) * fz(i)
!**CHARMM_ONLY**!##ENDIF
!
    enddo ! over nodes
   enddo ! over networks
  endif
!
!
  end subroutine pnm_ene
!
!**CHARMM_ONLY**!##ENDIF (pnm_main)
!======================================================================
end module pnm
