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
!
! code extension that allows users to define communicators manually
! required by string module
!
! NOTE: all routines in this module are meant to be executed on ALL processes
!
!**CHARMM_ONLY**!##IF MULTICOM
!
      module multicom
      use ivector ! container for vector of ints
!
! encapsulated communicator type
!
      type comm_type
       type (int_vector) :: nodes
! had to set types below to int4 to work correctly on 64bit gnu
       integer*4 :: me
       integer*4 :: size
       integer*4 :: comm_id
       integer*4 :: group_id
       integer*4 :: parent_id ! ID of the parent communicator ( N/A to first entry )
      end type comm_type
!
! routines
      public multicom_main ! parser
      public multicom_init ! initialization
      private multicom_add ! add communicator
      private multicom_list ! list communicators
      private multicom_set ! set a communicator to a communicator in this module; depends on parallel.fcm
      public multicom_parinit ! reinitialize parallel run (use after communicators have been changed)
      public multicom_cleanup
      private multicom_barrier ! call MPI_BARRIER on any of the known communicators
      public multicom_permute_string_ranks
!
! variables
!
      type (comm_type), allocatable, private, save :: communicators(:) ! array of communicators
      integer, save, private :: ncomm=0 ! number of defined communicators
      logical, save, private :: multicom_initialized=.false.
      integer, save, private :: first_me=-1, first_size=0
!
! parameters
!
      integer, parameter, private :: max_comm = 2048
!
! subroutines
      contains
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine multicom_main(comlyn, comlen)
! multicom command parser
       use parselist
       use multicom_aux
!
       use parser
       use output
       use mpi
!
       implicit none
!
       CHARACTER*(*) COMLYN
       INTEGER COMLEN
!
! local variables
       character(len=8) :: keyword
       character(len=15) :: whoami
       integer*4 :: ncpu, me
       integer :: bug, newcomm
       integer :: kstep, i, j, k, ind, klen
       integer :: nlocal=-1, nrep=-1
       logical :: qens, qstr
       type (int_vector) :: list
       character(len=200) :: msg___(10)=(/'','','','','','','','','',''/); integer :: i_
!
       data whoami /' MULTICOM_MAIN>'/
!
       call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, bug)
       call MPI_COMM_RANK(MPI_COMM_WORLD, me, bug)
!
       keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if ( keyword(1:4).eq.'INIT'(1:4) ) then
        if (me.eq.0) &
        call warning(whoami, 'REINITIALIZING USER-DEFINED COMMUNICATOR LIST.', 0)
        call multicom_init()
        call multicom_list()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif ( keyword(1:3).eq.'ADD'(1:3) ) then
!
! parse node list
        call trima(comlyn, comlen)
        call ilist_parse(list, comlyn)
!
        do i=1,list%last
         if (list%i(i).lt.0.or.list%i(i).ge.ncpu) then
          call warning(whoami, 'ONE OR MORE NODES OUT OF RANGE. NOTHING DONE.', 0)
          call int_vector_done(list)
          return
         endif
!
        enddo
!
        if (list%last.eq.0) then
         call warning(whoami, ' NO NODES SPECIFIED.', 0)
        else
!
         newcomm=multicom_add(list%i(1:list%last)+1) ! offsetting nodes at 0 and adding 1 here permits users to think of the list as ranks in COMM_WORLD
! ! conceptually, this is incorrect because the indices point into the a node array (indexed from 1)
! ! of the first communicator (on top of WORLD), but I am making this compromise for simplicity
         if (newcomm.lt.0) then
          call warning(whoami, ' COULD NOT ADD COMMUNICATOR.', 0)
         else
          write(msg___,222) whoami, newcomm;
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_),3);enddo;msg___='' ! print if verbosity >= 3
         endif
        endif
        call int_vector_done(list)
 222 FORMAT(/A, ' ADDED COMMUNICATOR ',I4,'.')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif ( keyword(1:4).eq.'BARR'(1:4) ) then
        keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn) ! communicator name, e.g. LOCAL
        call multicom_barrier(keyword)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif ( keyword(1:3).eq.'SET'(1:3) ) then
        keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn) ! communicator name, e.g. LOCAL
! process communicator list (essentially same code as for add; communicator indices start at 1)
! here, inode indexes communicators
        comlyn(comlen+1:)=''; comlyn=adjustl(comlyn); comlen=len_trim(comlyn)
        call ilist_parse(list, comlyn)
! check range
        do i=1, list%last
         if (list%i(i).lt.1.or.list%i(i).gt.ncomm) then
          call warning(whoami, ' ONE OR MORE COMMUNICATORS OUT OF RANGE. NOTHING DONE.', 0)
          call int_vector_done(list)
          return
         endif
        enddo ! check
!
        if (list%last.eq.0) then
         call warning(whoami, ' NO COMMUNICATORS SPECIFIED.', 0)
        else
         call multicom_set(keyword, list%i(1:list%last)) ! pass subset of array -- only the communnicators
! warn if LOCAL modified
         if (len(keyword).ge.5) then
          if (keyword(1:5).eq.'LOCAL'.or.keyword(1:5).eq.'GROUP') then ! keep group for backward compatibility
           if (first_me.eq.0) then
            write(msg___, '(2A)') whoami, &
     & ' WARNING: MPI_COMM_LOCAL CHANGED (REINITIALIZING PARALLEL RUN).'
            do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
           endif
           call multicom_parinit()
          endif ! local
         endif ! len
        endif ! ind
!
        call int_vector_done(list)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif ( keyword(1:4).eq.'LIST'(1:4) ) then
        call multicom_list()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif ( keyword(1:4).eq.'DONE'(1:4) ) then
        if (me.ge.0) &
     & call warning(whoami, 'REMOVING USER-DEFINED COMMUNICATOR LIST.', 0)
        call multicom_cleanup()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! now have to reinitialize communicators; this is not done here
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       elseif ( keyword(1:7).eq.'PARINIT'(1:7) ) then
!ccccccc call parallel initialization routines (use after comm_local has been modified)
        call multicom_parinit()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc in this section, provide "canned" communicator setups for various tasks, ]
!cccc e.g. ENSEMBLE/STRING modules
       elseif (.false. &
     & .or. ( keyword(1:4).eq.'ENSE'(1:4) ) & !!**CHARMM_ONLY**!##ENSEMBLE
     & .or. ( keyword(1:4).eq.'STRI'(1:4) ) & !!**CHARMM_ONLY**!##STRINGM
     & ) then
! setup a standard ensemble/string case
        qens=( keyword(1:4).eq.'ENSE'(1:4) )
        qstr=( keyword(1:4).eq.'STRI'(1:4) )
!
        klen=len(keyword)
        keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
        nlocal=atoi(keyword(1:klen))
        if (nlocal.le.0) then
         write(msg___,*)' INVALID NUMBER OF PROCESSORS PER REPLICA (',keyword,').';call warning(whoami, msg___(1), 0)
        else
         i=index(comlyn(1:comlen), 'BY');
         if (i.gt.0) then
          comlyn=comlyn(2+i:comlen)
          comlen=comlen-i-1
          keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
         else
          keyword=''
         endif
         nrep=atoi(keyword(1:i))
         if (nrep.le.0) then
          write(msg___,*)' INVALID NUMBER OF REPLICAS (',keyword,').';call warning(whoami, msg___(1), 0)
         elseif (ncpu.lt.nrep*nlocal) then
          call warning(whoami, ' REQUESTING TOO MANY PROCESSORS. NOTHING DONE.', 0)
         else
! set up for a (nlocal x nrep) run
! use existing routines
          if (qens) keyword='ENSEMBLE'
          if (qstr) keyword='STRING'
          if (first_me.eq.0) then
           write(msg___, 333) whoami, keyword, whoami, nrep, nlocal
 333 FORMAT(/A, ' WILL SET UP COMMUNICATORS FOR ',A, &
     & /,A,' ON ', I5, ' REPLICA(S) WITH ', I5, &
     & ' PROCESSORS PER REPLICA.')
           do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
          endif ! first_me
!
          call multicom_init() ! (re)initialize
! define nrep communicators (local)
          do i=1, nrep
           ind = multicom_add( (/ ((i-1)*nlocal+j, j=1,nlocal) /) )
          enddo
!
          call multicom_set('LOCAL', (/ (i+1, i=1, nrep)/)) ! set local communicator
! define ENSBL/STRNG
          ind = multicom_add( (/ ((j-1)*nlocal+1, j=1, nrep) /) )
          if (qens) call multicom_set('ENSBL', (/ nrep+2/) )
          if (qstr) call multicom_set('STRNG', (/ nrep+2/) )
!
          call multicom_parinit()
! warn about `unused' processors
          if (first_me.eq.0) then
           if (first_size.gt.nlocal*nrep) then
            write(msg___,444) whoami,first_size-nlocal*nrep, keyword
 444 FORMAT(A, I5, ' PROCESSORS WILL NOT BE USED FOR ', &
     & A,' CALCULATIONS.');
            do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
           endif ! first_size
          endif ! first_me
         endif ! nrep
        endif ! nlocal
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! disabled for CHARMM because it would lead to errors
       elseif ( keyword(1:5).eq.'PARSE'(1:5) ) then
!ccccccc set parser communicator to a local or global
        keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn) ! communicator name, e.g. LOCAL
!
        if (.not.multicom_initialized) then
         call warning(whoami, ' COMMUNICATOR MODULE NOT INITIALIZED. NOTHING DONE.', 0)
        else
!
         keyword=adjustl(keyword); klen=len_trim(keyword)
         if (klen.le.0) then
          call warning(whoami, ' COMMUNICATOR NAME NOT SPECIFIED. NOTHING DONE.', 0)
         else
!
! cycle through valid options
!
          klen=5
          bug=0
          if (keyword(1:5).eq.'WORLD') then
! this should not be necessary since 'GLOBAL' wraps 'WORLD'
           MPI_COMM_PARSER=MPI_COMM_WORLD;
           call MPI_COMM_SIZE(MPI_COMM_PARSER, SIZE_PARSER, bug)
           call MPI_COMM_RANK(MPI_COMM_PARSER, ME_PARSER, bug)
          elseif (keyword(1:5).eq.'LOCAL') then
           MPI_COMM_PARSER=MPI_COMM_LOCAL; SIZE_PARSER=SIZE_LOCAL; ME_PARSER=ME_LOCAL;
          elseif (keyword(1:5).eq.'GROUP') then ! for backward compatibility
           MPI_COMM_PARSER=MPI_COMM_LOCAL; SIZE_PARSER=SIZE_LOCAL; ME_PARSER=ME_LOCAL;
           keyword='GROUP'
          elseif (keyword(1:6).eq.'GLOBAL') then
           MPI_COMM_PARSER=MPI_COMM_GLOBAL; SIZE_PARSER=SIZE_GLOBAL; ME_PARSER=ME_GLOBAL;
           klen=6
          else
           write(msg___,*)keyword(1:klen),' IS NOT A VALID OPTION. NOTHING DONE.';call warning(whoami, msg___(1), 0)
           bug=1
          endif
         endif
        endif
!
        if (first_me.eq.0.and.bug.eq.0) then
         write(msg___, '(2A)') whoami, &
     & ' SETTING PARSER COMMUNICATOR TO "',keyword(1:klen),'".'; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
            write(msg___,*)' UNRECOGNIZED SUBCOMMAND: ',keyword;call warning(whoami, msg___(1), 0)
       endif
!
       end subroutine multicom_main
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine multicom_parinit()
       use mpi
       !**CHARMM_ONLY**! use ewald ! ugly hack to work around module detection bug in setmk.com
!
       implicit none
!
       character(len=18) :: whoami
!
       data whoami /' MULTICOM_PARINIT>'/
!
!
      ! nothing by default
!
       end subroutine multicom_parinit
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine multicom_init()
       use mpi
       implicit none
! locals
       integer*4 :: communicator
       integer :: i, j
       integer :: bug=0
       character(len=15) :: whoami
       logical :: found=.false., mpiinit
!
       data whoami /' MULTICOM_INIT>'/
!
       if (multicom_initialized) then
        call multicom_cleanup()
       endif
! allocate communicator array:
       allocate(communicators(max_comm))
       ncomm=0 ! number of communicators
       do i=1, max_comm
        communicators(i)%me=MPI_UNDEFINED
        communicators(i)%size=1
        communicators(i)%comm_id=MPI_COMM_NULL
        communicators(i)%group_id=MPI_GROUP_NULL
        communicators(i)%parent_id=-1
       enddo
!
       if (ncomm.eq.max_comm) then
        call warning(whoami, ' NUMBER OF ALLOWED COMMUNICATORS EXCEEDED.', 0)
        call multicom_cleanup()
        return
       endif
!
       call mpi_initialized(mpiinit,bug)
       if (.not.mpiinit) call mpi_init(bug)
       communicator=MPI_COMM_WORLD
!
       ncomm=ncomm+1
!
       call MPI_COMM_GROUP(communicator, communicators(ncomm)%group_id, &
     & bug)
       call MPI_GROUP_SIZE(communicators(ncomm)%group_id, i,j)
       call MPI_COMM_RANK(communicator, communicators(ncomm)%me, bug)
       call MPI_COMM_SIZE(communicator, communicators(ncomm)%size, bug)
       communicators(ncomm)%comm_id=communicator
       communicators(ncomm)%parent_id=0
       first_me=communicators(ncomm)%me
       first_size=communicators(ncomm)%size
!
! store procs:
       do i=0,first_size-1
        j=int_vector_uadd(communicators(ncomm)%nodes, i)
       enddo
!
       multicom_initialized=.true.
       end subroutine multicom_init
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       function multicom_add(nodes,comm)
! define a new communicator & return pointer into communicator array
! all communicators are defined from the absolute one defined in multicom_init
! it is currently MPI_COMM_WORLD; for a communicator parent_id=n means that the
! _values_ in the nodes array correspond to _indices_ into the parent's nodes array,
! except when n=0, in thich case the _values_ in the nodes array are actual ranks in the
! absolute communicator; these ranks are also group ranks, which can be used to create
! additional groups and communicators; in the current implementation, parent_id is 0
! for all communicators;
!
       use mpi
!
       implicit none
!
       integer, optional :: comm ! omission implies 1st communicator (WORLD)
       integer :: nodes(:)
! locals
       integer :: comm_id, node_id
       integer :: numnodes, nproc, i, j
       integer :: multicom_add
       integer :: bug=0
       integer*4, allocatable :: list(:)
!
       character(len=14) :: whoami
!
       data whoami /' MULTICOM_ADD>'/
!
       multicom_add=-1
       if (.not.multicom_initialized) call multicom_init()
!
       if (present(comm)) then
        comm_id=comm
       else
        comm_id=1
       endif
!
       if (comm_id.gt.ncomm.or.comm_id.lt.1) then
        call warning(whoami, ' INVALID PARENT COMMUNICATOR.', 0)
        return
       endif
!
! parse processor array
       if (ncomm.eq.max_comm) then
        call warning(whoami, ' NUMBER OF ALLOWED COMMUNICATORS EXCEEDED.', 0)
        return
       endif
!
       ncomm=ncomm+1
!
       do i=1,SIZE(nodes)
! determine node_id relative to first communicator (WORLD)
        j=comm_id
        node_id=nodes(i) ! node relative to the requested comm.
        do while (j.gt.0) ! want to get node index relative to the absolute communicator
         node_id=communicators(j)%nodes%i(node_id) ! node id in parent communicator of j
         j =communicators(j)%parent_id
        enddo
! save absolute node ID
        j=int_vector_uadd(communicators(ncomm)%nodes, node_id) ! relative to 1
       enddo
!
       nproc=communicators(ncomm)%nodes%last
       allocate(list(nproc)); ! stay with 4-bit ints for MPI
       list=communicators(ncomm)%nodes%i(1:nproc)
! define new group
       call MPI_GROUP_INCL(communicators(1)%group_id, &
     & nproc, list, & ! 4-bit integer array
     & communicators(ncomm)%group_id, bug)
!
! check for empty group
       if (communicators(ncomm)%group_id.eq.MPI_GROUP_EMPTY) then
        call warning(whoami, ' NO VALID NODES FOUND. NOTHING DONE.', 0)
        ncomm=ncomm-1
        return
       endif
!
! call MPI_GROUP_SIZE(communicators(ncomm)%group_id, i,j)
! create communicator based on new group; on some nodes this will return MPI_COMM_NULL
! although the corresponding group should never be NULL
       call MPI_COMM_CREATE(communicators(1)%comm_id, &
     & communicators(ncomm)%group_id, &
     & communicators(ncomm)%comm_id, bug)
!
! obtain size + rank from communicator
       if (communicators(ncomm)%comm_id.ne.MPI_COMM_NULL) then
        call MPI_COMM_RANK(communicators(ncomm)%comm_id, &
     & communicators(ncomm)%me, bug)
        call MPI_COMM_SIZE(communicators(ncomm)%comm_id, &
     & communicators(ncomm)%size, bug)
! since we expressed all nodes relative to absolute comm, parent ID is 0
        communicators(ncomm)%parent_id=0
       else
        communicators(ncomm)%size=1 ! for consistency with many communication routines, which bail when size=1
        communicators(ncomm)%parent_id=0 ! parent is 0 -- see above
        communicators(ncomm)%me=MPI_UNDEFINED
       endif
!
       multicom_add=ncomm
!
       if (allocated(list)) deallocate(list)
!
       end function multicom_add
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine multicom_cleanup()
! remove all traces of multicom, but leave MPI_COMM_WORLD intact
       use mpi
!
       implicit none
!
       integer :: i, bug
       integer*4 :: result
!
       if (allocated(communicators)) then
        do i=1,ncomm
         call int_vector_done(communicators(i)%nodes)
! destroy group & communicators
         if (communicators(i)%comm_id.ne.MPI_COMM_NULL) then ! it is an error to use MPI_COMM_NULL in compare
       call MPI_COMM_COMPARE(MPI_COMM_WORLD,communicators(i)%comm_id, &
     & result, bug)
! write(0,*) first_me, communicators(i)%comm_id, i
!
          if (result.ne.MPI_IDENT) &
     & call MPI_COMM_FREE(communicators(i)%comm_id, bug) ! keep world
         endif
!cccccccc
         if (communicators(i)%group_id.ne.MPI_GROUP_NULL) &
     & call MPI_GROUP_FREE(communicators(i)%group_id, bug)
        enddo
        ncomm=0
        deallocate(communicators)
       endif
!
! reset communicators
       call parcomm_reset() ! MPI_COMMs correct, but numnod is not updated
       call multicom_parinit()
!
       multicom_initialized=.false.
!
       end subroutine multicom_cleanup
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine multicom_list()
! list currently defined communicators
! use mpi
       use output
!
       implicit none
!
       integer :: i, bug
       integer*4 :: j,k ! aardvark
       character(len=15) :: whoami
 character(len=200) :: msg___(10)=(/'','','','','','','','','',''/); integer :: i_
!
       data whoami /' MULTICOM_LIST>'/
!
       if (.not.multicom_initialized) then
        call warning(whoami, ' COMMUNICATOR MODULE NOT INITIALIZED.', 0)
        return
       endif
!
       if (ncomm.lt.1) then
        call warning(whoami, ' INTERNAL ERROR: NO COMMUNICATORS DEFINED.', 0)
        return
       endif
!
       if (communicators(1)%me.eq.0) then
         write(msg___, 111) whoami ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___='' ! 0 writes
         do i=1, ncomm
           write(msg___, 112) whoami, i, communicators(i)%parent_id ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
           write(msg___, '(" MULTICOM_LIST> ", 10I5)') &
     & communicators(i)%nodes%i(1:communicators(i)%nodes%last) ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         enddo
       endif
 111 FORMAT (/A,' THE FOLLOWING COMMUNICATORS ARE DEFINED:')
 112 FORMAT (A, ' #', I5, ' NODE LIST (RELATIVE TO #',I5,'):')
!
! aardvaark
! j= MPI_UNDEFINED
! if (MPI_COMM_LOCAL.ne.MPI_COMM_NULL)
! & call MPI_COMM_RANK(MPI_COMM_LOCAL,j,bug)
! write(900+first_me,*) first_me, MPI_COMM_LOCAL, j
! if (SIZE_LOCAL.gt.1) CALL MPI_BARRIER(MPI_COMM_LOCAL,bug)
! aardvaark
       end subroutine multicom_list
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine multicom_barrier(comm_name)
       use multicom_aux
       use parser
       use output
       use mpi
!
       implicit none
!
       character(len=*) :: comm_name
! local vars
       character(len=18) :: whoami
       integer :: c, s, m, l, ierror
 character(len=200) :: msg___(10)=(/'','','','','','','','','',''/); integer :: i_
!
       data whoami /' MULTICOM_BARRIER>'/
!
       if (.not.multicom_initialized) then
        call warning(whoami, ' COMMUNICATOR MODULE NOT INITIALIZED. NOTHING DONE.', 0)
        return
       endif
!
       comm_name=adjustl(comm_name); l=len_trim(comm_name)
       if (l.le.0) then
        call warning(whoami, ' COMMUNICATOR NAME NOT SPECIFIED. NOTHING DONE.', 0)
        return
       endif
!
! cycle through existing communicator list
       l=5 ! length of communicator id string
!
       if (comm_name(1:5).eq.'WORLD') then
! this should not be necessary since 'GLOBAL' wraps 'WORLD'; primarily for debugging
        c=MPI_COMM_WORLD; s=1; m=1;
       elseif (comm_name(1:5).eq.'LOCAL') then
        c=MPI_COMM_LOCAL; s=SIZE_LOCAL; m=ME_LOCAL;
       elseif (comm_name(1:5).eq.'GROUP') then ! for backward compatibility
        c=MPI_COMM_LOCAL; s=SIZE_LOCAL; m=ME_LOCAL;
       elseif (comm_name(1:6).eq.'GLOBAL') then
        c=MPI_COMM_GLOBAL; s=SIZE_GLOBAL; m=ME_GLOBAL;
        l=6
       elseif (comm_name(1:5).eq.'ENSBL') then !!**CHARMM_ONLY**!##ENSEMBLE
        c=MPI_COMM_ENSBL; s=SIZE_ENSBL; m=ME_ENSBL; !!**CHARMM_ONLY**!##ENSEMBLE
       elseif (comm_name(1:5).eq.'STRNG') then !!**CHARMM_ONLY**!##STRINGM
        c=MPI_COMM_STRNG; s=SIZE_STRNG; m=ME_STRNG; !!**CHARMM_ONLY**!##STRINGM
       elseif (comm_name(1:6).eq.'PARSER') then
        c=MPI_COMM_PARSER; s=SIZE_PARSER; m=ME_PARSER;
        l=6
       else
        write(msg___,*)comm_name,' IS NOT IN USE. NOTHING DONE.';call warning(whoami, msg___(1), 0)
        return
       endif
!
! if (first_me.eq.0) then
       if (m.eq.0) then ! all roots print
         write(msg___, '(4A)') whoami, &
     & ' CALLING MPI_BARRIER ON COMMUNICATOR "',comm_name(1:l),'".'; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
       endif
!
       if (c.ne.MPI_COMM_NULL.and.m.ne.MPI_UNDEFINED.and.s.gt.1) &
     & call MPI_BARRIER(c, ierror)
!
       end subroutine multicom_barrier
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine multicom_set(comm_name, commlist)
! set various communicators based on the multicom 'database'
! note: all processors are meant to execute this routine
       use multicom_aux
       use output
       use mpi
!
       implicit none
!
       character(len=*) :: comm_name
       integer :: commlist(:)
! local vars
       integer :: llist
       integer :: comm_id, node_id, parent_id
       character(len=14) :: whoami
       integer :: numnodes, c, m, s
       integer :: itype, bug, lcomm, i, j
       logical :: found=.false., foundg=.false.
!
       integer :: flag(0:first_size-1), flagg(0:first_size-1)
 character(len=200) :: msg___(10)=(/'','','','','','','','','',''/); integer :: i_
!
       data whoami /' MULTICOM_SET>'/
!
       if (.not.multicom_initialized) then
        call warning(whoami, ' COMMUNICATOR MODULE NOT INITIALIZED. NOTHING DONE.', 0)
        return
       endif
!
! first check whether the communicator ids are valid
       llist=SIZE(commlist)
       do i=1, llist
        if (commlist(i).gt.ncomm.or.commlist(i).lt.1) then
         call warning(whoami, ' SOME COMMUNICATOR IDs INVALID. NOTHING DONE.', 0)
         return
        endif
       enddo
!
       flag=0
       flagg=0
!
       lcomm=len(comm_name)
!
! loop over all communicators. Note that no check for overlapping groups is made
! in the case a node belonging to multiple groups,
! the last communicator assignment overwrites previous ones
!
       do j=1, llist
        comm_id=commlist(j)
! define dummy vars to shorten code
! note that any of these might be undefined
        c=communicators(comm_id)%comm_id ! could be MPI_COMM_NULL on this node
        m=communicators(comm_id)%me ! could be MPI_COMM_UNDEFINED
        s=communicators(comm_id)%size
!
! detect if null communicator (unlikely, just bug check)
        found=c.eq.MPI_COMM_NULL
        call MPI_ALLREDUCE(found, foundg, 1, MPI_LOGICAL, MPI_LAND, &
     & communicators(1)%comm_id, bug)
        if (foundg) then
         call warning(whoami, ' UNEXPECTED ERROR: COMMUNICATOR HANDLE NULL ON ALL NODES.', 0)
         return
        endif
!
        numnodes=communicators(comm_id)%nodes%last
!
! write(0,*) communicators(comm_id)%nodes%i(1:numnodes)
        do i=1, numnodes
!
         node_id =communicators(comm_id)%nodes%i(i) ! node id in parent (possibly 0)
         parent_id=communicators(comm_id)%parent_id ! possibly 0
! write(700+ME_GLOBAL,*) parent_id ! aa
!
         do while (parent_id.ne.0)
          node_id =communicators(parent_id)%nodes%i(node_id) ! node_id in parent of parent
          parent_id=communicators(parent_id)%parent_id ! parent of parent (possibly 0)
         enddo
! write(0,*) numnodes, node_id, comm_id, parent_id
!
! set communicator
         if (first_me .eq. node_id) then
! check communicator list
          if (lcomm.ge.5) then
!
           flag(node_id)=node_id+1 ! this flag means that a communicator assignment is successful
!
           if (comm_name(1:5).eq.'LOCAL') then
            MPI_COMM_LOCAL=c; SIZE_LOCAL=s; ME_LOCAL=m;
           elseif (comm_name(1:5).eq.'GROUP') then ! for backward compatibility
            MPI_COMM_LOCAL=c; SIZE_LOCAL=s; ME_LOCAL=m;
!
           elseif (comm_name(1:5).eq.'ENSBL') then !!**CHARMM_ONLY**!##ENSEMBLE
            MPI_COMM_ENSBL=c; SIZE_ENSBL=s; ME_ENSBL=m !!**CHARMM_ONLY**!##ENSEMBLE
!
           elseif (comm_name(1:5).eq.'STRNG') then !!**CHARMM_ONLY**!##STRINGM
            MPI_COMM_STRNG=c; SIZE_STRNG=s; ME_STRNG=m !!**CHARMM_ONLY**!##STRINGM
!
           elseif (comm_name(1:6).eq.'PARSER') then
            MPI_COMM_PARSER=c; SIZE_PARSER=s; ME_PARSER=m
           else
            flag(node_id)=0
           endif
          endif ! lcomm
!
         endif ! first_me
!
        enddo ! i=1, numnodes
       enddo ! over commlist array (j)
!
!
       if (kind(itype).eq.8) then
        itype=MPI_INTEGER8
       else
        itype=MPI_INTEGER
       endif
! combine flagg; have to do some ad-hoc gymnastics because of integer size
       call MPI_ALLREDUCE(flag, flagg, &
     & first_size, itype, MPI_MAX, communicators(1)%comm_id, bug)
! flagg array is positive for nodes on which assignment worked, and 0
! for which it did not; for the latter nodes, explicitly set
! communicator to null, size to 1, and me to MPI_UNDEFINED:
       c=MPI_COMM_NULL; s=1; m=MPI_UNDEFINED
       if (flag(first_me).eq.0) then
         if (comm_name(1:5).eq.'LOCAL') then
          MPI_COMM_LOCAL=c; SIZE_LOCAL=s; ME_LOCAL=m;
         elseif (comm_name(1:5).eq.'GROUP') then ! for compatibility with older scripts
          MPI_COMM_LOCAL=c; SIZE_LOCAL=s; ME_LOCAL=m;
!
         elseif (comm_name(1:5).eq.'ENSBL') then !!**CHARMM_ONLY**!##ENSEMBLE
          MPI_COMM_ENSBL=c; SIZE_ENSBL=s; ME_ENSBL=m !!**CHARMM_ONLY**!##ENSEMBLE
!
         elseif (comm_name(1:5).eq.'STRNG') then !!**CHARMM_ONLY**!##STRINGM
          MPI_COMM_STRNG=c; SIZE_STRNG=s; ME_STRNG=m !!**CHARMM_ONLY**!##STRINGM
         elseif (comm_name(1:6).eq.'PARSER') then
          MPI_COMM_PARSER=c; SIZE_PARSER=s; ME_PARSER=m
         endif
       endif ! lcomm
!
       llist=sum(min(flagg,1)) ! # nodes with successful assignments
       found=(llist.gt.0)
       if (.not.found) then
        write(msg___,*)comm_name,' IS NOT IN USE ON REQUESTED NODES.';call warning(whoami, msg___(1), 0)
       else
        if (first_me.eq.0) then
         write(msg___, '(2A)') whoami, &
     & ' ASSIGNED "'//comm_name(1:5)//                               &
     & '" COMMUNICATOR ON THE FOLLOWING NODES:'; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         write(msg___, '(" MULTICOM_SET> ",20I5)') &
     & PACK(flagg, flagg.gt.0)-1 ! node numbering starts from 0
          do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
         if (llist.lt.first_size) & ! warn that on some nodes, communicator is still unassigned
     & write(msg___, '(2A,I5,A)') whoami, &
     & ' COMMUNICATOR "'//comm_name(1:5)//'" UNASSIGNED ON ',         &
     & first_size-llist, ' NODES.' ; do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
        endif ! first_me
       endif ! .not.found
!
       end subroutine multicom_set
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**CHARMM_ONLY**!##IF STRINGM
       subroutine multicom_permute_string_ranks(ranks)
       use multicom_aux
       use mpi
! assume that the string communicator exists (MPI_COMM_STRNG)
! loop over all communicators & find the string communicator
! create a new communicator using the nodes of the string comm.,
! using WORLD as parent. Delete the existing string communicator
! and set the new communicator to MPI_COMM STRNG.
!
       implicit none
       integer :: ranks(:)
!
       integer :: ranks_world(size(ranks)) ! ranks with respect to WORLD communicator
       integer :: i, j, k, comm, node_id, parent_id
       logical :: found
       integer :: itype, bug
       character(len=31) :: whoami
       data whoami /' MULTICOM_PERMUTE_STRING_RANKS>'/
!
       if (.not.multicom_initialized) then
        call warning(whoami, ' COMMUNICATOR MODULE NOT INITIALIZED. NOTHING DONE.', 0)
        return
       endif
!
       if (kind(comm).eq.8) then
        itype=MPI_INTEGER8
       else
        itype=MPI_INTEGER
       endif
!
       if (ncomm.lt.1) then
        call warning(whoami, ' INTERNAL ERROR: NO COMMUNICATORS DEFINED.', 0)
        return
       endif
!
       comm=0
       if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
        found=.false.
        do while (.not.found.and.comm.lt.ncomm)
         comm=comm+1
         found=(communicators(comm)%comm_id.eq.MPI_COMM_STRNG)
        enddo
        if (.not.found) comm=0
       endif
! reduce
       call MPI_ALLREDUCE(comm, i, 1, itype, MPI_BOR,MPI_COMM_GLOBAL, bug)
       comm=i
       if (comm.eq.0) then
        call warning(whoami, ' INTERNAL ERROR: MPI_COMM_STRNG NOT FOUND.', 0)
        return
       endif
!
! add new communicator
! write(800+ME_GLOBAL,*) comm
! write(800+ME_GLOBAL,*) ranks-1
! write(800+ME_GLOBAL,*)
! & communicators(comm)%nodes%i(1:size(ranks))
!
!
       j=multicom_add(ranks, comm)
! assing MPI_COMM_STRNG to new communicator
       do i=1, communicators(j)%nodes%last
!
         node_id =communicators(j)%nodes%i(i) ! node_id in j's parent
         parent_id=communicators(j)%parent_id ! j's parent
!
         do while (parent_id.ne.0) ! want to get node index relative to the global communicator (WORLD)
          node_id =communicators(parent_id)%nodes%i(node_id) ! node id in parent of parent
          parent_id=communicators(parent_id)%parent_id
         enddo
! set communicator
         if (first_me .eq. node_id ) then
          MPI_COMM_STRNG=communicators(j)%comm_id
          SIZE_STRNG= communicators(j)%size
          ME_STRNG= communicators(j)%me
         endif ! first_me
!
       enddo ! i=1, numnodes
! destroy old communicator and group
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! do not do any checks -- assume rest of code works correctly
       if (communicators(comm)%comm_id.ne.MPI_COMM_NULL) & ! comm will be null on slave nodes
     & call MPI_COMM_FREE(communicators(comm)%comm_id, bug)
       call MPI_GROUP_FREE(communicators(comm)%group_id, bug) ! group will be defined on all nodes
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! overwrite old communicator by new communicator;
       communicators(comm)%me =communicators(j)%me
       communicators(comm)%size =communicators(j)%size
       communicators(comm)%comm_id =communicators(j)%comm_id
       communicators(comm)%group_id =communicators(j)%group_id
       communicators(comm)%parent_id=communicators(j)%parent_id
       call int_vector_done(communicators(comm)%nodes)
       do i=1,communicators(j)%nodes%last
        k=int_vector_uadd(communicators(comm)%nodes, &
     & communicators(j)%nodes%i(i))
       enddo
! remove old string communicator
       call int_vector_done(communicators(j)%nodes)
       ncomm=ncomm-1
! write(800+ME_GLOBAL,*)
! & communicators(comm)%nodes%i(1:size(ranks)) !aa
! close(800+ME_GLOBAL)
!
       end subroutine multicom_permute_string_ranks
!**CHARMM_ONLY**!##ENDIF
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module multicom
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**CHARMM_ONLY**!##ENDIF
