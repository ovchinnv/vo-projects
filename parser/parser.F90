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
      module parser
      implicit none
      private
! read, parse & store input file;
!
! define a derived type to store simulation parameters
      type params
       character(len=200), dimension(:), pointer :: tag, val
       integer, dimension(:), pointer :: tlen, vlen
       integer :: length ! length of the vector
       integer :: last ! index of last element
       logical :: initialized=.false. ! has the vector been initialized
      end type params
!
       character(len=26), parameter :: lower='abcdefghijklmnopqrstuvwxyz'
       character(len=26), parameter :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       character(len=10), parameter :: digits='0123456789'
       integer :: i
       character, parameter :: lower2(26)=(/ (lower(i:i),i=1,26)/)
       character, parameter :: upper2(26)=(/ (upper(i:i),i=1,26)/)
       character, parameter, public :: digits2(10)=(/ (digits(i:i),i=1,10)/)
       character, parameter :: tilde='~'
       character, parameter :: decimal='.'
       character, parameter :: underscore='_'
       character, parameter, public :: hyphen='-'
       character, parameter :: slash='/'
       character, parameter, public :: tab=char(9)
       character, parameter, public :: space(4) = (/' ', tab, ',',';'/)
       character, parameter :: comment(4) = (/'*','!','#','%'/)
       character, parameter :: equals(1) = (/'='/)
       character, parameter :: vector_open(1)=(/'('/)
       character, parameter :: vector_close(1)=(/')'/)
       character(len=200), parameter :: allowed=upper//lower//digits//decimal//underscore//hyphen//slash//tilde
!
      integer, parameter, private :: expand_incr=100
!
      logical, save :: parser_initialized=.false. ! set to true after parse_file is called successfully; private
      type (params), save :: parameters
!
      public atoi ! convert string to integer
      public atof ! convert string to double
      public atol ! convert string to logical
      public atoiv ! convert string to a vector of ints
      public atofv ! convert string to a vector of double
! public atolv ! convert string to a vector of logical (not yet)
      private params_getval ! return tag value
      private params_getval_nocase ! return tag value (case insensitive version)
      public getval ! return tag value
      public getval_nocase ! return tag value (case insensitive version)
      public getval_nocase_upper! return tag value in uppercase (case insensitive version)
      public existtag ! return tag value
      public existtag_nocase ! return tag value (case insensitive version)
      public parse_file ! read input file and store all parameters
      public parser_done
      public list_params ! list parameters
      public adjustleft
      public numword
! CHARMM-compatibility routines
      public pop_string
      public find_tag
      public remove_tag
      public get_remove_parameter
      public toupper
      public tolower
      public itoa
      public ftoa
      public ltoa
!
      contains
!******************************************** implement data routines*********************
       subroutine params_init( v )
       type (params) :: v
! if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%tag(expand_incr),v%val(expand_incr),v%tlen(expand_incr),v%vlen(expand_incr))
!
       v%tag=''; v%val=''; v%tlen=0; v%vlen=0
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine params_init
!ccccc
       subroutine params_done( v )
       type (params) :: v
       if (associated(v%tag)) deallocate(v%tag)
       if (associated(v%val)) deallocate(v%val)
       if (associated(v%tlen)) deallocate(v%tlen)
       if (associated(v%vlen)) deallocate(v%vlen)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine params_done
!
       subroutine params_expand( v )
       type (params) :: v
       integer :: newlength
       character(len=200), dimension(:), pointer :: ntag, nval
       integer, dimension(:), pointer :: ntlen, nvlen
!
       if (.not.v%initialized) then
        call params_init(v)
       else
! assume length is valid
        newlength=v%length+expand_incr ! temporary storage space
        allocate(ntag(newlength),nval(newlength),ntlen(newlength),nvlen(newlength)) ! copy old data
        ntag(1:v%length)=v%tag; nval(1:v%length)=v%val; ntlen(1:v%length)=v%tlen; nvlen(1:v%length)=v%vlen ! deallocate old array
        deallocate(v%tag,v%val,v%tlen,v%vlen) ! deallocate old data
        v%tag=>ntag; v%val=>nval; v%tlen=>ntlen; v%vlen=>nvlen
        v%length=newlength
       endif
       end subroutine params_expand
!ccccc
       function params_add(v,newtag,newval,ltag,lval) ! add a new element to the list (not necessarily unique)
! and return its index
       use output
       type (params) :: v
       integer :: params_add, ltag, lval
       character(len=*) :: newtag, newval
       integer :: j
!
       if (.not.v%initialized) call params_init(v)
! add element to the list
       if (v%last.eq.v%length) call params_expand(v)
       j=v%last+1
       v%tag(j)=newtag(1:ltag); v%tlen(j)=ltag
       v%val(j)=newval(1:lval); v%vlen(j)=lval
       v%last=j
       params_add=j
       end function params_add
!ccccc
       function params_uadd(v,newtag,newval,ltag,lval) ! add a UNIQUE new element to the list and return its index
! if the element already exists, overwrite and warn
       use output, only: warning
       type (params) :: v
       integer :: j, params_uadd, ltag, lval
       character(len=*) :: newtag, newval
!
       if (.not.v%initialized) call params_init(v)
       do j=1,v%last
        if (v%tag(j).eq.newtag(1:ltag)) then
! found element
         params_uadd=j
         call warning('PARAMS_UADD','Parameter "'//newtag(1:ltag)//'" is already present and has the value '&
& //v%val(j)(1:v%vlen(j))//&
& '. Will overwrite.',0)
         v%val(j)=newval(1:lval); v%vlen(j)=lval
         return
        endif
       enddo
! add element to the list: use regular routine
       params_uadd=params_add(v,newtag,newval,ltag,lval)
!
       end function params_uadd
!
       function params_getval( v,atag )
       type (params) :: v
       character(len=200) :: params_getval
       character(len=*) :: atag
       integer :: j
!
       if (.not.parser_initialized) call parser_init()
       params_getval=''
       if (.not.v%initialized) then
!
       else
        do j=1,v%last
         if (v%tag(j).eq.atag) then
! found element
          params_getval=v%val(j)(1:v%vlen(j))
          return
         endif
        enddo
       endif
!
       end function params_getval
!
       function params_getval_nocase( v,atag )
       type (params) :: v
       character(len=200) :: params_getval_nocase, tag1, tag2
       character(len=*) :: atag
       integer :: j
!
       if (.not.parser_initialized) call parser_init()
       params_getval_nocase=''
       if (.not.v%initialized) then
! nothing
       else
        tag2=atag; call toupper(tag2)
        do j=1,v%last
         tag1=v%tag(j);
         call toupper(tag1)
         if (tag1.eq.tag2) then
! found element
          params_getval_nocase=v%val(j)(1:v%vlen(j))
          return
         endif
        enddo
       endif
!
       end function params_getval_nocase
!****************************************** end of data routines **************************************
       character(len=200) function getval(atag)
       use output, only: warning
       character(len=200) :: value
       character(len=*) :: atag
!
       value=params_getval(parameters,atag)
       if (len_trim(value).eq.0) call warning('GETVAL','Parameter "'//trim(atag)//'" not found.',-1)
       getval=value
!
       end function getval
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       character(len=200) function getval_nocase(atag)
       use output, only: warning
       character(len=200) :: value
       character(len=*) :: atag
!
       value=params_getval_nocase(parameters,atag)
       if (len_trim(value).eq.0) call warning('GETVAL_NOCASE','Parameter "'//trim(atag)//'" not found.',-1)
       getval_nocase=value
!
       end function getval_nocase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       character(len=200) function getval_nocase_upper(atag)
       use output, only: warning
       character(len=200) :: value
       character(len=*) :: atag
!
       value=params_getval_nocase(parameters,atag)
       if (len_trim(value).eq.0) call warning('GETVAL_NOCASE_UPPER','Parameter "'//trim(atag)//'" not found.',-1)
       call toupper(value)
       getval_nocase_upper=value
!
       end function getval_nocase_upper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function existtag(atag)
       use output, only: error
       logical :: existtag
       character(len=200) :: value
       character(len=*) :: atag
!
       existtag=.false.
       if (.not.parser_initialized) call parser_init()
       value=params_getval(parameters,atag)
       if (len_trim(value).gt.0) existtag=.true.
!
       end function existtag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function existtag_nocase(atag)
       logical :: existtag_nocase
       character(len=200) :: value
       character(len=*) :: atag
!
       existtag_nocase=.false.
       if (.not.parser_initialized) call parser_init()
       value=params_getval_nocase(parameters,atag)
       if (len_trim(value).gt.0) existtag_nocase=.true.
!
       end function existtag_nocase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine parser_init()
! initialize params structure
       call params_done(parameters)
       call params_init(parameters)
       parser_initialized=.true.
       end subroutine parser_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine parser_done()
       call params_done(parameters)
       parser_initialized=.false.
       end subroutine parser_done
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       subroutine parse_file(fid &
     & )
       use output
!
       implicit none
!
!
       integer :: fid ! input file handle
       integer :: i, j
       character(len=10), parameter :: whoami = 'PARSE_FILE'
!
       integer :: allowed_len, lower_len, upper_len, digits_len
       character :: allowed2(200),c
       logical :: qtag, qeq, qval, qerror
       logical :: qvector ! flag that is active when we are reading a vector (certain default behavior is overridden)
!
       integer :: ioerr, ierr
       integer :: l=0, ltag=0, lval=0
       character(len=400) :: cmdline
       character(len=200) :: tag, val
!
       allowed_len=len_trim(allowed)
       do i=1,allowed_len
        allowed2(i)=allowed(i:i)
       enddo
!
! do work
!
       qvector=.false.
       qerror=.false.
       if (.not. parser_initialized) call parser_init()
       call message(whoami, 'Reading input file.')
       do while (.true.)




        read(fid,'(A)',IOSTAT=ioerr) cmdline ! if running in parallel, then only the root node is passed a valid handle

!
        if (ioerr.eq.0) then



! write(0,*) cmdline
! read from the line
         call adjustleft(cmdline)
         l=len_trim(cmdline)
! write(0,*) cmdline(1:l)
! add comment character to know when to stop below
         if (l.lt.200) l=l+1
         cmdline(l:l)='*'
         if (any(comment.eq.cmdline(1:1))) l=1 ! skip lines that are comments

         i=0
         qtag=.true. ; ltag=0 ! each line is required to begin with a tag
         qval=.false.; lval=0
         qeq=.false.
         qvector=.false. ! this implies that a vector entry cannot span multiple lines
!
         do while (l.gt.1)
           i=i+1
           c=cmdline(i:i)
! write(0,*) tag(1:ltag), val(1:lval),i
           if (any(comment.eq.c.or.i.eq.l)) then ! end of command line
            if ((qtag.and.i.gt.1).or.qeq) then
             call warning(whoami,'Unexpected end of line',0)
             qerror=.true.; exit
            elseif (qval) then
             if (lval.eq.0) then
              call warning(whoami, 'Missing value for parameter "'//tag(1:ltag)//'". Skipping line',0)
              qerror=.true.; exit
            elseif (qvector) then
              call warning(whoami, 'Missing closing bracket for vector parameter "'//tag(1:ltag)//'". Skipping line',0)
              qerror=.true.; exit
             else ! store tag and val
              call message(whoami, tag(1:ltag)//' <= '//val(1:lval))
              j=params_uadd(parameters,tag,val,ltag,lval)
              ltag=0; lval=0
              exit
             endif
            else
             exit
            endif ! qtag
!
           elseif (any(space.eq.c)) then ! completed a tag or value (since lines cannot begin with blanks -- see above)
!
            if (qtag) then
             qtag=.false.
             qeq=.true.
            elseif (qeq) then
! nothing
            elseif (qval) then
             if (qvector) then ! if we are in the process of reading a vector, then the separator (space) tags do not apply;
! otherwise, a tag/value pair has been completed, so attempt to add pair and reinitialize for new tag
              lval=lval+1; val(lval:lval)=' ' ! substitute a space for all separators
             else
!
              qval=.false.
              qtag=.true.
              if (lval.eq.0) then
               call warning(whoami, 'Missing value for parameter "'//tag(1:ltag)//'". Skipping line',0)
               qerror=.true.; exit
              else ! store tag and val
               call message(whoami, tag(1:ltag)//' <= '//val(1:lval))
               j=params_uadd(parameters,tag,val,ltag,lval)
               ltag=0; lval=0
               cmdline=cmdline(i:l) ! remove tag/val pair from string
               call adjustleft(cmdline)
               l=len_trim(cmdline)
               i=0
              endif ! lval
              ltag=0
             endif ! qvector
            endif ! qtag
!
           elseif (any(equals.eq.c)) then ! completed a tag
!
            if (qtag) then
             qtag=.false.
             if (ltag.eq.0) then
              call warning(whoami, 'Missing parameter name. Skipping line',0)
              qerror=.true.; exit
             endif
!
             qval=.true. ; lval=0
             cmdline=cmdline(i+1:l) ! remove 'tag=' and leading spaces pair from string
             call adjustleft(cmdline)
             l=len_trim(cmdline)
             i=0
!
             qeq=.false.
            elseif (qeq) then
             qeq=.false.
             qval=.true. ; lval=0
             cmdline=cmdline(i+1:l) ! remove 'tag=' and leading spaces pair from string
             call adjustleft(cmdline)
             l=len_trim(cmdline)
             i=0
            elseif (qval) then
             if (lval.eq.0) then
              call warning(whoami, 'Missing value for parameter "'//tag(1:ltag)//'". Skipping line',0)
              qerror=.true.; exit
             elseif (qvector) then
              call warning(whoami, 'Missing closing bracket for vector parameter "'//tag(1:ltag)//'". Skipping line',0)
              qerror=.true.; exit
             else ! store tag and val
              call adjustleft(val(1:lval)); lval=len_trim(val(1:lval));
              call message(whoami, tag(1:ltag)//' <= '//val(1:lval))
              j=params_uadd(parameters,tag,val,ltag,lval)
              ltag=0; lval=0
              cmdline=cmdline(i:l) ! remove tag/val pair from string
              call adjustleft(cmdline) ! remove tag/val pair from string
              l=len_trim(cmdline)
              i=0
             endif ! lval
             qtag=.true. ; ltag=0
             qval=.false.
            endif ! qtag
!
           elseif (any(vector_open.eq.c)) then ! encountered an opening vector bracket
            if (qtag) then ! opening bracket not allowed
             call warning(whoami, 'Unexpected opening bracket when processing parameter "'//tag(1:ltag)//'". Skipping line',0)
             qerror=.true.; exit
            elseif (qeq) then
             call warning(whoami, 'Unexpected opening bracket when processing parameter "'//tag(1:ltag)//'". Skipping line',0)
             qerror=.true.; exit
            elseif (qval) then
             if (qvector) then ! already read one bracket
             call warning(whoami, &
                           'Unexpected opening bracket when processing vector parameter "'//tag(1:ltag)//'". Skipping line',0)
             qerror=.true.; exit
             else
              qvector=.true. ! the only legitimate occurrence
             endif
            endif
!
           elseif (any(vector_close.eq.c)) then ! encountered a closing vector bracket
! write(0,*) '**', qtag, qeq, qval, qvector, tag, ltag, val, lval; !aa
            if (qtag) then ! closing bracket not allowed
             call warning(whoami, 'Unexpected closing bracket when processing parameter 1"'//tag(1:ltag)//'". Skipping line',0)
! write(0,*) cmdline(i:l) ! aa
             qerror=.true.; exit
            elseif (qeq) then
             call warning(whoami, 'Unexpected closing bracket when processing parameter "'//tag(1:ltag)//'". Skipping line',0)
             qerror=.true.; exit
            elseif (qval) then
             qvector=.false. ! the only legitimate occurrence
! peek ahead to make sure we have a spacer, if not, issue a warning and skip line
             c=cmdline(i+1:i+1)
             if ( (i.lt.l) .and. (all(space.ne.c)) .and. (all(comment.ne.c))) then
              call warning(whoami, 'Error reading vector value for parameter "'//tag(1:ltag)//'". Skipping line',0)
              qerror=.true.; exit
! otherwise, attempt to add tag/value pair
             else
              if (lval.eq.0) then
               call warning(whoami, 'Missing vector value for parameter "'//tag(1:ltag)//'". Skipping line',0)
               qerror=.true.; exit
              endif
              call adjustleft(val(1:lval)); lval=len_trim(val(1:lval));
              if (lval.eq.0) then
               call warning(whoami, 'Missing vector value for parameter "'//tag(1:ltag)//'". Skipping line',0)
               qerror=.true.; exit
              else ! add pair
               call adjustleft(val(1:lval)); lval=len_trim(val(1:lval));
               call message(whoami, tag(1:ltag)//'[VECTOR] <= ('//val(1:lval)//')')
               j=params_uadd(parameters,tag,val,ltag,lval)
               ltag=0; lval=0
               cmdline=cmdline(i+1:l) ! remove tag/val pair from string
               call adjustleft(cmdline) ! remove tag/val pair from string
               l=len_trim(cmdline)
               i=0
              endif ! lval
              qtag=.true. ; ltag=0
              qval=.false.
             endif ! (i.lt.l)
            else ! no imput mode (qtag/qeq/qval) : this should never happen
             call warning(whoami, 'Misplaced closing bracket in input. Skipping line.',0)
             qerror=.true.; exit
            endif
!
           elseif (any(allowed2(1:allowed_len).eq.c)) then ! check that the characters are allowed
!
            if (qtag) then
             if (all(lower2(1:26).ne.c).and.all(upper2(1:26).ne.c)) then
              if (ltag.eq.0) then
               call warning(whoami, 'Parameter names must start with a letter. Skipping line.',0)
               qerror=.true.; exit
              elseif (all(digits2(1:10).ne.c).and.underscore.ne.c) then
               call warning(whoami, 'Illegal character in parameter name. Skipping line.',0)
               qerror=.true.; exit
              endif
             endif
             ltag=ltag+1; tag(ltag:ltag)=c
            elseif (qeq) then
             call warning(whoami, 'Missing value for parameter "'//tag(1:ltag)//'". Skipping line',0)
             qerror=.true.; exit
            elseif (qval) then
             lval=lval+1; val(lval:lval)=c
            else
             call warning(whoami, 'Unknown error. Skipping line',0)
             qerror=.true.; exit
            endif
!
           else
             call warning(whoami, 'UNRECOGNIZED CHARACTER "'//c//'". SKIPPING LINE.',0)
             exit
           endif ! character loop
         enddo ! while l.gt.1
        else ! end of file
         exit
        endif
       enddo ! over all lines in the file
!



       if (qerror) then
         call warning(whoami, 'ERROR(S) FOUND IN INPUT.',0)
       else
         call message(whoami, 'Input file read.')
       endif
!
       end subroutine parse_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine list_params( &



       & )
!
       use output
       character(len=11), parameter :: whoami = 'LIST_PARAMS'
       integer :: i



!
       if (.not.parser_initialized) call parser_init()
!



        call message(whoami,'THE FOLLOWING PARAMETERS ARE DEFINED')
        call message(whoami,'====================================')
        do i=1,parameters%last
         call message(whoami,tab//parameters%tag(i)(1:parameters%tlen(i))//' = "'//parameters%val(i)(1:parameters%vlen(i))//'"')
        enddo
        call message(whoami,'====================================')



       end subroutine list_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! auxiliary functions (they need not be part of this module)
       function atoi(a, invalid)
! NOTE: no overflow check yet
       use output, only: warning
       integer, optional, intent(in) :: invalid
       integer :: atoi, i, l, j, k, sgn, base, missing
       character(len=*), intent(in) :: a
       character(len=len(a)) :: b
       character(len=4), parameter :: whoami = 'ATOI'
       integer :: flag(10)
! convert string to integer
       if (present(invalid)) then ; missing=invalid ; else ; missing=-999 ; endif
       i=0
!
       b=a
       call adjustleft(b)
       l=len_trim(b)
       if (l.ge.1) then
        if (b(1:1).eq.hyphen) then
         sgn=-1
         b(1:l-1)=b(2:l); l=l-1;
         if (l.eq.0) i = -missing ; ! only a hyphen present : will multiply by -1 and quit
        else
         sgn=1
        endif
       else
        i=missing; sgn=1
       endif ! l.ge.1
!
       base=1
       do j=l,1,-1
        where(digits2.eq.b(j:j)); flag=1 ; elsewhere; flag=0 ; endwhere; k=sum(maxloc(flag))-1
        if (all(flag.eq.0)) then
         call warning(whoami, 'ERROR CONVERTING STRING "'//b(1:l)//'" TO INTEGER.',-1)
         i=missing; sgn=1;
         exit
        else
         i=i+base*k
         base=base*10
        endif
       enddo
       atoi=i*sgn
       end function atoi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atoiv(a,n)
       use output, only: warning
       character(len=*) :: a
       character(len=200), allocatable :: b(:)
       character(len=5), parameter :: whoami = 'ATOIV'
       integer :: n, i, j
       integer, dimension(n) :: atoiv
!
       i=numword(a)
       if (n.le.0) then
        call warning(whoami,' Vector has nonpositive dimension',-1)
        return
       elseif (i.ne.n) then
        call warning(whoami,' Vector dimension mismatch',-1)
       endif
       i=min(i,n)
       atoiv=0
       allocate(b(i))
       read(a,*) b
       do j=1,i
        atoiv(j)=atoi(b(j))
       enddo
       deallocate(b)
!
       end function atoiv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atof(a,invalid)
       use output, only: warning
       real*8, optional, intent(in) :: invalid
       real*8 :: atof, missing
       real*8 :: f
       integer :: i, l, j, k, sgn, base
       character(len=*) :: a
       character(len=len(a)) :: b
       character(len=4), parameter :: whoami = 'ATOF'
       logical :: fraction
       integer :: flag(10)
! convert string to floating point number
       if (present(invalid)) then ; missing=invalid ; else ; missing=-99999. ; endif
       f=0;
       fraction=.false.
!
       b=a
       call adjustleft(b)
       l=len_trim(b)
       if (l.ge.1) then
        if (b(1:1).eq.hyphen) then
         sgn=-1
         b(1:l-1)=b(2:l); l=l-1;
         if (l.eq.0) f = -missing ; ! only a hyphen present : will multiply by -1 and quit
        else
         sgn=1
        endif
       else
        f=missing; sgn=1
       endif ! l.ge.1
!
       base=0
       do j=l,1,-1
        if (b(j:j).eq.decimal) then
         if (fraction) then ! two decimal points are invalid
          call warning(whoami, 'ERROR CONVERTING STRING "'//b(1:l)//'" TO REAL.',-1)
          f=missing; sgn=1;
          exit
         else
          fraction=.true.
          do while (base.gt.0)
           f=f/10.
           base=base-1
          enddo
         endif
        else
         where(digits2.eq.b(j:j)); flag=1 ; elsewhere; flag=0 ; endwhere; k=sum(maxloc(flag))-1
         if (all(flag.eq.0)) then
          call warning(whoami, 'ERROR CONVERTING STRING "'//b(1:l)//'" TO REAL.',-1)
          f=missing; sgn=1;
          exit
         else
          f=f+1.0d0*(10.0d0**base)*k
          base=base+1
         endif
        endif
       enddo
       atof=f*sgn
       end function atof
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atofv(a,n)
       use output, only: warning
       character(len=*) :: a
       character(len=200), allocatable :: b(:)
       character(len=5), parameter :: whoami = 'ATOFV'
       integer :: n, i, j
       real*8, dimension(n) :: atofv
!
       i=numword(a)
       if (n.le.0) then
        call warning(whoami,' Vector has nonpositive dimension',-1)
        return
       elseif (i.ne.n) then
        call warning(whoami,' Vector dimension mismatch',-1)
       endif
       i=min(i,n)
       atofv=0
       allocate(b(i))
       read(a,*) b
       do j=1,i
        atofv(j)=atof(b(j))
       enddo
       deallocate(b)
!
       end function atofv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atol(a)
       use output, only: error
       logical :: atol
       character(len=*) :: a
       character(len=200) :: b
       character(len=4), parameter :: whoami = 'ATOL'
       integer :: l
!
       b=a
       call adjustleft(b)
       l=len_trim(b)
       select case(a)
!
        case('true', '.true.', '.TRUE.', 'TRUE', 'YES', 'ON', 'yes', 'on', 'y', 'Y');
         atol=.true.
        case('false', '.false.', '.FALSE.', 'FALSE', 'NO', 'OFF', 'no', 'off', 'n', 'N');
         atol=.false.
        case default
         call error(whoami, 'ERROR CONVERTING STRING "'//b(1:l)//'" TO BOOLEAN.',-1)
        atol=.false.
       end select
!
       end function atol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function itoa(i)
   integer :: i
   character(len=((ceiling(log10(1.0*abs(i)+1.)))+1)) :: itoa
   character(len=80) :: b
   write(b,*) i
   b=adjustl(b)
   itoa=b(1:len_trim(b))
   end function itoa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ftoa(f)
   real*8 :: f
   character(len=20) :: ftoa
   write(ftoa,'(G20.10)') f
   end function ftoa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ltoa(l)
   logical :: l
   character(len=3) :: ltoa
   if (l) then ; ltoa='YES' ; else ; ltoa='NO ' ; endif
   end function ltoa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine toupper(a)
       character(len=*) :: a
       integer :: i,j
!
       do j=1, len_trim(a)
        do i=1,26; if (lower2(i).eq.a(j:j)) then; a(j:j)=upper2(i); exit; endif; enddo
       enddo
!
       end subroutine toupper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine tolower(a)
       character(len=*) :: a
       integer :: i,j
!
       do j=1, len_trim(a)
        do i=1,26; if (upper2(i).eq.a(j:j)) then; a(j:j)=lower2(i); exit; endif; enddo
       enddo
!
       end subroutine tolower
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine uses a custom definition of white space
       subroutine adjustleft(a, space)
        character(len=*) :: a
        character, optional :: space(:)
        integer :: l, i, j
        l=len(a)
        if (l.gt.0) then
         i=1
         do while (i.le.l)
          if (all(space.ne.a(i:i))) exit
          i=i+1
         enddo
! move string left
         j=1
         do while (j.le.l-i+1)
          a(j:j)=a(i+j-1:i+j-1)
          j=j+1
         enddo
! pad with blanks
         do while (j.le.l)
          a(j:j)=' '
          j=j+1
         enddo
        endif
! write(0,*) 'AL***:',a
!
       end subroutine adjustleft
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function numword(a) ! returns the number of blank-separated words in a string
   character(len=*) :: a
   character, parameter :: tab=char(9)
   character(len=2), parameter :: space = ' '//tab
   character, parameter :: space2(2) = (/' ',tab/)
   integer :: n, numword,i,j,l

   l=len_trim(a);
   n=0
   if (l.eq.0.or.any(a(1:l).eq.space2)) return ! this is string comparison: any applies to space2
   j=1
   i=1
   do while (j.gt.0)
    j=scan(a(i:l),space)
! write(0,*) i,j,n
! pause
    if (j.gt.1) then
     n=n+1
    endif
    i=i+j
   enddo
   if (i.le.l) then
    if (all(a(i:l).ne.space2)) n=n+1
   endif
   numword=n
  end function numword
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pop_string(a,n) result(b) ! returns the next white-space-delimited word in a string, and removes it from the string
   character(len=*) :: a
   character(len=len(a)) :: b
   character, parameter :: tab=char(9)
   character(len=2), parameter :: space = ' '//tab
   character, parameter :: space2(2) = (/' ',tab/)
   integer, optional, intent(inout) :: n
   integer :: j,l

   if (present(n)) a(max(0,n):)='' ! erase string beyond length n
!
   call adjustleft(a,space2)
   l=len_trim(a);
   j=scan(a(1:l),space)
   if (j.eq.0) j=l+1 ! there must be only one word
   b=a(1:j-1)
   a=a(j:l)
   if (present(n)) n=l-j+1
  end function pop_string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function find_tag(a,tag,n)
  character(len=*), intent(inout) :: a ! because we might be removing white space
  character(len=*) :: tag
  integer, optional, intent(inout) :: n
  integer :: ltag, i, j, k, find_tag
  character(len=len(tag)) :: copy
!
  character(len=2), parameter :: space = ' '//tab
  character, parameter :: space2(2) = (/' ',tab/)
!
  copy=tag; call adjustleft(copy,space2); ltag=len_trim(copy)
  if (present(n)) a(max(0,n):)='' ! erase string beyond length n
  call adjustleft(a,space2); n=len_trim(a);
  if (ltag.gt.0.and.n.ge.ltag) then ! proceed only if tag nontrivial and could fit inside string
! note: we are only looking for matches that correspond to a new word (not mid-word)
   if (a(1:ltag).eq.copy(1:ltag)) then ! string begins with tag
    j=1
   else ! tag might be somewhere in the middle
    j=-1
    do i=1,size(space2)
     k=index(a(1:n), space(i:i)//copy(1:ltag) ); ! pre-pend an instance of white space to tag
     if (k.gt.0) j=min(j,k) ! take the first occurrence
    enddo
    j=j+1
   endif
  else
   j=0
  endif
!
  find_tag=j
  end function find_tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function remove_tag(a,tag,n) result(j)
  character(len=*) :: a
  character(len=*) :: tag
  integer, optional, intent(inout) :: n
  integer :: ltag, i, j, k, l
  character(len=len(a)) :: right
  character(len=len(tag)) :: copy
!
  character(len=2), parameter :: space = ' '//tab
  character, parameter :: space2(2) = (/' ',tab/)
!
  copy=tag; call adjustleft(copy,space2); ltag=len_trim(copy)
  if (present(n)) a(max(0,n):)='' ! erase string beyond length n
  call adjustleft(a,space2); l=len_trim(a);
  if (ltag.gt.0.and.l.ge.ltag) then ! proceed only if tag nontrivial and could fit inside string
! note: we are only looking for matches that correspond to a new word (not mid-word)
   if (a(1:ltag).eq.copy(1:ltag)) then ! string begins with tag
    j=1
   else ! tag might be somewhere in the middle
    j=-1
    do i=1,size(space2)
     k=index(a(1:l), space(i:i)//copy(1:ltag) ); ! pre-pend an instance of white space to tag
     if (k.gt.0) j=min(j,k) ! take the first occurrence
    enddo
    j=j+1
   endif
! remove tag and parameter from string, if tag found
   if (j.gt.0) then
    right=a(j+ltag:)
    a(j:)=right
    if (present(n)) n=len_trim(a)
   endif
  endif
!
  end function remove_tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_remove_parameter(a,tag,n) result(b)
  character(len=*) :: a
  character(len=*) :: tag
  integer, optional, intent(inout) :: n
  integer :: ltag, i, j, k, l
  character(len=len(a)) :: b, right
  character(len=len(tag)) :: copy
!
  character(len=2), parameter :: space = ' '//tab
  character, parameter :: space2(2) = (/' ',tab/)
!
  b=''
  copy=tag; call adjustleft(copy,space2); ltag=len_trim(copy)
  if (present(n)) a(max(0,n):)='' ! erase string beyond length n
  call adjustleft(a,space2); l=len_trim(a);
  if (ltag.gt.0.and.l.ge.ltag) then ! proceed only if tag nontrivial and could fit inside string
! note: we are only looking for matches that correspond to a new word (not mid-word)
   if (a(1:ltag).eq.copy(1:ltag)) then ! string begins with tag
    j=1
   else ! tag might be somewhere in the middle
    j=-1
    do i=1,size(space2)
     k=index(a(1:l), space(i:i)//copy(1:ltag) ); ! pre-pend an instance of white space to tag
     if (k.gt.0) j=min(j,k) ! take the first occurrence
    enddo
    j=j+1
   endif
! remove tag and parameter from string, if tag found
   if (j.gt.0) then
    right=a(j+ltag:)
    b=pop_string(right)
    a(j:)=right
    if (present(n)) n=len_trim(a)
   endif
  endif
!
  end function get_remove_parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module parser
