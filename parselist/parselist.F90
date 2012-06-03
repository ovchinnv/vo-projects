/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
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
      module parselist
!
      contains
       subroutine ilist_parse(list, str)
! parses a _unique_ list of ints stored in a character string
        use ivector
        use parser
!
        implicit none
!
        integer :: imin, imax
        character(len=*) :: str
        type (int_vector) :: list
! locals
        integer :: i, j, k, strl, inode, jnode, kstep
        integer, parameter :: missing=-99999999
        character(len=11) :: whoami
!
        data whoami /' PARSE_LIST>'/
!
        call int_vector_init(list)
!
        inode=missing; jnode=missing;
        strl=len(str)
!
        do
         strl=min(max(0,strl),len(str));str(strl+1:)='';call adjustleft(str,(/' ',tab/));strl=len_trim(str)
         i=find_tag(str, 'THRU', strl)
         if (i.gt.1.or.i.eq.0) then
! first, add previous node number (no error if invalid)
            if (inode.ne.missing) j=int_vector_uadd(list,inode)
         elseif (i.eq.1) then ! have to deal with THRU
            jnode=atoi(get_remove_parameter(str, 'THRU', strl), missing)
! check for 'STEP' keyword
            strl=min(max(0,strl),len(str));str(strl+1:)='';call adjustleft(str,(/' ',tab/));strl=len_trim(str)
            i=find_tag(str, 'STEP', strl)
            if (i.eq.1) then ! STEP keyword present
             kstep=atoi(get_remove_parameter(str, 'STEP', strl), missing)
!! allow (almost) any value of STEP for flexibility
             if (kstep.eq.missing) then
              call warning(whoami, ' INVALID STEP VALUE SPECIFIED. USING STEP=1.', 0)
              kstep=1
             endif
            else ! STEP not specified -- using 1
             kstep=1
            endif
!
            if (inode.ne.missing.and.jnode.ne.missing) then
             do k=inode, jnode, kstep
               j=int_vector_uadd(list,k)
             enddo
            else
              call warning(whoami, ' INVALID RANGE SPECIFIED. SKIPPING ENTRY.', 0)
            endif
         endif ! i
! read the next number
         if (strl.lt.1) exit
!
         inode=atoi(pop_string(str,strl)) ; strl=len_trim(str)
        enddo
!
       end subroutine ilist_parse
      end module parselist
