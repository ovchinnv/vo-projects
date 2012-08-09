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
      module freeio
!
       character, parameter, private :: comment(4) = (/'*', '#', '%', '!'/)
!
       contains
!
       subroutine atomid_coor_read(fid,r) ! read coordinates in list-directed format
       use output, only: error, warning, message
       use parser, only: adjustleft, numword
       use psf
       implicit none
       integer :: fid ! input file handle
       real*8 :: r(:,:)
!
       integer :: i, j
       character(len=16), parameter :: whoami = 'ATOMID_COOR_READ'
!
       integer :: ioerr
       integer :: l=0, atomid, natom
       real*8 :: x,y,z
       character(len=400) :: cmdline
       logical :: flags(atoms%last)
       logical :: found=.false.
!
       natom=size(r,2)
       flags=.false.
!
       i=0
!
       do while (.true.)
        read(fid,'(A)',IOSTAT=ioerr) cmdline ! if running in parallel, then only the root node is passed a valid handle
        if (ioerr.eq.0) then
         call adjustleft(cmdline)
         l=len_trim(cmdline)
! ignore lines that begin with comment symbols
         if (any(comment.eq.cmdline(1:1))) cycle
! skip empty lines
         if (l.eq.0) cycle
!
         if (numword(cmdline(1:l)).lt.4) then
          call warning(whoami, cmdline(1:l),0)
          call warning(whoami, 'LINE TOO SHORT. SKIPPING.',0)
          cycle
         endif
!
         read(cmdline,*) atomid, x, y, z
!! match atom coordinate entry with structure
         if (atomid.lt.1) then
          call warning(whoami, cmdline(1:len_trim(cmdline)),0)
          call warning(whoami, 'NEGATIVE ATOMID READ. ABORT. SOME COORDINATES UNDEFINED',0)
          elseif (natom.lt.atomid) then
          call error(whoami, 'COORDINATE ARRAY HAS INCORRECT DIMENSIONS. ABORT.',-1)
          return
         endif
! find index
         found=.false.
         do j=atomid, atoms%last
          if (atoms%atomid(j).eq.atomid) then
           found=.true.
           exit
          endif
         enddo
         if (.not.found) then
          do j=atomid-1, 1
           if (atoms%atomid(j).eq.atomid) then
            found=.true.
            exit
           endif
          enddo
         endif
!
         if (.not.found) then
          call warning(whoami, cmdline(1:len_trim(cmdline)),0)
          call warning(whoami, 'CANNOT FIND ATOMID IN STRUCTURE FILE. SKIPPING LINE.',0)
         else
          r(:,j)=(/x,y,z/);
          flags(atomid)=.true.
         endif
         i=i+1 ! increment atom count
!
        else ! ioerr -- EOF
         exit
        endif
!
       enddo ! while(.true.)
!
       if (i.ne.atoms%last) call warning(whoami, 'NUMBER OF VALID ENTRIES IN COORDINATE FILE INCONSISTENT WITH STRUCTURE.',0)
       if (.not.all(flags)) call warning(whoami, 'SOME COORDINATES WERE MISSING.',0)
!
       call message(whoami, 'Coordinate file read.')
!
       end subroutine atomid_coor_read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       subroutine atomid_coor_write(fid,r)
       use psf
       use output, only: warning, message
       implicit none
       integer :: fid
       real*8 :: r(:,:)
       integer :: natom, i
       character(len=13), parameter :: fmt='(I10,3E25.15)'
       character(len=17), parameter :: whoami = 'ATOMID_COOR_WRITE'
!
       natom=size(r,2)
!
       if (natom.ne.atoms%last) then
        call warning(whoami, 'COORDINATE ARRAY HAS INCORRECT DIMENSIONS. ABORT.',-1)
       else
        write(fid,'(A)') '* FREE FORMAT COORDINATE FILE WRITTEN BY DYNAMO PROGRAM'
        do i=1, natom
         write(fid, fmt) atoms%atomid(i), r(1:3,i), atoms%segid(i)
        enddo
        call message(whoami, 'Coordinate file writen.')
       endif
!
       end subroutine atomid_coor_write
!
      end module freeio
