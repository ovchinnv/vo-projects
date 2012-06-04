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
module pdbio
 contains
 subroutine pdb_read(fid,r,O,B)
 use psf
 use parser, only: atoi, numword, adjustleft, toupper
 use output, only: error, warning, message
 implicit none
 character, parameter :: comment(1)=(/'*'/)
 integer :: fid
 real*8 :: r(:,:) ! coordinates
 real*8, optional :: B(:) ! B-factor column
 real*8, optional :: O(:) ! occupancy column
 character(len=8), parameter :: whoami = 'PDB_READ'
 character(len=200) :: cmdline
 character(len=8) :: keyword
 integer :: n=-1
 integer :: natom
! character(len=100) :: fmt= 'A6,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3,2F6.2,1X,3X,2X,A4'
 character(len=100), parameter :: fmt='(A6,I5,1X,A4,1X,A3,2X,A4,4X,3F8.3,2F6.2,1X,3X,2X,A4)' ! resid is a string
 integer :: atomid, i, ioerr, j
 real*8 :: x,y,z
 character(len=8) :: aname, resname, segid, resid
 logical :: flags(atoms%last)
 logical :: found=.false.
 logical :: qo, qb ! B/occupancy flags
 real*8 :: occ, bf
!
!
 natom=size(r,2)
 flags=.false.
!
 qb=present(B)
 qo=present(O)
!
 if (qb) then
  if (size(B).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF B-FACTOR ARRAY',-1)
  endif
 endif
!
 if (qo) then
  if (size(O).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF OCCUPANCY ARRAY',-1)
  endif
 endif
! remove comments at the beginning, if any
 do while (.true.)
  read(fid,'(A)',IOSTAT=ioerr) cmdline
  if (ioerr.eq.0) then
   if (any(comment.eq.cmdline(1:1))) cycle ! note: this will not get rid of 'REMARK' fields in PDBs
   exit
  else
   call error(whoami, 'Unexpected end of file',-1)
   return
  endif
 enddo
!
 i=0
 do while (ioerr.eq.0)
! process command line
! process command line
! only lines that begin with 'ATOM'/'HETATM' are processed; 'TER ' or 'END ' indicates end of read
  keyword(1:6)=cmdline(1:6);
  call toupper(keyword)
!
  select case(keyword(1:6))
   case('ATOM  ', 'HETATM');
   read(cmdline,fmt) keyword, atomid, aname, resname, resid, x, y, z, occ, bf, segid
   call adjustleft(cmdline)
   call adjustleft(segid)
   call adjustleft(resid)
   call adjustleft(resname)
   call adjustleft(aname)
! match atom coordinate entry with structure
   if (atomid.lt.1) then
    call warning(whoami, cmdline(1:len_trim(cmdline)),0)
    call warning(whoami, 'Nonpositive atom ID read. Skipping entry. Some coordinates may be undefined',0)
! try to read next line
    read(fid,'(A)',IOSTAT=ioerr) cmdline
    cycle
   elseif (natom.lt.atomid) then
    call error(whoami, 'Coordinate array has incorrect dimensions. Abort.',-1)
    return
   endif
! find index of the atom in structure
   found=.false.
   do j=atomid, atoms%last ! first, a forward search, assuming ordering; usually atomids will be [1..natom]
    if (atoms%atomid(j).eq.atomid) then
     found=.true.
     exit
    endif
   enddo
   if (.not.found) then ! try a reverse search, in case the file is disordered
    do j=min(atomid-1,atoms%last), 1, -1 ! make sure indices stay within bounds
     if (atoms%atomid(j).eq.atomid) then
      found=.true.
      exit
     endif
    enddo
   endif
!
   if (.not.found.or.((atoms%segid(j).ne.segid).or.(atoms%resid(j).ne.resid).or.&
          (atoms%resname(j).ne.resname).or.(atoms%aname(j).ne.aname))) then
    call warning(whoami, cmdline(1:len_trim(cmdline)),0)
    call warning(whoami, 'Cannot find corresponding entry in structure file. Skipping line.',0)
! try to read next line
    read(fid,'(A)',IOSTAT=ioerr) cmdline
    cycle
   else
    r(:,j)=(/x,y,z/);
    if (qb) B(j)=bf; ! b-factor
    if (qo) O(j)=occ; ! occupancy
    flags(atomid)=.true.
   endif
!
   i=i+1 ! increment atom count
! try to read next line
   read(fid,'(A)',IOSTAT=ioerr) cmdline
!**********************************************************************************
   case('TER   ','END   ');
   exit ! loop over lines
   case default ! keyword not recognized; assume that we can continue (we are not enforcing the PDB standard)
! try to read next line
   read(fid,'(A)',IOSTAT=ioerr) cmdline
  end select
!
 enddo ! while()
!
! the number of lines processed is i; were all of them valid atom entries?
 if (i.ne.atoms%last) call warning(whoami, 'NUMBER OF ATOMS IN COORDINATE FILE INCONSISTENT WITH STRUCTURE.',0)
 if (.not.all(flags)) call warning(whoami, 'SOME COORDINATES WERE MISSING.',0)
!
 call message(whoami, 'PDB file read.')
!
 end subroutine pdb_read
!
!**********************************************************************************
!
 subroutine pqr_read(fid,r,O,B)
! format very similar to PDB; whitespace-delimited free-form instead of fixed form
! occupancy contains charge; b-factor column contains radius
! no segid field; APBS PQR has a chain ID, which is currently unsupported
 use psf
 use parser, only: atoi, numword, adjustleft, toupper
 use output, only: error, warning, message
 implicit none
 character, parameter :: comment(1)=(/'*'/)
 integer :: fid
 real*8 :: r(:,:) ! coordinates
 real*8, optional :: B(:) ! B-factor column
 real*8, optional :: O(:) ! occupancy column
 character(len=8), parameter :: whoami = 'PQR_READ'
 character(len=200) :: cmdline
 character(len=8) :: keyword
 integer :: n=-1
 integer :: natom
 character(len=100) :: fmt='*'
 integer :: atomid, i, ioerr, j
 real*8 :: x,y,z
 character(len=8) :: aname, resname, resid
 logical :: flags(atoms%last)
 logical :: found=.false.
 logical :: qo, qb ! B/occupancy flags
 real*8 :: occ, bf
!
!
 natom=size(r,2)
 flags=.false.
!
 qb=present(B)
 qo=present(O)
!
 if (qb) then
  if (size(B).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF B-FACTOR ARRAY',-1)
  endif
 endif
!
 if (qo) then
  if (size(O).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF OCCUPANCY ARRAY',-1)
  endif
 endif
! remove comments at the beginning, if any
 do while (.true.)
  read(fid,'(A)',IOSTAT=ioerr) cmdline
  if (ioerr.eq.0) then
   if (any(comment.eq.cmdline(1:1))) cycle
   exit
  else
   call error(whoami, 'UNEXPECTED END OF FILE',-1)
   return
  endif
 enddo
!
 i=0
 do while (.true.)
! process command line
! only lines that begin with 'ATOM'/'HETATM' are processed; 'TER ' or 'END ' indicates end of read
  keyword(1:6)=cmdline(1:6);
  call toupper(keyword)
!
  select case(keyword(1:6))
   case('ATOM  ', 'HETATM');
   read(cmdline,fmt) keyword, atomid, aname, resname, resid, x, y, z, occ, bf
   call adjustleft(cmdline)
   call adjustleft(resid)
   call adjustleft(resname)
   call adjustleft(aname)
! match atom coordinate entry with structure
   if (atomid.lt.1) then
    call warning(whoami, cmdline(1:len_trim(cmdline)),0)
    call warning(whoami, 'NEGATIVE ATOM ID READ. ABORT. SOME COORDINATES MAY BE UNDEFINED',0)
   elseif (natom.lt.atomid) then
    call error(whoami, 'COORDINATE ARRAY HAS INCORRECT DIMENSIONS. ABORT.',-1)
    return
   endif
! find index of the atom in structure
   found=.false.
   do j=atomid, atoms%last ! first, a forward search, assuming ordering
    if (atoms%atomid(j).eq.atomid) then
     found=.true.
     exit
    endif
   enddo
   if (.not.found) then ! try a reverse search, in case the file is disordered
    do j=atomid-1, 1
     if (atoms%atomid(j).eq.atomid) then
      found=.true.
      exit
     endif
    enddo
   endif
!
   if (.not.found.or.((atoms%resid(j).ne.resid).or.&
          (atoms%resname(j).ne.resname).or.(atoms%aname(j).ne.aname))) then
    call warning(whoami, cmdline(1:len_trim(cmdline)),0)
    call warning(whoami, 'CANNOT FIND CORRESPONDING ENTRY IN STRUCTURE FILE. SKIPPING LINE.',0)
   else
    r(:,j)=(/x,y,z/);
    if (qb) B(j)=bf; ! b-factor
    if (qo) O(j)=occ; ! occupancy
    flags(atomid)=.true.
   endif
!
   i=i+1 ! increment atom count
! try to read next line
   read(fid,'(A)',IOSTAT=ioerr) cmdline
   if (ioerr.ne.0) exit
!**********************************************************************************
   case('TER   ','END   ');
   exit ! loop over lines
  end select
!
 enddo ! while(.true.)
!
! the number of lines processed is i; were all of them valid atom entries?
 if (i.ne.atoms%last) call warning(whoami, 'NUMBER OF ATOMS IN COORDINATE FILE INCONSISTENT WITH STRUCTURE.',0)
 if (.not.all(flags)) call warning(whoami, 'SOME COORDINATES WERE MISSING.',0)
!
 call message(whoami, 'PQR file read.')
!
 end subroutine pqr_read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine pdb_write(fid, r, O, B)
 use psf
 use output, only: error, message
 implicit none
 integer :: fid
 real*8 :: r(:,:)
 real*8, optional :: B(:) ! B-factor column
 real*8, optional :: O(:) ! occupancy column
 integer :: natom, i, resnum
 character(len=100) :: fmt='A6,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3,2F6.2,1X,3X,2X,A4'
 character(len=9), parameter :: whoami = 'PDB_WRITE'
 logical :: qo, qb ! B/occupancy flags
 real*8 :: occ, bf
!
 natom=size(r,2)
!
 qb=present(B)
 qo=present(O)
!
 if (qb) then
  if (size(B).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF B-FACTOR ARRAY',-1)
  endif
 endif
!
 if (qo) then
  if (size(O).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF OCCUPANCY ARRAY',-1)
  endif
 endif
!
 if (natom.ne.atoms%last) then
   call error(whoami, 'COORDINATE ARRAY HAS INCORRECT DIMENSIONS. ABORT.',-1)
 else
  write(fid,'(A)') 'REMARK: PDB FILE WRITTEN BY DYNAMO PROGRAM'
!
  bf=0d0
  occ=0d0
  do i=1, natom
   read(atoms%resid(i),*) resnum ! set residue numbers to resid (ad-hoc, but we do not store any info on them)
   if (qb) bf= B(i); ! b-factor
   if (qo) occ=O(i); ! occupancy
   write(fid, fmt) 'ATOM  ',atoms%atomid(i),atoms%aname(i),atoms%resname(i),atoms%resid(i),r(1:3,i),occ,bf,atoms%segid(i)
  enddo
!
 endif
!
 call message(whoami, 'PDB file writen.')
!
 end subroutine pdb_write
 subroutine pqr_write(fid, r, O, B)
! format very similar to PDB; whitespace-delimited free-form instead of fixed form
! occupancy contains charge; b-factor column contains radius
! no segid field; APBS PQR has a chain ID, which is currently unsupported
 use psf
 use output, only: error, message
 implicit none
 integer :: fid
 real*8 :: r(:,:)
 real*8, optional :: B(:) ! B-factor column
 real*8, optional :: O(:) ! occupancy column
 integer :: natom, i, resnum
 character(len=100) :: fmt='A6,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3,2F6.2,1X,3X,2X,A4'
 character(len=9), parameter :: whoami = 'PQR_WRITE'
 logical :: qo, qb ! B/occupancy flags
 real*8 :: occ, bf
!
 natom=size(r,2)
!
 qb=present(B)
 qo=present(O)
!
 if (qb) then
  if (size(B).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF B-FACTOR ARRAY',-1)
  endif
 endif
!
 if (qo) then
  if (size(O).ne.natom) then
    call error(whoami, 'INCONSISTENT SIZE OF OCCUPANCY ARRAY',-1)
  endif
 endif
!
 if (natom.ne.atoms%last) then
   call error(whoami, 'COORDINATE ARRAY HAS INCORRECT DIMENSIONS. ABORT.',-1)
 else
  write(fid,'(A)') 'REMARK: PQR FILE WRITTEN BY DYNAMO PROGRAM'
!
  bf=0d0
  occ=0d0
  do i=1, natom
   read(atoms%resid(i),*) resnum ! set residue numbers to resid (ad-hoc, but we do not store any info on them)
   if (qb) bf= B(i); ! b-factor
   if (qo) occ=O(i); ! occupancy
   write(fid, fmt) 'ATOM  ',atoms%atomid(i),atoms%aname(i),atoms%resname(i),atoms%resid(i),r(1:3,i),occ,bf
  enddo
!
 endif
!
 call message(whoami, 'PQR file writen.')
!
 end subroutine pqr_write
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module pdbio
