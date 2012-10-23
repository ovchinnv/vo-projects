/*#define __WRN(__WHO,__MSG) write(0,*) 'WARNING FROM: ',__WHO,': ',__MSG*/
/*#define __PRINT(__MSG) write(0,'(A)') __MSG*/
/*#define __PRINT(__MSG) call plainmessage(__MSG)*/
/*#define __PRINTL(__MSG,__LEVEL) call plainmessage(__MSG,__LEVEL)*/
/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
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
module confcons
 !**CHARMM_ONLY**! __DEP_KINDS
 implicit none
 private
!
 !
 ! in this module :
 ! the CONFormational CONSistency module
 ! is applicable to groups of atoms that have
 ! rotational symmetry because some od the constituent
 ! atoms are indistinguishable, e.g. H1, H2, and H3 in
 ! a methyl group. In specifying the coordinates of a molecule
 ! with such a group, the atoms can be arranged in three completely
 ! equivalent ways, relative to a reference conformation.
 ! These are H1/H2/H3, H2/H3/H1, H3/H1/H2. They are related by
 ! rotations around the bond that involves the parent atom
 ! (e.g. methyl carbon). Note that the conformation H1/H3/H2
 ! is not allowed because it corresponds to a change of (pseudo)
 ! chirality, which is assumed to be fixed (see chirality module)
 ! The purpose of this module is to identify all such symmetry groups
 ! in a molecular coordinate set and choose the rotation that
 ! brings the coordinates into the closest proximity (in the sense
 ! of a certain RMSD ) to the coordinates in a reference set.
 ! This is done according to a set of rules, which can be
 ! specified by the user; a default set for the charmm22 topology
 ! is defined in this module
!
!
!
  type confcons_rule
   character(len=8) :: resname
   character(len=8), pointer :: test_atoms(:)
   character(len=8), pointer :: orient_atoms(:)
   character(len=8), pointer :: orient_atom_permutations(:,:)
   character(len=8), pointer :: other_atom_permutations(:,:)
! contains ! some compilers lack this feature
! final :: confcons_rule_finalize
  end type confcons_rule
!
  interface operator(.eq.)
   module procedure compare_rules
  end interface operator(.eq.)
!
  type(confcons_rule), pointer :: confcons_rules(:), confcons_rules_bkp(:)
  integer, parameter :: rules_expand = 10 ! reallocation size increment
  integer :: num_rules=0
!
  logical :: confcons_initialized=.false., qprint=.true., qfix=.false.
!
  character(len=80), pointer :: confcons_rules_data(:)
  character(len=80), target :: confcons_rules_charmm22(5)=&
& [character(len=80) :: &
& 'ARG 1 CD 1 CZ    2 1 NH1 NH2         2 HH11 HH21 HH12 HH22', &
& 'ASP 1 CA 1 CG    2 1 OD1 OD2         0', &
& 'GLU 1 CB 1 CD    2 1 OE1 OE2         0', &
& 'PHE 1 CA 0       2 2 CD1 CD2 CE2 CE1 2 HD1 HD2 HE1 HE2', &
& 'TYR 1 CA 0       2 2 CD1 CD2 CE2 CE1 2 HD1 HD2 HE1 HE2' &
&]
! other rules
  character(len=80), target :: hydrogen_confcons_rules_charmm22(17)=&
& [character(len=80) :: &
& 'ALA 1 HA 1 CB   3 1 HB1  HB2  HB3  0', &
!
& 'ARG 1 NE 1 NH1  2 1 HH11 HH12      0', &
& 'ARG 1 NE 1 NH2  2 1 HH21 HH22      0', &
!
& 'ASN 1 CB 1 ND2  2 1 HD21 HD22      0', &
& 'GLN 1 CG 1 NE2  2 1 HE21 HE22      0', &
!
& 'ILE 1 CB 1 CD   3 1 HD1  HD2  HD3  0', &
& 'ILE 1 CA 1 CG2  3 1 HG21 HG22 HG23 0', &
!
& 'LEU 1 CB 1 CD1  3 1 HD11 HD12 HD13 0', &
& 'LEU 1 CB 1 CD2  3 1 HD21 HD22 HD23 0', &
!
& 'LYS 1 CD 1 NZ   3 1 HZ1  HZ2  HZ3  0', &
& 'MET 1 CG 1 CE   3 1 HE1  HE2  HE3  0', &
& 'THR 1 CA 1 CG2  3 1 HG21 HG22 HG23 0', &
!
& 'VAL 1 CA 1 CG1  3 1 HG11 HG12 HG13 0', &
& 'VAL 1 CA 1 CG2  3 1 HG21 HG22 HG23 0', &
!
& 'ALAD 1 OL 1 CL 3 1 HL1 HL2 HL3 0', &
& 'ALAD 1 HR 1 CR 3 1 HR1 HR2 HR3 0', &
& 'ALAD 1 HA 1 CB 3 1 HB1 HB2 HB3 0' &
&]
!
! subroutines listing
  public confcons_main
  private confcons_init
  private confcons_finalize
  public confcons_check
  private confcons_read_rules
  private compare_rules
  private confcons_rule_finalize
  private confcons_rules_finalize
  private confcons_rules_initialize
!**CHARMM_ONLY**! public ifindc
!
 contains
!!============================================
 function compare_rules(r1,r2) result(q)
 type(confcons_rule), intent(in) :: r1, r2
 logical :: q
 q=.false.
 if (r1%resname.ne.r2%resname) return
 if (.not.(associated(r1%test_atoms).eqv.associated(r2%test_atoms))) return
 if (.not.(associated(r1%orient_atoms).eqv.associated(r2%orient_atoms))) return
 if (.not.(associated(r1%orient_atom_permutations).eqv.associated(r2%orient_atom_permutations))) return
 if (.not.(associated(r1%other_atom_permutations).eqv.associated(r2%other_atom_permutations))) return
!
 if (size(r1%test_atoms).ne.size(r2%test_atoms)) return
 if (size(r1%orient_atoms).ne.size(r2%orient_atoms)) return
 if (size(r1%orient_atom_permutations).ne.size(r2%orient_atom_permutations)) return
 if (size(r1%other_atom_permutations).ne.size(r2%other_atom_permutations)) return
!
 if (any(r1%test_atoms.ne.r2%test_atoms)) return
 if (any(r1%orient_atoms.ne.r2%orient_atoms)) return
 if (any(r1%orient_atom_permutations.ne.r2%orient_atom_permutations)) return
 if (any(r1%other_atom_permutations.ne.r2%other_atom_permutations)) return
!
 q=.true. ! if we are here, then all tests passed
 end function compare_rules
!============================================
 subroutine confcons_rule_finalize(r)
 type(confcons_rule) :: r
 if (associated(r%test_atoms)) deallocate(r%test_atoms)
 if (associated(r%orient_atoms)) deallocate(r%orient_atoms)
 if (associated(r%orient_atom_permutations)) deallocate(r%orient_atom_permutations)
 if (associated(r%other_atom_permutations)) deallocate(r%other_atom_permutations)
 end subroutine confcons_rule_finalize
!============================================
 subroutine confcons_rules_finalize(rules)
 type(confcons_rule), pointer :: rules(:), r
 integer :: i
 do i=1,size(rules)
  r=>rules(i)
  call confcons_rule_finalize(r)
 enddo
 end subroutine confcons_rules_finalize
!============================================
 subroutine confcons_rules_initialize(rules)
 type(confcons_rule), pointer :: rules(:), r
 integer :: i
 do i=1,size(rules)
  r=>rules(i)
  nullify(r%test_atoms,r%orient_atoms,r%orient_atom_permutations,r%other_atom_permutations)
 enddo
 end subroutine confcons_rules_initialize
!============================================
 subroutine confcons_init()
 if (confcons_initialized) then
   call confcons_finalize()
 else
   nullify(confcons_rules, confcons_rules_bkp)
 endif
 confcons_rules_data=>confcons_rules_charmm22
 allocate(confcons_rules(rules_expand)) ; num_rules=0
 call confcons_rules_initialize(confcons_rules)
 confcons_initialized=.true.
 end subroutine confcons_init
!============================================
 subroutine confcons_finalize()
 if (confcons_initialized) then
   if (associated(confcons_rules)) then
    call confcons_rules_finalize(confcons_rules)
    deallocate(confcons_rules)
    num_rules=0
   endif
   if (associated(confcons_rules_bkp)) then
    call confcons_rules_finalize(confcons_rules_bkp)
    deallocate(confcons_rules_bkp)
   endif
 endif
 confcons_initialized=.false.
 end subroutine confcons_finalize
!============================================
 subroutine confcons_read_rules(ifile)
 use parser
 use output
!
 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
 integer, optional :: ifile
 logical :: qread
!
 character(len=21), parameter :: whoami=' CONFCONS_READ_RULES>'
 integer, parameter :: maxrulelength=200
 character(len=maxrulelength) :: line, line_bkp
 integer :: ioerr, i, j, k, l, m, n, ii, l_bkp
!
 type(confcons_rule), pointer :: r
!
 character, parameter :: comment(2) = (/'*','!'/)
!
 if (present(ifile)) then ; qread=ifile.gt.0 ; else ; qread=.false. ; endif
!
 if (qread) then
  if (qprint) then
   write(msg___,'(2A,I4)') whoami, ' READING RULES FROM UNIT ',ifile
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  endif
  do while (.true.)
! skip comments
   read(ifile,'(A)',IOSTAT=ioerr) line
!aa
! write(0,*) ioerr,line
!aa
   if (ioerr.eq.0) then
    l=len(line)
    l=min(max(0,l),len(line));line(l+1:)='';call adjustleft(line,(/' ',tab/));l=len_trim(line)
    if (any(comment.eq.line(1:1))) cycle ! skip comments
    if (l.eq.0) exit
    if (line(1:3).eq.'END') exit
! read rule
! reallocate rules array, if too small
    if (num_rules.eq.size(confcons_rules)) then ! reallocate
     if (associated(confcons_rules_bkp)) then
      call confcons_rules_finalize(confcons_rules_bkp)
      deallocate(confcons_rules_bkp)
     endif
     allocate(confcons_rules_bkp(num_rules+rules_expand))
     call confcons_rules_initialize(confcons_rules_bkp) ! nullify backup pointers
     confcons_rules_bkp(1:num_rules)=confcons_rules ! copy rules to new location (pointer assignment; data stays)
! cannot "finalize" confcons because that will destroy the data, not just the pointers
! call confcons_rules_finalize(confcons_rules) ! destroy original rules array
     deallocate(confcons_rules) ; confcons_rules=>confcons_rules_bkp ! point to new rules array
! cannot initialize below because it will nullify _the same_ pointers that confcons_rules points to
! call confcons_rules_initialize(confcons_rules_bkp) ;
     nullify(confcons_rules_bkp) ! nullify backup pointers
    endif ! end of reallocation
    r=>confcons_rules(num_rules+1)
! read rule using string processing routines
    line_bkp=line ; l_bkp=l ! save line because it is destroyed by string routines
! residue name
    r%resname=pop_string(line,l) ; l=len_trim(line)
! test atoms
    j=atoi(pop_string(line,l)) ; l=len_trim(line)
    if (j.ge.0) then ; allocate(r%test_atoms(j)) ; do i=1,j ; r%test_atoms(i)=pop_string(line,l) ; l=len_trim(line) ; enddo
    else ; ioerr = 1 ; endif
! orientation atoms
    j=atoi(pop_string(line,l)) ; l=len_trim(line)
    if (j.ge.0) then ; allocate(r%orient_atoms(j)) ; do i=1,j ; r%orient_atoms(i)=pop_string(line,l) ; l=len_trim(line) ; enddo
    else ; ioerr = 1 ; endif
! permutation pairs (part of orientation set)
    k=atoi(pop_string(line,l)) ; l=len_trim(line); ! number of atoms inside group (permutations)
    if (k.gt.0) then
     j=atoi(pop_string(line,l)) ; l=len_trim(line) ! number of orientation groups
     if (j.gt.0) then
      allocate(r%orient_atom_permutations(k,j))
      do i=1,j ; do m=1,k ; r%orient_atom_permutations(m,i)=pop_string(line,l) ; l=len_trim(line) ; enddo ; enddo
     else ; ioerr = 1 ; endif ! j>0
!
     j=atoi(pop_string(line,l)) ; l=len_trim(line) ! number of other groups that are to be permuted (not part of orientation)
     if (j.ge.0) then ! this can be zero if there are no other atoms
      allocate(r%other_atom_permutations(k,j))
      do i=1,j ; do m=1,k ; r%other_atom_permutations(m,i)=pop_string(line,l) ; l=len_trim(line) ; enddo ; enddo
     else ; ioerr = 1 ; endif ! j>0 ?
    else ; ioerr = 1 ; endif ! k>0 ?
!
    if (ioerr.eq.0) then
     num_rules=num_rules+1 ! successful read
! check for duplicate rules
     do j=1,num_rules-1
      if (r .eq. confcons_rules(j)) then
! if (any(r .eq. confcons_rules(1:num_rules))) then
       write(msg___,'(A)') 'IGNORING DUPLICATE OR CONFLICTING RULE: '//line_bkp(1:l_bkp)
       write(msg___,*)msg___;call warning(whoami, msg___(1), 0)
       num_rules=num_rules-1
       exit
      endif
     enddo
    else
     write(msg___,*)'DETECTED ERRORS IN RULES FILE.';call warning(whoami, msg___(1), 0)
     write(msg___, '(A)') 'SKIPPING LINE:'//line_bkp(1:l_bkp)
     write(msg___,*)msg___;call warning(whoami, msg___(1), 0)
     continue
    endif ! read error
   else ! end of file error
     exit
   endif
  enddo ! over all lines in the file
!
  if (qprint) then
   write(msg___,'(A,I4,A,I4)') whoami, num_rules, ' RULES READ FROM UNIT',ifile
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  endif
 else ! qread
! read internal file
! reallocate rules array, if too small
  ii=size(confcons_rules_data)+num_rules ! total number of rules
  if (ii.gt.size(confcons_rules)) then ! reallocate
   if (associated(confcons_rules_bkp)) then
    call confcons_rules_finalize(confcons_rules_bkp)
    deallocate(confcons_rules_bkp)
   endif
   allocate(confcons_rules_bkp(ii)) ! new size
   call confcons_rules_initialize(confcons_rules_bkp) ! nullify pointers inside bkp_rules
   confcons_rules_bkp(1:num_rules)=confcons_rules ! copy rules from rule array (pointer assignment; data in place)
! cannot "finalize" confcons because that will destroy the data, not just the pointers
! call confcons_rules_finalize(confcons_rules) ! destroy original rules array
   deallocate(confcons_rules) ; confcons_rules=>confcons_rules_bkp ! point to new rules array
! cannot initialize below because it will nullify _the same_ pointers that confcons_rules points to
! call confcons_rules_initialize(confcons_rules_bkp) ;
   nullify(confcons_rules_bkp) ! nullify backup rules array
  endif ! end of reallocation
!
  n=num_rules
  do ii=1,size(confcons_rules_data)
   read(confcons_rules_data(ii),'(A)') line
   l=len(line)
   l=min(max(0,l),len(line));line(l+1:)='';call adjustleft(line,(/' ',tab/));l=len_trim(line)
   n=n+1
   r=>confcons_rules(n)
! read rule using string processing routines
   r%resname=pop_string(line,l) ; l=len_trim(line)
   j=atoi(pop_string(line,l)) ; l=len_trim(line) ; allocate(r%test_atoms(j)) ; do i=1,j ; r%test_atoms(i)=pop_string(line,l) ; l=len_trim(line) ; enddo
   j=atoi(pop_string(line,l)) ; l=len_trim(line) ; allocate(r%orient_atoms(j)) ; do i=1,j ; r%orient_atoms(i)=pop_string(line,l) ; l=len_trim(line) ; enddo
!
   k=atoi(pop_string(line,l)) ; l=len_trim(line) ! number of atoms inside group (permutations)
   j=atoi(pop_string(line,l)) ; l=len_trim(line) ! number of orientation groups
   allocate(r%orient_atom_permutations(k,j))
   do i=1,j ; do m=1,k ; r%orient_atom_permutations(m,i)=pop_string(line,l) ; l=len_trim(line) ; enddo ; enddo
!
   j=atoi(pop_string(line,l)) ; l=len_trim(line) ! number of other groups that are to be permuted (not part of orientation)
   allocate(r%other_atom_permutations(k,j))
   do i=1,j ; do m=1,k ; r%other_atom_permutations(m,i)=pop_string(line,l) ; l=len_trim(line) ; enddo ; enddo
!
! the user can load the same internal file more than once (e.g. HYDR)
! therefore we should chek for duplicated here also
   do j=1,n-1
    if (r .eq. confcons_rules(j)) then
! if (any(r .eq. confcons_rules(1:num_rules))) then
     write(msg___, '(A)') 'IGNORING DUPLICATE OR CONFLICTING RULE: '//confcons_rules_data(ii)
     write(msg___,*)msg___;call warning(whoami, msg___(1), 0)
     n=n-1
     exit
    endif
   enddo
!
! aa
! write(0,*) n,':'
! write(0,*) r%resname
! write(0,*) r%test_atoms
! write(0,*) r%orient_atoms
! write(0,*) r%orient_atom_permutations
! write(0,*) r%other_atom_permutations
!aa
  enddo ! internal file
  num_rules=n
 endif
end subroutine confcons_read_rules
!============================================
subroutine confcons_main(comlyn,comlen,islct)
 use parser
 use output
 use multicom_aux !**CHARMM_ONLY**! !##MULTICOM
!
 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
 character(len=*) :: comlyn
 integer :: comlen, l, j
 integer, intent(in) :: islct(:)
!
 type(confcons_rule), pointer :: r
 integer :: errnum, ifile
 integer :: ngroups, npermute
 logical :: qcheck
 character(len=15), parameter :: whoami=' CONFCONS_MAIN>'
!
 qcheck=.false. ! whether to run checker
!
 l=comlen
 999 continue
 if (l.le.1 .or. remove_tag(comlyn,'HELP',comlen).gt.0) then
  write(msg___,'(2A)') whoami, ' CONFORMATION CONSISTENCY CHECKER ( VO / KARPLUS GROUP / HARVARD U. 2008-12 )',&
  & whoami, ' ____________________________________________________________________________________',&
  & whoami, ' DESCRIPTION: FIND AND ELIMINATE APPARENT INCONSISTENCY IN THE LABELING OF',&
  & whoami, ' DESCRIPTION: ATOMS IN RESIDUES WITH SYMMETRY ACCORDING TO PRESCRIBED RULES',&
  & whoami, ' SYNTAX : coor confcons <COMMANDS> <ATOM SELECTION>',&
  & whoami, ' COMMANDS CAN BE ONE OR MORE OF THE FOLLOWING:',&
  & whoami, '   INIT  - initialize/reinitialize',&
  & whoami, '   FINA  - deallocate all arrays',&
  & whoami, '   CHECK - check structure',&
  & whoami, '   FIX   - check structure and attempt to fix errors',&
  & whoami, '   NOFX  - do not attempt to correct errors',&
  & whoami, '   NOPR  - disable output',&
  & whoami, '   PRIN  - enable output',&
  & whoami, '   HELP  - print this screen',&
  & whoami, '   RULE  - print rules that are currently defined',&
  & whoami, '   HYDR  - include hydrogen atom groups in consistency check (excluded by default)',&
  & whoami, '   READ <UNIT> <ADD> - read rules from unit (or from input stream if unit is omitted)',&
  & whoami, '    "ADD" will append rules to the existing rules',&
  & whoami, ' ATOM SELECTION IS OPTIONAL (ALL ATOMS WILL BE SELECTED BY DEFAULT)',&
  & whoami, ' ____________________________________________________________________________________'
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  return
 endif
!
 if (remove_tag(comlyn,'FINA',comlen).gt.0) then ; call confcons_finalize() ; return ; endif
 if ((remove_tag(comlyn,'INIT',comlen).gt.0).or.(.not. confcons_initialized)) then ;
  call confcons_init() ;
  qprint=ME_LOCAL.eq.0
!
  call confcons_read_rules()
 endif
!
 if (remove_tag(comlyn,'NOPR',comlen).gt.0) then
  qprint=.false.
 elseif (remove_tag(comlyn,'PRIN',comlen).gt.0) then
  qprint=ME_LOCAL.eq.0
 endif
!
 if (remove_tag(comlyn,'HYDR',comlen).gt.0) then
  confcons_rules_data=>hydrogen_confcons_rules_charmm22
  if (remove_tag(comlyn,'ADD',comlen).le.0) num_rules=0 ! overwrite rules
  call confcons_read_rules()
 endif
!
 if (remove_tag(comlyn,'READ',comlen).gt.0) then
  if (remove_tag(comlyn,'ADD',comlen).le.0) num_rules=0 ! overwrite rules
  ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
  if (ifile .eq. -1) then
   write(msg___,*)' RULES UNIT NUMBER NOT SPECIFIED. WILL NOT READ.';call warning(whoami, msg___(1), 0)
  else
   call confcons_read_rules(ifile)
  endif
 endif ! read
!
 if (remove_tag(comlyn,'RULE',comlen).gt.0) then
  write(msg___,'(2A)') whoami, ' __________________________________________________________________________'
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  write(msg___,'(2A)') whoami, ' THE FOLLOWING RULES ARE CURRENTLY DEFINED:'
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  write(msg___,'(2A)') whoami, ' __________________________________________________________________________'
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  do j=1,num_rules
   r=>confcons_rules(j)
   write(msg___,'(A," ",I3," : ",A)') whoami, j, r%resname
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   ngroups=size(r%orient_atom_permutations,2)
   npermute=size(r%orient_atom_permutations,1)
   if (ngroups.eq.1) then
    write(msg___,'(A,'//itoa(npermute)//'A)') '    ORIENTATION GROUPS      : ',  r%orient_atom_permutations
    do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   elseif (ngroups.gt.1) then
    write(msg___,'(A,'//itoa(npermute)//'A,'//itoa(ngroups-1)//'(" / ",'//itoa(npermute)//'A))') '    ORIENTATION GROUPS      : ',  r%orient_atom_permutations
    do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   endif
!
   if (size(r%orient_atoms).gt.0) then
    write(msg___,'(A,'//itoa(size(r%orient_atoms))//'A)' )             '    OTHER ORIENTATION ATOMS : ',  r%orient_atoms
    do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   endif
   if (size(r%orient_atoms).gt.0) then
    write(msg___,'(A,'//itoa(size(r%test_atoms))//'A)' )               '    TEST ATOMS              : ',  r%test_atoms
    do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   endif
!
   ngroups=size(r%other_atom_permutations,2)
   if (ngroups.eq.1) then
    write(msg___,'(A,'//itoa(npermute)//'A)') '    OTHER PERMUTATION GROUPS: ',  r%orient_atom_permutations
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   elseif (ngroups.gt.1) then
    write(msg___,'(A,'//itoa(npermute)//'A,'//itoa(ngroups-1)//'(" / ",'//itoa(npermute)//'A))') '    OTHER PERMUTATION GROUPS: ',  r%other_atom_permutations
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   endif
  enddo
  write(msg___,'(2A)') whoami, ' __________________________________________________________________________'
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 endif ! debug
!
 if (remove_tag(comlyn,'NOFX',comlen).gt.0) then
  qfix=.false.
 elseif (remove_tag(comlyn,'FIX',comlen).gt.0) then ! whether to move atoms to fix the inconsistency
  qfix=.true.
  qcheck=.true.
 endif
!
 if ((remove_tag(comlyn,'CHCK',comlen).gt.0) .or. &
 & (remove_tag(comlyn,'CHECK',comlen).gt.0) .or.&
 & (remove_tag(comlyn,'CHEC',comlen).gt.0) .or.&
 & (remove_tag(comlyn,'CHK',comlen).gt.0) .or. &
 & qcheck) &
 & errnum=confcons_check(islct)
!
 if (comlen.eq.l) then
  comlyn='HELP '//comlyn
  comlen=min(max(0,comlen),len(comlyn));comlyn(comlen+1:)='';call adjustleft(comlyn,(/' ',tab/));comlen=len_trim(comlyn)
  goto 999
 endif
!
end subroutine confcons_main
!============================================
function confcons_check(islct) result(errnum)
 use parser
 use output
 use multicom_aux
 use system, only : r, rcomp
 use psf
 use constants
!




  use bestfit, only : eig3s, RMSBestFit, rmsd, norm3, veccross3
!
 integer :: islct(:)
!
 type(confcons_rule), pointer :: rr
 integer :: i, j, k, l, m, errnum, ires, iseg
 integer :: itest, iorient, iorient_permute, iother_permute, itotal, ipermute, npermute, ngroups
 integer :: istart, iend
 integer :: i0, i1
 integer, allocatable :: ind(:)
 real*8, allocatable :: r1(:,:), r2(:,:), ow(:)
 real*8 :: u(3,3), rcom1(3), rcom2(3) ! rotation matrix, center-of-mass coordinates
 real*8 :: msd0, msd1, msd2
 logical :: flagged
 character(len=16), parameter :: whoami=' CONFCONS_CHECK>'

 integer :: ierror

!
 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
!
 errnum=0 ! number of errors
!
! other code to be written
 call mpi_bcast(r,psf_natom(),MPI_REAL,0,MPI_COMM_LOCAL,ierror)
!
! set error count
 !**CHARMM_ONLY**! call setmsi('CONFERR',errnum)
!
!
! print summary
!
 if (qprint) then
  write(msg___,'(9A)') whoami,' ________________________________________________________________________________'
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  if (qfix) then
   if (errnum.eq.1) then
    write(msg___,'(A,I5,A)') whoami, errnum, ' INCONSISTENCY WAS FOUND AND CORRECTED.'
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   elseif (errnum.gt.1) then
    write(msg___,'(A,I5,A)') whoami, errnum, ' INCONSISTENCIES WERE FOUND AND CORRECTED.'
    do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   endif
  else
   if (errnum.eq.1) then
    write(msg___,'(A,I5,A)') whoami, errnum, ' INCONSISTENCY WAS FOUND.'
    do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   elseif (errnum.gt.1) then
    write(msg___,'(A,I5,A)') whoami, errnum, ' INCONSISTENCIES WERE FOUND.'
    do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
   endif
  endif
 endif
!
end function confcons_check
!
!
!
end module confcons
!
