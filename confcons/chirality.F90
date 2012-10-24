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
module chirality
 !**CHARMM_ONLY**! __DEP_KINDS
 implicit none
 private
!
 !
 ! in this module :
 ! chirality is determined using an dihedral angle
 ! rule specified in a `rules' table (below) or in an
 ! externally provided file (for added flexibility:
 ! it is clear that different force fields and/or different
 ! molecules will have specific corresponding rules);
 ! the rules table lists the residue name,
 ! the corresponding dihedral, the threshold in degrees for
 ! passing the test, the direction in which
 ! the threshold can be exceeded, the atom name whose position
 ! can be reflected around a symmetry plane (normally, a proton,
 ! because nothing is attached to it except for the parent heavy atom)
 ! and the atom which lies in the symmetry plane whose distance
 ! vector to the previous atom is normal to the symmetry plane
 ! (e.g. the parent atom of the hydrogen)
!
!
  type dihedral_rule
   character(len=8) :: resname, aname(6)
   real*8 :: threshold
   integer :: sgn
  end type dihedral_rule
!
  interface operator(.eq.)
    module procedure compare_rules
  end interface operator(.eq.)
!
  type(dihedral_rule), pointer :: dihedral_rules(:), dihedral_rules_bkp(:)
  integer, parameter :: rules_expand = 10 ! reallocation size increment
  integer :: num_rules=0
!
  logical :: chirality_initialized=.false., qprint=.true., qflip=.false.
!
  character(len=80), pointer :: dihedral_rules_data(:)
  character(len=80), target :: dihedral_rules_charmm22(23)=&
& [character(len=80) :: &
& 'ALA N C CB HA 0. -1 HA CA',& ! the dihedral specified should be less than 0
& 'ARG N C CB HA 0. -1 HA CA',&
& 'ASN N C CB HA 0. -1 HA CA',&
& 'ASP N C CB HA 0. -1 HA CA',&
& 'CYS N C CB HA 0. -1 HA CA',&
& 'GLN N C CB HA 0. -1 HA CA',&
& 'GLU N C CB HA 0. -1 HA CA',&
& 'HIS N C CB HA 0. -1 HA CA',&
& 'HSC N C CB HA 0. -1 HA CA',&
& 'HSD N C CB HA 0. -1 HA CA',&
& 'HSE N C CB HA 0. -1 HA CA',&
& 'HSP N C CB HA 0. -1 HA CA',&
& 'ILE N C CB HA 0. -1 HA CA',&
& 'LEU N C CB HA 0. -1 HA CA',&
& 'LYS N C CB HA 0. -1 HA CA',&
& 'MET N C CB HA 0. -1 HA CA',&
& 'PHE N C CB HA 0. -1 HA CA',&
& 'PRO N C CB HA 0. -1 HA CA',&
& 'SER N C CB HA 0. -1 HA CA',&
& 'THR N C CB HA 0. -1 HA CA',&
& 'THR OG1 CA CG2 HB 0. 1 HB CB',& ! T has a chiral center at the CB atom
& 'TYR N C CB HA 0. -1 HA CA',&
& 'VAL N C CB HA 0. -1 HA CA'&
&]
!
  character(len=80), target :: pseudo_dihedral_rules_charmm22(42)=& ! these rules enforce pseudo-chirality, which concerns equivalent atoms
& [character(len=80) :: &
& 'ALA CA HB3 HB2 HB1 0. 2 HB1 HB2',& ! proton ordering within methyl groups
& 'ASN CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'ASP CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'ARG CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'ARG CB CD HG2 HG1 0. 2 HG1 HG2',&
& 'ARG CG NE HD2 HD1 0. 2 HD1 HD2',&
& 'CYS CA SG HB2 HB1 0. 2 HB1 HB2',&
& 'GLN CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'GLN CB CD HG2 HG1 0. 2 HG1 HG2',&
& 'GLU CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'GLU CB CD HG2 HG1 0. 2 HG1 HG2',&
& 'GLY N  C  HA2 HA1 0. 2 HA1 HA2',&
& 'HIS CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'HSD CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'HSE CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'ILE CA CG1 CG2 HB 0. -1 HB CB',&
& 'ILE CB HG23 HG22 HG21 0. 2 HG21 HG22',&
& 'ILE CG1 HD3 HD2 HD1 0. 2 HD1 HD2',&
& 'ILE CB CD HG12 HG11 0. 2 HG11 HG12',&
& 'LEU CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'LEU CB CD2 CD1 HG 0. -1 HG CG',&
& 'LEU CG HD13 HD12 HD11 0. 2 HD11 HD12',&
& 'LEU CG HD23 HD22 HD21 0. 2 HD21 HD22',&
& 'LYS CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'LYS CB CD HG2 HG1 0. 2 HG1 HG2',&
& 'LYS CG CE HD2 HD1 0. 2 HD1 HD2',&
& 'LYS CE HZ3 HZ2 HZ1 0. 2 HZ1 HZ2',&
& 'LYS CD NZ HE2 HE1 0. 2 HE1 HE2',&
& 'MET CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'MET CB SD HG2 HG1 0. 2 HG1 HG2',&
& 'MET SD HE3 HE2 HE1 0. 2 HE1 HE2',&
& 'PHE CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'PRO CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'PRO CB CD HG2 HG1 0. -2 HG1 HG2',&
& 'PRO N  CG HD2 HD1 0. 2 HD1 HD2',&
& 'SER CA OG HB2 HB1 0. 2 HB1 HB2',&
& 'THR CB HG23 HG22 HG21 0. 2 HG21 HG22',&
& 'TRP CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'TYR CA CG HB2 HB1 0. 2 HB1 HB2',&
& 'VAL CA CG2 CG1 HB 0. -1 HB CB',&
& 'VAL CB HG13 HG12 HG11 0. 2 HG11 HG12',&
& 'VAL CB HG23 HG22 HG21 0. 2 HG21 HG22'&
&]
!
  logical, parameter :: qdebug=.false.
! subroutines listing
  public chirality_main
  private compute_dihedral
  private chirality_init
  private chirality_finalize
  public chirality_check ! I am making this public so that checks can be called from within CHARMM
                         ! it is obviously the programmer's responsibility to make sure everything is set up for this call
  private chirality_read_rules
  private compare_rules
!
 contains
!
!============================================
 function compare_rules(r1,r2)
 type(dihedral_rule), intent(in) :: r1, r2
 logical :: compare_rules
 compare_rules=r1%resname.eq.r2%resname.and.all(r1%aname(1:4).eq.r2%aname(1:4)) ! stop here VO .and.abs(r1%threshold-r2%threshold).lt.1d-6.and.r1%sgn.eq.r2%sgn
 end function compare_rules
!============================================
 subroutine chirality_main(comlyn,comlen,islct)
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
 type(dihedral_rule), pointer :: r
 integer :: errnum, ifile
 logical :: qcheck
 character(len=16), parameter :: whoami=' CHIRALITY_MAIN>'
!
 qcheck=.false. ! whether to run checker
!
 l=comlen
 999 continue
 if (l.le.1 .or. remove_tag(comlyn,'HELP',comlen).gt.0) then
  if &
& (ME_LOCAL.eq.0)&
& then
  write(msg___,'(2A)') whoami, ' CHIRALITY CHECKER ( VO / KARPLUS GROUP / HARVARD U. 2012 )'&
& ,whoami, ' _______________________________________________________________________________'&
& ,whoami, ' DESCRIPTION: FIND AND CORRECT CHIRALITY ERRORS ACCORDING'&
& ,whoami, ' DESCRIPTION: TO PRESCRIBED RULES, APPLIED TO SELECTED RESIDUES'&
& ,whoami, ' SYNTAX : coor chirality <COMMANDS> <ATOM SELECTION>'&
& ,whoami, ' COMMANDS CAN BE ONE OR MORE OF THE FOLLOWING:'&
& ,whoami, '   INIT  - initialize/reinitialize'&
& ,whoami, '   FINA  - deallocate all arrays'&
& ,whoami, '   CHECK - check structure'&
& ,whoami, '   FIX   - attempt to correct errors during checking'&
& ,whoami, '   NOFX  - do not attempt to correct errors'&
& ,whoami, '   NOPR  - disable output'&
& ,whoami, '   PRIN  - enable output'&
& ,whoami, '   HELP  - print this screen'&
& ,whoami, '   RULE  - print rules that are currently defined'&
& ,whoami, '   PSEU  - also check for pseudo-chirality errors (i.e. atom ordering)'&
& ,whoami, '   READ <UNIT> <ADD> - read rules from unit (or from input stream if unit is omitted)'&
& ,whoami, '    "ADD" will append rules to the existing rules'&
& ,whoami, ' ATOM SELECTION IS OPTIONAL (ALL ATOMS WILL BE SELECTED BY DEFAULT)'&
& ,whoami, ' _______________________________________________________________________________'
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  endif ! me==0
  return
 endif
!
 if (remove_tag(comlyn,'FINA',comlen).gt.0) then ; call chirality_finalize() ; return ; endif
 if ((remove_tag(comlyn,'INIT',comlen).gt.0).or.(.not. chirality_initialized)) then ;
  call chirality_init() ;
  qprint=ME_LOCAL.eq.0
!
  call chirality_read_rules()
 endif
!
 if (remove_tag(comlyn,'NOPR',comlen).gt.0) then
  qprint=.false.
 elseif (remove_tag(comlyn,'PRIN',comlen).gt.0) then
  qprint=ME_LOCAL.eq.0
!
 endif
!
 if (remove_tag(comlyn,'PSEU',comlen).gt.0) then
  dihedral_rules_data=>pseudo_dihedral_rules_charmm22
  if (remove_tag(comlyn,'ADD',comlen).le.0) num_rules=0 ! overwrite rules
  call chirality_read_rules()
 endif
!
 if (remove_tag(comlyn,'READ',comlen).gt.0) then
  if (remove_tag(comlyn,'ADD',comlen).le.0) num_rules=0 ! overwrite rules
  ifile=atoi(get_remove_parameter(comlyn, 'UNIT', comlen), -1)
  if (ifile .eq. -1) then
   call warning(whoami, ' RULES UNIT NUMBER NOT SPECIFIED. WILL NOT READ.', 0)
  else
   call chirality_read_rules(ifile)
  endif
 endif ! read
!
 if (remove_tag(comlyn,'RULE',comlen).gt.0) then
!
  if &




& (ME_LOCAL.eq.0)&

& then
!
  write(msg___,'(2A)') whoami, ' __________________________________________________________________________';
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  write(msg___,'(2A)') whoami, ' THE FOLLOWING RULES ARE CURRENTLY DEFINED:';
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  write(msg___,'(2A)') whoami, ' __________________________________________________________________________';
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  do j=1,num_rules
   r=>dihedral_rules(j)
   write(msg___,'(A," ",I3," ",5A,G12.5,I3," ",2A)') whoami, j, r%resname, r%aname(1:4), r%threshold, r%sgn, r%aname(5:6)
   ;do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  enddo
  write(msg___,'(2A)') whoami, ' __________________________________________________________________________';
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  endif ! me==0
 endif ! RULES
!
 if (remove_tag(comlyn,'NOFX',comlen).gt.0) then
  qflip=.false.
 elseif (remove_tag(comlyn,'FIX',comlen).gt.0) then ! whether to flip the atom specified in the rule to change chirality
  qflip=.true.
 endif
!
 if ((remove_tag(comlyn,'CHCK',comlen).gt.0) .or. &
 & (remove_tag(comlyn,'CHECK',comlen).gt.0).or. &
 & (remove_tag(comlyn,'CHEC',comlen).gt.0) .or. &
 & (remove_tag(comlyn,'CHK',comlen).gt.0) .or. &
 & qcheck) &
 & errnum=chirality_check(islct)
!
 if (comlen.eq.l) then
  comlyn='HELP '//comlyn
  comlen=min(max(0,comlen),len(comlyn));comlyn(comlen+1:)='';call adjustleft(comlyn,(/' ',tab/));comlen=len_trim(comlyn)
  goto 999
 endif
 end subroutine chirality_main
!============================================
 subroutine chirality_init()
 if (chirality_initialized) then
   call chirality_finalize()
 else
   nullify(dihedral_rules, dihedral_rules_bkp)
 endif
 dihedral_rules_data=>dihedral_rules_charmm22
 allocate(dihedral_rules(rules_expand)) ; num_rules=0
 chirality_initialized=.true.
 end subroutine chirality_init
!============================================
 subroutine chirality_finalize()
 if (chirality_initialized) then
   if (associated(dihedral_rules)) deallocate(dihedral_rules)
   num_rules=0
   if (associated(dihedral_rules_bkp)) deallocate(dihedral_rules_bkp)
 endif
 chirality_initialized=.false.
 end subroutine chirality_finalize
!============================================
 subroutine chirality_read_rules(ifile)
 use parser
 use output

 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_

 integer, optional :: ifile
 logical :: qread
!
 character(len=22), parameter :: whoami=' CHIRALITY_READ_RULES>'
 integer, parameter :: maxrulelength=200
 character(len=maxrulelength) :: line
 integer :: ioerr, i, j, k, l
!
 type(dihedral_rule), pointer :: r
!
 character, parameter :: comment(2) = (/'*','!'/)
!
 if (present(ifile)) then ; qread=ifile.gt.0 ; else ; qread=.false. ; endif
!
 if (qread) then
  if (qprint) then
   write(msg___,'(2A,I4)') whoami, ' READING DIHEDRAL RULES FROM UNIT ',ifile;
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  endif
  do while (.true.)
! skip comments
   read(ifile,'(A)',IOSTAT=ioerr) line
   if (ioerr.eq.0) then
    l=len(line)
    l=min(max(0,l),len(line));line(l+1:)='';call adjustleft(line,(/' ',tab/));l=len_trim(line)
    if (any(comment.eq.line(1:1))) cycle ! skip comments
    if (l.eq.0) exit
    if (line(1:3).eq.'END') exit
! read rule
! reallocate rules array, if too small
    if (num_rules.eq.size(dihedral_rules)) then ! reallocate
     if (associated(dihedral_rules_bkp)) deallocate(dihedral_rules_bkp)
     allocate(dihedral_rules_bkp(num_rules+rules_expand))
     dihedral_rules_bkp(1:num_rules)=dihedral_rules
     deallocate(dihedral_rules) ; dihedral_rules=>dihedral_rules_bkp ; nullify(dihedral_rules_bkp)
    endif
    r=>dihedral_rules(num_rules+1)
    read(line,*,IOSTAT=ioerr) r%resname, r%aname(1:4), r%threshold, r%sgn, r%aname(5:6)
    if (ioerr.eq.0) then
     num_rules=num_rules+1 ! successful read
! check for duplicate rules
     do j=1,num_rules-1
      if (r .eq. dihedral_rules(j)) then
! if (any(r .eq. dihedral_rules(1:num_rules))) then
       write(msg___,'(A)') 'IGNORING DUPLICATE OR CONFLICTING RULE: '//line(1:l); 
       call warning(whoami, msg___(1), 0)
       num_rules=num_rules-1
       exit
      endif
     enddo
    else
     call warning(whoami, 'DETECTED ERRORS IN RULES FILE.', 0)
     write(msg___, '(A)') 'SKIPPING LINE:'//line(1:l); 
     call warning(whoami, msg___(1), 0)
     continue
    endif ! read error
   else ! end of file error
     exit
   endif
  enddo ! over all lines in the file
!
  if (qprint) then
   write(msg___,'(A,I4,A,I4)') whoami, num_rules, 'DIHEDRAL RULES READ FROM UNIT',ifile ;
   do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  endif
 else ! qread
! read internal file
! reallocate rules array, if too small
  l=size(dihedral_rules_data)+num_rules
  if (l.gt.size(dihedral_rules)) then ! reallocate
   if (associated(dihedral_rules_bkp)) deallocate(dihedral_rules_bkp)
   allocate(dihedral_rules_bkp(l))
   dihedral_rules_bkp(1:num_rules)=dihedral_rules
   deallocate(dihedral_rules) ; dihedral_rules=>dihedral_rules_bkp ; nullify(dihedral_rules_bkp)
  endif
  j=num_rules
  do i=1,size(dihedral_rules_data)
   j=j+1
   r=>dihedral_rules(j)
   read(dihedral_rules_data(i),*) &
& r%resname, &
& r%aname(1:4), &
& r%threshold, &
& r%sgn,&
& r%aname(5:6)
! the user can load the same internal file more than once (e.g. PSEU)
! therefore we should chek for duplicated here also
   do k=1,j-1
    if (r .eq. dihedral_rules(k)) then
       write(msg___,'(A)') 'IGNORING DUPLICATE OR CONFLICTING RULE: '//dihedral_rules_data(i)
       call warning(whoami, msg___(1), 0)
     j=j-1
     exit
    endif
   enddo
!
  enddo ! internal file
  num_rules=j
 endif
end subroutine chirality_read_rules
!============================================
 function chirality_check(islct_, r__) result(errnum)
 use parser
 use output
 use system, only : r, rcomp, m, bfactor, occupancy
 use psf
 use constants, only : pi
 use multicom_aux !**CHARMM_ONLY**! !##MULTICOM



!
 character(len=200) :: msg___(20)=(/'','','','','','','','','','','','','','','','','','','',''/); integer :: i_
!
 integer, optional, intent(in) :: islct_(:)
 integer, pointer :: islct(:)
 real*8, optional :: r__(:,:)
 real*8, pointer :: r_(:,:)
!
 character, parameter :: op(-1:1) = (/'>',' ','<'/) ! trick to write the correct inequality
 type(dihedral_rule), pointer :: rr
 integer :: i,j,ind(6), errnum, ires, iseg
 real*8 :: phi, d
 logical :: flagged
 integer :: ierror
!

 integer :: natom

!

!
 character(len=17), parameter :: whoami=' CHIRALITY_CHECK>'
!
 if (present(islct_)) then
  allocate(islct(size(islct_))); islct=islct_
 else
  allocate(islct(natom)) ; islct=1
 endif
!
 if (present(r__)) then
   allocate(r_(size(r__,1),size(r__,2))); r_=r__
 else
   allocate(r_(size(r(1,:)),3)) ; r_(:,1)=r(1,:) ; r_(:,2)=r(2,:) ; r_(:,3)=r(3,:)
 endif
!
!
!
 errnum=0 ! number of errors
!
! charmm-specific code
!
! other code to be written
! populate main coordinates, if needed
 if (present(r__)) then
  r__=r_
 else
  r(1,:)=r_(:,1) ; r(2,:)=r_(:,2) ; r(3,:)=r_(:,3)
 endif
!
! set error count
 !**CHARMM_ONLY**! call setmsi('CHIERR',errnum)
!
! print summary
!
 if (qprint) then
  write(msg___,'(9A)') whoami,' ________________________________________________________________________________' ;
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
  if (qflip.and.errnum.gt.0) then
   write(msg___,'(A,I5,A)') whoami, errnum, ' VIOLATIONS WERE FOUND AND CORRECTED. SYSTEM MAY NEED TO BE MINIMIZED.'
  else
   write(msg___,'(A,I5,A)') whoami, errnum, ' VIOLATIONS WERE FOUND.'
  endif
  do i_=1,size(msg___);if(msg___(i_)=='')exit;call plainmessage(msg___(i_));enddo;msg___=''
 endif
!
 if (associated(r_)) deallocate(r_)
 if (associated(islct)) deallocate(islct)
!
end function chirality_check
!============================================================================
 function compute_dihedral(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4) result(theta)
 real*8 :: theta, costh, sinth,&
& x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4,&
& dx12, dy12, dz12, dx23, dy23, dz23, dx34, dy34, dz34,&
& vx, vy, vz, vn, ux, uy, uz, un, wx, wy, wz, wn
 dx12=x2-x1; dy12=y2-y1; dz12=z2-z1;
 dx23=x3-x2; dy23=y3-y2; dz23=z3-z2;
 dx34=x4-x3; dy34=y4-y3; dz34=z4-z3;
! note:
! costh = [ (d34 x d32) . (d23 x d21) / |(d34 x d23)| |(d23 x d12)| ] =
! = cos [ (u . v) / |u| |v| ] (definition of u and v)
 ux=dy23*dz34-dy34*dz23;
 uy=dz23*dx34-dz34*dx23;
 uz=dx23*dy34-dx34*dy23;
 un=sqrt(ux*ux+uy*uy+uz*uz);
!
 vx=dy12*dz23-dy23*dz12;
 vy=dz12*dx23-dz23*dx12;
 vz=dx12*dy23-dx23*dy12;
 vn=sqrt(vx*vx+vy*vy+vz*vz);
!
 wx=dy23*vz-vy*dz23;
 wy=dz23*vx-vz*dx23;
 wz=dx23*vy-vx*dy23;
 wn=sqrt(wx*wx+wy*wy+wz*wz);
!
 if (un.eq.0d0) un=1d0
 if (vn.eq.0d0) vn=1d0
 if (wn.eq.0d0) wn=1d0
 costh=(ux*vx+uy*vy+uz*vz)/(un*vn)
 sinth=(wx*ux+wy*uy+wz*uz)/(wn*un)
 theta=atan2(sinth, costh)
!
 end function compute_dihedral
!
end module chirality
