! **********************************************************************!
! This source file was was generated automatically from a master source !
! code tree, which may or may not be distributed with this code, !
! because it is up to the distributor, and not up to me. !
! If you edit this file (rather than the master source file) !
! your changes will be lost if another pull from the master tree occurs.!
! In case you are wondering why, this approach makes it possible for !
! me to have the same master source code interfaced with different !
! applications (some of which are written in a way that is quite far !
! from being object-oriented) at the source level. !
! **********************************************************************!
module system
 use parser
 use output
 use ch_param
 use psf
 use files
 private
 integer, public, save :: natom =-1
 integer, public, parameter :: ndim=3 ! number of dimensions
 real*8, public, pointer, save, dimension(:,:) :: r, vr, fr ! positions, velocities, forces
 real*8, public, pointer, save, dimension(:) :: m, q, radius ! mass, charges & radii (duplicated from structure & topology files)
 real*8, public, pointer, save, dimension(:) :: bfactor, occupancy
 real*8, public, save :: BondE, AngleE, DiheE, ImprE, KinE(2)! , CmapE, ElecE, VdWE, PotE, KinE, AuxE ! energy terms
 logical, public :: system_parameters_initialized=.false.
 logical, public :: system_structure_initialized=.false.
 logical, public :: system_coordinates_initialized=.false.
 logical, public :: system_velocities_initialized=.false.
 logical, public :: system_ok=.false.
 public system_read_parameters
 public system_list_parameters
 public system_read_structure
 public system_check
 public system_center ! bring molecule to the center (of geometry, since masses are not, generally, available)
 public system_get_vw_radius ! populate radius array (above) from nonbonded parameters
 public system_printe
 public system_read_coordinates
 public system_write_dcd
 public system_write_coordinates
 public system_read_velocities
 public system_init_velocities
 public system_compute
 public system_done
 contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_read_parameters(filename) ! read parameter file
 implicit none
 character(len=*) :: filename
 character(len=80) :: fname
 character(len=10) :: paramtype
 character(len=22) , parameter :: whoami='SYSTEM_READ_PARAMETERS'
 integer :: flen
 integer :: fid=-1
!
 fname=filename
 call adjustleft(fname)
 flen=len_trim(fname)
 if (flen.gt.0) then
   call files_open(fid,fname(1:flen),'FORMATTED','READ')
   if (fid.le.0) then
    call warning(whoami, 'Cannot open parameter file. Abort.',-1)
    return
   endif
 else
  call warning(whoami, 'Parameter file name not specified. Abort.',-1)
  return
 endif
!
! call parser
 paramtype=getval('paramtype')
 call toupper(paramtype)
  select case(paramtype)
   case('CHARMM','XPLOR')
    call message(whoami, 'Reading parameters in CHARMM format.')
    call parse_ch_param(fid)
   case default
    call error(whoami, 'Unknown parameter file format. Abort.',-1)
    return
  end select
!
 system_parameters_initialized=.true.
 call files_close(fid)
!
 end subroutine system_read_parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_list_parameters()
 implicit none
 character(len=10) :: paramtype
 character(len=22) , parameter :: whoami='SYSTEM_LIST_PARAMETERS'
!
 paramtype=getval('paramtype')
 call toupper(paramtype)
  select case(paramtype)
   case('CHARMM','XPLOR')
    call list_ch_params()
   case default
    call error(whoami, 'Unknown parameter file format. Abort.',-1)
  end select
!
 end subroutine system_list_parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_read_structure(filename)
 use files
 implicit none
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: fname
 character(len=10) :: paramtype
 character(len=21) , parameter :: whoami='SYSTEM_READ_STRUCTURE'
 integer :: flen
 integer :: fid=-1
!
 fname=filename
 call adjustleft(fname)
 flen=len_trim(fname)
 if (flen.gt.0) then
   call files_open(fid,fname(1:flen),'FORMATTED','READ')
   if (fid.le.0) then
    call warning(whoami, 'Cannot open structure file. Abort.',-1)
    return
   endif
 else
  call error(whoami, 'Structure file name not specified. Abort.',-1)
  return
 endif
! call parser
 paramtype=getval_nocase('paramtype')
 call toupper(paramtype)
  select case(paramtype)
   case('CHARMM');
    call message(whoami, 'Reading structure in CHARMM format (PSF) from file '//fname(1:flen)//'.')
    call psf_read(fid)
    call psf_info()
   case('XPLOR');
    call message(whoami, 'Reading structure in X-PLOR format (PSF) from file '//fname(1:flen)//'.')
    call psf_read(fid,.true.)
    call psf_info()
! call psf_print()
   case default;
    call error(whoami, 'Unknown structure file format. Abort.',-1)
    return
  end select
!
 call files_close(fid)
!
 system_structure_initialized=.true.
!
!allocate local arrays & copy mass and charge from structure info
 end subroutine system_read_structure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_read_coordinates(filename)
 use charmmio
 use freeio
 use pdbio
 use files
!
 implicit none
!



!
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: fname
 character(len=10) :: coortype
 character(len=23), parameter :: whoami='SYSTEM_READ_COORDINATES'
 integer :: flen
 integer :: fid=-1
 real*8 :: unknown = -99999999
!
 integer*4 :: me
 me=0
!
 fname=filename
 call adjustleft(fname)
 flen=len_trim(fname)
 if (flen.gt.0) then
  if (me.le.0) then ! if me=-1, assume valid serial mode
   call files_open(fid,fname(1:flen),'FORMATTED','READ')
   if (fid.le.0) then
    call warning(whoami, 'Cannot open coordinate file. Abort.',-1)
    return
   endif
  endif
 else
  call warning(whoami, 'Coordinate file name not specified. Abort.',-1)
  return
 endif
!
 if (.not.system_structure_initialized) then
  call error(whoami, 'Structure not initialized. Cannot proceed.',-1)
  return
 endif
!
!allocate local arrays & copy mass and charge from structure info
!
 natom=atoms%last
 allocate(r(ndim,natom), vr(ndim, natom), fr(ndim, natom), m(natom), q(natom), radius(natom), bfactor(natom), occupancy(natom))
 m=atoms%mass(1:natom)
 q=atoms%charge(1:natom)
! initialize coordinates
 r=unknown; vr=unknown; fr=0d0
 occupancy=unknown;
 bfactor=unknown;
!
! call parser
 coortype=getval('coortype')
 call toupper(coortype)
  select case(coortype)
   case('CHARMM')
    call message(whoami, 'Reading coordinates in CHARMM format from file "'//fname(1:flen)//'".')
    if (me.le.0) call ch_coor_read(fid,r)
   case('ATOMID')
    call message(whoami, 'Reading coordinates in free format by atomid from file "'//fname(1:flen)//'".')
    if (me.le.0) call atomid_coor_read(fid,r)
   case('PDB')
    call message(whoami, 'Reading coordinates in PDB format from file "'//fname(1:flen)//'".')
    if (me.le.0) call PDB_read(fid,r,occupancy,bfactor)
   case('PQR')
    call message(whoami, 'Reading coordinates in PQR format from file "'//fname(1:flen)//'".')
    if (me.le.0) call PQR_read(fid,r,q,radius)
    call warning(whoami, 'PQR format lacks SEGID identifier. This may cause coordinate errors.',0)


! call coor_check()
   case default
    call error(whoami, 'UNKNOWN COORDINATE FILE FORMAT. ABORT.',-1)
    return
  end select
!
 if (me.le.0) call files_close(fid)
!
 system_coordinates_initialized=.true.
!
 end subroutine system_read_coordinates
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_write_coordinates(filename, which)
 use charmmio
 use freeio
 use parser



 implicit none
!
 character(len=*) :: filename, which
 character(len=len(which)) :: which2
 character(len=80) :: fname, tag
 character(len=10) :: coortype
 character(len=24) , parameter :: whoami='SYSTEM_WRITE_COORDINATES'
 integer :: flen
 integer :: fid=100, l
 real*8, pointer :: rout(:,:)
!
 integer*4 :: me
 me=0
!
 which2=which
 call toupper(which2)
 select case(which2)
  case('COOR','CORD','C','COORDINATES','COORD'); tag='coordinates'; rout=>r
  case('VEL','V','VELO','VELOCITIES'); tag='velocities'; rout=>vr
  case('FC','F','FORCE','FORCES'); tag='forces'; rout=>fr
  case default;
   call error(whoami, 'INVALID ARRAY REQUESTED FOR OUTPUT. ABORT.',0)
   return
 end select
 l=len_trim(tag)
!
 fname=filename
 call adjustleft(fname)
 flen=len_trim(fname)
 if (flen.gt.0) then
  if (me.le.0) open(fid, file=fname(1:flen), status='UNKNOWN', form='FORMATTED')
 else
  call error(whoami, 'FILE NAME NOT SPECIFIED. ABORT.',-1)
  return
 endif
!
 if (.not.system_structure_initialized) then
  call error(whoami, 'STRUCTURE NOT INITIALIZED. CANNOT PROCEED.',-1)
  return
 endif
!
! call parser
 coortype=getval('coortype')
 call toupper(coortype)
  select case(coortype)
   case('CHARMM')
    call message(whoami, 'Writing '//tag(1:l)//' in CHARMM format to file "'//fname(1:flen)//'".')
    if (me.le.0) call ch_coor_write(fid,rout)
   case('ATOMID')
    call message(whoami, 'Writing '//tag(1:l)//' in free format by atomid to file "'//fname(1:flen)//'".')
! if (me.eq.0) call atomid_coor_write(fid,rout)
   case default
    call error(whoami, 'UNKNOWN FILE FORMAT. ABORT.',-1)
    return
  end select
!
 close(fid)
!
 end subroutine system_write_coordinates
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_write_dcd(fid,addheader)
 use charmmio


 implicit none
 integer :: fid
 logical :: addheader
 character(len=16), parameter :: whoami='SYSTEM_WRITE_DCD'
!
 integer*4 :: me
 me=0
!
 if (me.le.0) call dcd_write(fid,r,addheader)
 end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_read_velocities(fname)
 implicit none
 character(len=*) :: fname
 end subroutine system_read_velocities
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_init_velocities()
 use rng
 use output
 use constants
 implicit none
 character(len=22), parameter :: whoami='SYSTEM_INIT_VELOCITIES'
 character(len=20) :: keyword
 real*8 :: itemp
 integer :: i, l
 integer :: c=9 ! random channel
!
 if (.not.system_structure_initialized) then ! masses may not be available
  call error(whoami, 'STRUCTURE NOT INITIALIZED. CANNOT PROCEED.',-1)
  return
 endif
! call parser
 if (existtag('init_temp')) then
  keyword=getval('init_temp')
 elseif (existtag('temperature')) then
  keyword=getval('temperature')
 else
  call error(whoami, 'INITIALIZATION TEMPERATURE NOT SPECIFIED. ABORT.',0)
  return
 endif
 itemp=atof(keyword)
 l=len_trim(keyword)
!
 call message(whoami, 'Initializing velocities at temperature '//keyword(1:l)//' K.')
 call randomg_vector(vr, ndim*natom, c)
 itemp=sqrt(kboltzmann*itemp)
 do i=1,3
  vr(i,:)=vr(i,:)/sqrt(m(:))*itemp
 enddo
!
 system_velocities_initialized=.true.
 end subroutine system_init_velocities
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! NOTE: need to take into account generic params (X)
 subroutine system_check()
 use ivector
 use tlist
 use psfatom, only: unknown
 implicit none
 character(len=8) :: a1, a2, a3, a4
 integer :: i, j, i1, i2, i3, i4, j1, j2, j3, j4, k, l
 integer, pointer :: jj(:)
 character(len=12) , parameter :: whoami='SYSTEM_CHECK'
!
 if (.not.system_parameters_initialized) then
   call error(whoami, 'Parameters not initialized. Cannot proceed.',-1)
  return
 endif
!
 if (.not.system_structure_initialized) then
   call error(whoami, 'Structure not initialized. Cannot proceed.',-1)
  return
 endif
!
! first, look up the atom type in the topology file and store in the structure
 do i=1,atoms%last ! structure
! see if the atom type is unknown
  if (atoms%type(i).eq.unknown) then
   i1=atoms%typeid(i)
   do j=1, tpar%last ! param
    if (tpar%typeid(j).eq.i1) then
     atoms%type(i)=tpar%a1(j)
     exit
    endif
   enddo
   if (atoms%type(i).eq.unknown) then
    call error(whoami, 'Invalid type_id '//itoa(i1)//' (Missing parameters?)',0)
    return
   endif
  endif
 enddo
!
! for each bond defined in blist, look for the corresponding parameter entry
 call message(whoami, 'Gathering bond parameters.')
 do i=1, blist%last
  j1=blist%ind(1,i); j2=blist%ind(2,i) ! bond indices
! find atom types that correspond to the type id's
  a1='';
  a2='';
  do j=1, atoms%last
   if (atoms%atomid(j).eq.j1) a1=atoms%type(j)
   if (atoms%atomid(j).eq.j2) a2=atoms%type(j)
  enddo
!
  if (a1.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(j1)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a2.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(j2)//' (THIS IS STRANGE)',0)
   return
  endif
! now look the bond parameter list for a match to a1 -- a2
  j=getbpar_ind(a1,a2)
  if (j.lt.0) return
! set bond parameter index
  blist%ind(3,i)=j
 enddo
! done with bonds
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 call message(whoami, 'Gathering angle parameters.')
 do i=1, alist%last
  j1=alist%ind(1,i); j2=alist%ind(2,i); j3=alist%ind(3,i); ! angle indices
! find atom types that correspond to the type id's
  a1='';
  a2='';
  a3='';
  do j=1, atoms%last
   if ( atoms%atomid(j).eq. j1 ) a1=atoms%type(j)
   if ( atoms%atomid(j).eq. j2 ) a2=atoms%type(j)
   if ( atoms%atomid(j).eq. j3 ) a3=atoms%type(j)
  enddo
!
  if (a1.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(j1)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a2.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(j2)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a3.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(j2)//' (THIS IS STRANGE)',0)
   return
  endif
! now look up the angle parameter list for a match to a1 -- a2 -- a3
  j=getapar_ind(a1,a2,a3)
  if (j.lt.0) return
! set angle parameter index
  alist%ind(4,i)=j
 enddo
! done with angles
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 call message(whoami, 'Gathering dihedral angle parameters.')
 do i=1, dlist%last
  i1=dlist%ind(1,i); i2=dlist%ind(2,i); i3=dlist%ind(3,i); i4=dlist%ind(4,i);! dihedral indices
! find atom types that correspond to the type id's
  a1=''; a2=''; a3=''; a4='';
  do j=1, atoms%last
   if ( atoms%atomid(j).eq. i1 ) a1=atoms%type(j)
   if ( atoms%atomid(j).eq. i2 ) a2=atoms%type(j)
   if ( atoms%atomid(j).eq. i3 ) a3=atoms%type(j)
   if ( atoms%atomid(j).eq. i4 ) a4=atoms%type(j)
  enddo
!
  if (a1.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i1)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a2.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i2)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a3.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i3)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a4.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i4)//' (THIS IS STRANGE)',0)
   return
  endif
! now look up the dihedral parameter list for a match to a1 -- a2 -- a3 -- a4
! we may have several dihedral parameters with different multiplicities -- get them all in an array:
  jj=>getdpar_ind(a1,a2,a3,a4)
  l=jj(1) ! first element is the number of following dihedral entries
  if (l.lt.1) return
! set dihedral parameter index
  dlist%ind(5,i)=jj(2)
! if additional entries are present, add additional dihedral lines to dlist:
  do k=3,l+1
   j=toplist_add(dlist,5,(/ i1, i2, i3, i4, jj(k) /)) ! dihedral entry
  enddo
  deallocate(jj)
 enddo
! done with dihedrals
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 call message(whoami, 'Gathering improper dihedral angle parameters.')
 do i=1, ilist%last
  i1=ilist%ind(1,i); i2=ilist%ind(2,i); i3=ilist%ind(3,i); i4=ilist%ind(4,i);! dihedral indices
! find atom types that correspond to the type id's
  a1=''; a2=''; a3=''; a4='';
  do j=1, atoms%last
   if ( atoms%atomid(j).eq. i1 ) a1=atoms%type(j)
   if ( atoms%atomid(j).eq. i2 ) a2=atoms%type(j)
   if ( atoms%atomid(j).eq. i3 ) a3=atoms%type(j)
   if ( atoms%atomid(j).eq. i4 ) a4=atoms%type(j)
  enddo
!
  if (a1.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i1)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a2.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i2)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a3.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i3)//' (THIS IS STRANGE)',0)
   return
  endif
  if (a4.eq.'') then
   call error(whoami, 'NO TYPE FOR ATOM '//itoa(i4)//' (THIS IS STRANGE)',0)
   return
  endif
! now look up the dihedral parameter list for a match to a1 -- a2 -- a3 -- a4
! we may have several dihedral parameters with different multiplicities -- get them all in an array:
  jj=>getipar_ind(a1,a2,a3,a4)
  l=jj(1) ! first element is the number of following dihedral entries
  if (l.lt.1) return
! set dihedral parameter index
  ilist%ind(5,i)=jj(2)
! if additional entries are present, add additional dihedral lines to dlist (note, in CHARMM, multiple impropers are illegal)
  do k=3,l+1
   j=toplist_add(ilist,5,(/ i1, i2, i3, i4, jj(k) /)) ! dihedral entry
  enddo
  deallocate(jj)
 enddo
! done with improper dihedrals
!
  system_ok=.true.
 end subroutine system_check
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_get_vw_radius() ! set radius to the VW radii in parameter file
 use psfatom, only: unknown
 implicit none
 integer :: i, j, i1
!
 character(len=20) , parameter :: whoami='SYSTEM_GET_VW_RADIUS'
!
 if (.not.system_parameters_initialized) then
   call error(whoami, 'Parameters not initialized. Cannot proceed.',-1)
  return
 endif
!
 if (.not.system_structure_initialized) then
   call error(whoami, 'Structure not initialized. Cannot proceed.',-1)
  return
 endif
!
 if (.not.system_coordinates_initialized) then
   call error(whoami, 'Coordinates not initialized. Cannot proceed.',-1)
  return
 endif
!
! note code duplication from system_check()
 do i=1,atoms%last ! structure
! see if the atom type is unknown
  if (atoms%type(i).eq.unknown) then
   i1=atoms%typeid(i)
   do j=1, tpar%last ! param
    if (tpar%typeid(j).eq.i1) then
     atoms%type(i)=tpar%a1(j)
     exit
    endif
   enddo
   if (atoms%type(i).eq.unknown) then
    call error(whoami, 'Invalid type_id '//itoa(i1)//' (Missing parameters?)',0)
    return
   endif
  endif
! now can copy the radius
  radius(i)=tpar%rmino2(j) ! note: rmin_over_two is not necessarily the optimal radius
 enddo
!
 end subroutine system_get_vw_radius
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_compute(forces)
 implicit none
 logical, optional :: forces
 logical :: f
 character(len=14) , parameter :: whoami='SYSTEM_COMPUTE'
!
 include 'interface.h'
!
 if (.not.system_ok) then
   call error(whoami, 'SYSTEM NOT INITIALIZED. ABORT.',0)
  return
 endif
!
 if (present(forces)) then
  f=forces
 else
  f=.true.
 endif
!
 fr=0d0
! (1) -- bonded interactions
 call compute_bonds3(BondE, blist, bpar, r, fr, f)
 call compute_angles3(AngleE, alist, apar, r, fr, f)
 call compute_dihes3(diheE, dlist, dpar, r, fr, f)
 call compute_dihes3(imprE, ilist, ipar, r, fr, f)
! (1i) CMAP
! (2) nonbonded interactions
 end subroutine system_compute
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_done()
 implicit none
 character(len=10) :: paramtype
 system_ok=.false.
!
 if (system_parameters_initialized) then
  system_parameters_initialized=.false.
! call parser
  if (existtag('paramtype')) then
   paramtype=getval('paramtype')
   call toupper(paramtype)
   select case(paramtype)
    case('CHARMM','XPLOR')
     call ch_param_done()
   end select
  endif
 endif ! initialized
! structure
 if (system_structure_initialized) then
  system_structure_initialized=.false.
! call parser
  if (existtag('paramtype')) then
   paramtype=getval('paramtype')
   call toupper(paramtype)
   select case(paramtype)
    case('CHARMM','XPLOR')
     call psf_done()
   end select
  endif
 endif ! initialized
! coordinates
 natom=-1
 deallocate(r,vr,fr,m,q,radius,bfactor,occupancy)
 system_coordinates_initialized=.false.
 system_velocities_initialized=.false.
!
 end subroutine system_done
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine system_printe()
 use output, only: fout, message
 use parser, only: tab
 use stats


 implicit none
 character(len=13) , parameter :: whoami='SYSTEM_PRINTE'
!
 integer*4 :: me
 me=0
!
 KinE=calcKE(vr,m)
 if (me.eq.0) then
  call message(whoami, '');
  write(fout,'(A)') tab//'=================== ENERGY ==================';
  write(fout,'(1(A,F10.5))') tab//'TotalE: ', BondE+AngleE+DiheE+ImprE+KinE(1) 
  write(fout,'(2(A,F10.5))') tab//'BondE: ', BondE, tab//'AngleE: ', AngleE 
  write(fout,'(2(A,F10.5))') tab//'DiheE: ', DiheE, tab//' ImprE: ', ImprE 
  write(fout,'(2(A,F10.5))') tab//'KinE: ', KinE(1), tab//'Temp: ', KinE(2) 
  write(fout,'(A)') tab//'=============================================';
 endif
!
 end subroutine system_printe
end module system
