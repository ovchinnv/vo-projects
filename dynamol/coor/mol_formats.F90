/*COORDINATES AND MASSES:*/
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
module mol_formats
 public
!
  integer, parameter :: pdb=1, charmm=2, free=3, pqr=4
  integer, parameter :: num_mol_format=4
  character(len=20), parameter :: mol_format_name(num_mol_format)=&
  & [character(len=20) :: 'PDB', 'CHARMM', 'FREE', 'PQR']
  integer, parameter :: ocol=1, bcol=2, wcol=3, xcol=4, ycol=5, zcol=6
  integer, parameter :: num_mol_cols=6
  character(len=20), parameter :: mol_col_name(num_mol_cols)=&
  & [character(len=20) :: 'OCCUPANCY', 'BETA', 'WEIGHT', 'X', 'Y', 'Z']
!
end module mol_formats
