/*COORDINATES AND MASSES:*/
/*
#ifdef __IMPNONE
#undef __IMPNONE
#endif
#define __IMPNONE
*/
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
! CV_TYPES.MOD
!
! LISTING OF COLLECTIVE VARIABLE TYPES
!
! naming convention: type name is the same as the corresponding module name
! without the `cv_' prefix
!**CHARMM_ONLY**!##IF STRINGM
      module cv_types
       integer, parameter :: &
! posi_x=1, posi_y=2, posi_z=3, ! position (obsolete)
     & posi_com_x=4, posi_com_y=5, posi_com_z=6, & ! com-postion
     & dist_com=7, & ! distance
     & angle_com=8, & ! angle
     & dihe_com=9, & ! dihedral angle
     & anglvec=10, & ! angle between two vectors (possibly in different frames)
     & qcomp_1=11, qcomp_2=12, qcomp_3=13, qcomp_4=14, & ! components of orientation quaternion
     & rmsd=15, & ! RMS distance between a set of atoms
     & drmsd=16, & ! difference in the RMSDs from two structures
     & proj=17, & ! projection onto the vector connecting two reference structures; simulation structure used for orientation;
     & cvrms=99 ! rtmd-like variable (compatibility limited to steered dynamics as of 7.2010): z=sqrt( 1/N sum^N_i (z_i - z^0_i)^2 )
      end module cv_types
!
!**CHARMM_ONLY**!##ENDIF
