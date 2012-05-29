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
module multicom_aux
!**CHARMM_ONLY**!##IF MULTICOM
! Basic communicator scheme
!
! MPI_COMM_LOCAL - main communicator ; in general,
! not the same as MPI_COMM_WORLD
! SIZE_LOCAL - number of processors in corresponding communicator
! ME_LOCAL - rank in this communicator
!
! MPI_COMM_GLOBAL - global communicator (involves all nodes)
! SIZE_GLOBAL - number of processors in corresponding communicator
! ME_GLOBAL - rank in communicator
!
  use mpi ! define constants MPI_COMM_NULL, MPI_UNDEFINED
 integer, save :: MPI_COMM_LOCAL =MPI_COMM_NULL, ME_LOCAL=MPI_UNDEFINED, SIZE_LOCAL=MPI_UNDEFINED
 integer, save :: MPI_COMM_GLOBAL =MPI_COMM_NULL, ME_GLOBAL=MPI_UNDEFINED, SIZE_GLOBAL=MPI_UNDEFINED
!!**CHARMM_ONLY**!##ENSEMBLE
! Additional communicators
  integer, save :: MPI_COMM_ENSBL =MPI_COMM_NULL, ME_ENSBL=MPI_UNDEFINED, SIZE_ENSBL=MPI_UNDEFINED !!**CHARMM_ONLY**!##ENSEMBLE
  integer, save :: MPI_COMM_STRNG =MPI_COMM_NULL, ME_STRNG=MPI_UNDEFINED, SIZE_STRNG=MPI_UNDEFINED !!**CHARMM_ONLY**!##STRINGM
  integer, save :: MPI_COMM_DMOL =MPI_COMM_NULL, ME_DMOL=MPI_UNDEFINED, SIZE_DMOL=MPI_UNDEFINED


!
!**CHARMM_ONLY**!##ENDIF
!
end module multicom_aux
