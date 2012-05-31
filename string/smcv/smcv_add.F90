/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
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
! THIS FILE CONTAINS ROUTINES FOR ADDING NEW CV
! THERE IS A SEPARATE ROUTINE FOR EACH CV TYPE
!
!**CHARMM_ONLY**!##IF STRINGM
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_add(COMLYN,COMLEN)
      use cv_types
      use cv_common
      use cv_quaternion, only: quat_initialized, quat, quat_init
      use sm_config, only: qt_send_displ, qt_send_count, &
     & imap_displ, imap_count
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      implicit none
!
 character(len=80) :: msg___
!
      character(len=*) :: COMLYN
      integer :: COMLEN
! local
      integer i, j, k
      character(len=20) :: keyword
      character(len=10) :: whoami
      data whoami /' SMCV_ADD>'/
!
      keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn) ! directive
! COM positions
      if (( keyword(1:10).eq.'POSI_COM_X'(1:10) )) then
       call smcv_posi_com_add(comlyn, comlen, posi_com_x) ! note: the same routine for all position coordinates
      else if (( keyword(1:10).eq.'POSI_COM_Y'(1:10) )) then
       call smcv_posi_com_add(comlyn, comlen, posi_com_y)
      else if (( keyword(1:10).eq.'POSI_COM_Z'(1:10) )) then
       call smcv_posi_com_add(comlyn, comlen, posi_com_z)
! dihedral_COM
      else if (( keyword(1:8).eq.'DIHE_COM'(1:8) )) then
       call smcv_dihe_com_add(comlyn, comlen)
! angle_COM
      else if (( keyword(1:9).eq.'ANGLE_COM'(1:9) )) then
       call smcv_angle_com_add(comlyn, comlen)
! anglvec
      else if (( keyword(1:7).eq.'ANGLVEC'(1:7) )) then
       call smcv_anglvec_add(comlyn, comlen)
! distance_COM
      else if (( keyword(1:8).eq.'DIST_COM'(1:8) )) then
       call smcv_dist_com_add(comlyn, comlen)
! reference frame (not really a CV, but processed here)
      else if (( keyword(1:5).eq.'FRAME'(1:5) )) then
       call smcv_frame_add(comlyn, comlen)
! orientation quaternion
      else if (( keyword(1:10).eq.'QUATERNION'(1:10) )) then
       call smcv_quaternion_add(comlyn, comlen)
!
! (re)compute quaternion index limits (for parallelization) after each addition
!
       if (SIZE_LOCAL.gt.0) then
        if (.not.quat_initialized) call quat_init() ! make sure frames%num_frames is defined
        j=ceiling(1.0d0*quat%num_quat/SIZE_LOCAL) ! max. number of frames assigned to slave node
        k=ceiling(1.0d0*cv%amap%last/SIZE_LOCAL) ! max. number of amap indices assigned to slave node
!
        do i=1,SIZE_LOCAL
         qt_send_displ(i)=min((i-1)*j,quat%num_quat-1) ! cannot exceed num_cv
         qt_send_count(i)=max(0,min(j,quat%num_quat-j*(i-1))) ! how many CV I will send to CPU i
! atom map partitioning (for parallel computation of M
!
         imap_displ(i)=min((i-1)*k,cv%amap%last-1)
         imap_count(i)=max(0,min(k,cv%amap%last-k*(i-1)))
        enddo
       endif ! SIZE
! rmsd from a target structure
      else if (( keyword(1:4).eq.'RMSD'(1:4) )) then
       call smcv_rmsd_add(comlyn, comlen, rmsd)
! difference in the rmsd from two target structure (same routine!)
      else if (( keyword(1:5).eq.'DRMSD'(1:5) )) then
       call smcv_rmsd_add(comlyn, comlen, drmsd)
! normalized projection onto the vector connecting two structures aligned with the simulation structure
      else if (( keyword(1:4).eq.'PROJ'(1:4) )) then
       call smcv_rmsd_add(comlyn, comlen, proj)
! rms sum of CVs
      else if (( keyword(1:5).eq.'CVRMS'(1:5) )) then
       call smcv_cvrms_add(comlyn, comlen) ! rtmd-like variable (compatibility limited to steered dynamics as of 7.2010):
! z=sqrt( 1/N sum^N_i (z_i - z^0_i)^2 )
      else
        write(msg___,*)'UNRECOGNIZED CV TYPE: ',keyword;write(0,*) 'WARNING FROM: ',whoami,': ',msg___
      endif
!
      end subroutine smcv_add
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccc CV-SPECIFIC CODE ccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_posi_com_add(comlyn, comlen, cvtype)
      use cv_posi_com
      use ivector
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
!
      implicit none
 character(len=80) :: msg___
!
      character(len=*) comlyn
      integer :: comlen
      integer :: cvtype
! locals
      character(len=19) :: whoami
      character(len=11) :: cv_name
!
 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__
!
      integer :: i, j, nslct
      logical :: ok
      real*8 :: k, gam, w
      integer :: frame
      type (int_vector) :: posi_com_list ! for storing atom indices
!
      data whoami/' SMCV_POSI_COM_ADD>'/
      data cv_name/'POSI_COM'/
!
      select case (cvtype)
       case (posi_com_x); cv_name=' POSI_COM_X'
       case (posi_com_y); cv_name=' POSI_COM_Y'
       case (posi_com_z); cv_name=' POSI_COM_Z'
       case default;
        write(0,*) 'WARNING FROM: ',whoami,': ','UNKNOWN CV TYPE. NOTHING DONE.'
        return
      end select
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction coefficient
      frame=atoi(get_remove_parameter(comlyn, 'FRAM', comlen), 0) ! coordinate frame index for this position variable
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! now process atom selections;
! expecting 1 atom group
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      if (associated(iselct)) then ; nslct=size(iselct) ; else ; nslct=0 ; endif
!
      if(nslct.eq.0) then
       write(0,*) 'WARNING FROM: ',whoami,': ','ZERO ATOMS SELECTED'
      endif
!
      call int_vector_init(posi_com_list)
      do i=1,size(iselct)
       j=int_vector_add(posi_com_list,iselct(i))
       if (j.le.0) then ; write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD ATOM INDEX TO COM LIST' ; endif
      enddo
      if (associated(iselct)) deallocate(iselct)
!
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
       write(0,'(A)') msg___
!
       if (frame.ge.1) then
        write(msg___,'(A,I3)') whoami//' RELATIVE TO LOCAL FRAME ',frame
       else
        write(msg___,'(A)') whoami//' RELATIVE TO THE ABSOLUTE FRAME'
       endif
       write(0,'(A)') msg___
!
      endif
!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' AND GAMMA =',F7.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' WEIGHT =',F7.3,' AND GAMMA =',F7.3)
!
! now attempt to add CV
      ok=cv_posi_com_add(cvtype,posi_com_list,k,gam,w,max(frame,0)) ! no mass weighting; disallow negative frame indices
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD CV'
      endif
! deallocate lists
      call int_vector_done(posi_com_list)
! done!
      end subroutine smcv_posi_com_add
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_dihe_com_add(comlyn, comlen)
      use cv_dihe_com
      use ivector
!
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
!
      implicit none
!
      character(len=*) :: comlyn
      integer :: comlen
! locals
      character(len=19) :: whoami
      character(len=8) :: cv_name
!
 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__
 character(len=80) :: msg___
!
      integer :: i, j, atom_group, nslct
      logical :: ok
      real*8 :: k, gam, w
      type (int_vector), dimension(4) :: dihe_com_list ! for storing atom indices
!
      data whoami/' SMCV_DIHE_COM_ADD>'/
      data cv_name/'DIHE_COM'/
!
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction coeff.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! now process atom selections;
! expecting 4 atom selections that specify each atom group in succession;
! process each selection sequentially
      do atom_group=1,4
!
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      if (associated(iselct)) then ; nslct=size(iselct) ; else ; nslct=0 ; endif
!
       IF(NSLCT.EQ.0) THEN
        CALL WRNDIE(0,whoami,'ZERO ATOMS SELECTED')
        RETURN
       endif
!
       call int_vector_init(dihe_com_list(atom_group))



       do i=1,size(iselct)
        j=int_vector_add(dihe_com_list(atom_group),iselct(i))
        if (j.le.0) then ; write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD ATOM INDEX TO COM LIST' ; endif
       enddo
       if (associated(iselct)) deallocate(iselct)

      enddo ! loop over dihe_com selections
!
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
       write(0,'(A)') msg___
      endif
!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' AND GAMMA =',F7.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' WEIGHT =',F7.3,' AND GAMMA =',F7.3)
!
! now attempt to add CV
      ok=cv_dihe_com_add(dihe_com_list,k,gam,w) ! no mass weighting
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD CV'
      endif
! deallocate lists
      do i=1,4; call int_vector_done(dihe_com_list(i)); enddo
! done!
      end subroutine smcv_dihe_com_add
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_angle_com_add(comlyn, comlen)
      use cv_angle_com
      use ivector
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
      implicit none
!
 character(len=80) :: msg___
!
      character(len=*) :: comlyn
      integer :: comlen
! locals
      character(len=20) :: whoami
      character(len=9) :: cv_name
!



 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__

!
      integer :: i, j, atom_group, nslct
      logical :: ok
      real*8 :: k, gam, w
      type (int_vector), dimension(3) :: angle_com_list ! for storing atom indices
!
      data whoami/' SMCV_ANGLE_COM_ADD>'/
      data cv_name/'ANGLE_COM'/
!
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! now process atom selections;
! expecting 3 atom selections that specify each atom group in succession;
! process each selection sequentially
      do atom_group=1,3
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      if (associated(iselct)) then ; nslct=size(iselct) ; else ; nslct=0 ; endif
!
       if (nslct.eq.0) then
        write(0,*) 'WARNING FROM: ',whoami,': ','ZERO ATOMS SELECTED'
       endif
!
       call int_vector_init(angle_com_list(atom_group))





       do i=1,size(iselct)
        j=int_vector_add(angle_com_list(atom_group),iselct(i))
        if (j.le.0) then ; write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD ATOM INDEX TO COM LIST' ; endif
       enddo
       if (associated(iselct)) deallocate(iselct)

      enddo ! loop over angle_com selections
!
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
      endif
      write(0,'(A)') msg___
!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' AND GAMMA =',F7.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' WEIGHT =',F7.3,' AND GAMMA =',F7.3)
!
! now attempt to add CV
      ok=cv_angle_com_add(angle_com_list,k,gam,w) ! no mass weighting
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD CV'
      endif
! deallocate lists
      do i=1,3; call int_vector_done(angle_com_list(i)); enddo
! done!
      end subroutine smcv_angle_com_add
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_dist_com_add(comlyn, comlen)
      use cv_dist_com
      use ivector
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
!
      implicit none
!
      character(len=*) comlyn
      integer :: comlen
! locals
      character(len=19) :: whoami
      character(len=8) :: cv_name
!
 character(len=80) :: msg___



 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__

!
      integer :: i, j, atom_group, nslct
      logical :: ok
      real*8 :: k, gam, w
      type (int_vector), dimension(2) :: dist_com_list ! for storing atom indices
!
      data whoami/' SMCV_DIST_COM_ADD>'/
      data cv_name/'DIST_COM'/
!
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction coeff.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! now process atom selections;
! expecting 2 atom selections that specify each atom group in succession;
! process each selection sequentially
      do atom_group=1,2
!
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      if (associated(iselct)) then ; nslct=size(iselct) ; else ; nslct=0 ; endif
!
       if (nslct.eq.0) then
        write(0,*) 'WARNING FROM: ',whoami,': ','ZERO ATOMS SELECTED'
       endif
!
       call int_vector_init(dist_com_list(atom_group))







       do i=1,size(iselct)
        j=int_vector_add(dist_com_list(atom_group),iselct(i))
        if (j.le.0) then ; write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD ATOM INDEX TO COM LIST' ; endif
       enddo
       if (associated(iselct)) deallocate(iselct)

      enddo ! loop over dist_com selections
!
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
       write(0,'(A)') msg___
      endif
!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' AND GAMMA =',F7.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' WEIGHT =',F7.3,' AND GAMMA =',F7.3)
!
! now attempt to add CV
      ok=cv_dist_com_add(dist_com_list,k,gam,w) ! no mass weighting
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD CV'
      endif
! deallocate lists
      do i=1,2; call int_vector_done(dist_com_list(i)); enddo
! done!
      end subroutine smcv_dist_com_add
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_anglvec_add(comlyn, comlen)
      use cv_anglvec
      use ivector
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
!
      implicit none
!
      character(len=*) :: comlyn
      integer :: comlen
! locals
      character(len=18) :: whoami
      character(len=8) :: cv_name
      character(len=8) :: key
!
 character(len=80) :: msg___
!



 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__

!
      integer :: i, ipt, j, l, nslct
      logical :: ok
      real*8 :: k, gam, w
      type (int_vector), dimension(4) :: atom_list ! for storing atom indices
      integer :: f1, f2
      real*8 :: p(4,3) ! for point definition
      logical :: qp1=.false., qp2=.false., qp3=.false., qp4=.false.
!
      data whoami/' SMCV_ANGLVEC_ADD>'/
      data cv_name/'ANGLVEC'/
!
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! check for frame specification, so that comlyn is cleaned up
      f1=atoi(get_remove_parameter(comlyn, 'FR1', comlen), 0)
      f2=atoi(get_remove_parameter(comlyn, 'FR2', comlen), 0)
!
! initialize atom arrays
      do i=1,4; call int_vector_init(atom_list(i)); enddo
!
      p=0d0 ! initialize points
      do l=1,4
! now process vector specifications (expecting four points)
       key=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
       if (( key(1:2).eq.'P1'(1:2) )) then
        ipt=1 ; qp1=.true.
       elseif (( key(1:2).eq.'P2'(1:2) )) then
        ipt=2 ; qp2=.true.
       elseif (( key(1:2).eq.'P3'(1:2) )) then
        ipt=3 ; qp3=.true.
       elseif (( key(1:2).eq.'P4'(1:2) )) then
        ipt=4 ; qp4=.true.
       else
        write(msg___,*)'VECTOR DEFINITION ERROR. EXPECTED "P[1-4]", GOT "',key,'"';write(0,*) 'WARNING FROM: ',whoami,': ',msg___
       endif
!
       comlen=min(max(0,comlen),len(comlyn));comlyn(comlen+1:)='';call adjustleft(comlyn,(/' ',tab/));comlen=len_trim(comlyn)
       if (find_tag(comlyn, 'SELE', comlen).eq.1) then ! next word is select
!ccccccccccccc process atom selection
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      if (associated(iselct)) then ; nslct=size(iselct) ; else ; nslct=0 ; endif
!
        IF(NSLCT.EQ.0) THEN
         write(0,*) 'WARNING FROM: ',whoami,': ','ZERO ATOMS SELECTED'
        endif
!







       do i=1,size(iselct)
        j=int_vector_add(atom_list(ipt),iselct(i))
        if (j.le.0) then ; write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD ATOM INDEX TO COM LIST' ; endif
       enddo
       if (associated(iselct)) deallocate(iselct)

!
       else ! specify point manually
        do i=1,3; p(ipt,i)=atof(pop_string(comlyn,comlen)) ; comlen=len_trim(comlyn); enddo
       endif
      enddo
! check that all four points have been added
      if (.not.(qp1.and.qp2.and.qp3.and.qp4)) then
       write(0,*) 'WARNING FROM: ',whoami,': ','SOME POINTS WERE NOT DEFINED. NOTHING DONE'
       return
      endif
!
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
       write(0,'(A)') msg___
      endif
!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' AND GAMMA =',F7.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' WEIGHT =',F7.3,' AND GAMMA =',F7.3)
!
! now attempt to add CV
      ok=cv_anglvec_add(atom_list,p,f1,f2,k,gam,w) ! no mass weighting
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD CV'
      endif
! deallocate lists
      do i=1,4; call int_vector_done(atom_list(i)); enddo
! done!
      end subroutine smcv_anglvec_add
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 'frame' is not really a CV but processed here for convenience
      subroutine smcv_frame_add(comlyn, comlen)
      use cv_frames
      use ivector
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
!
      implicit none
 character(len=80) :: msg___
!
      character(len=*) :: comlyn
      integer :: comlen
! locals
      character(len=16) :: whoami



 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__

      integer :: i, j, nslct
      logical :: ok
      type (int_vector) :: frame_list ! for storing atom indices
!
      data whoami/' SMCV_FRAME_ADD>'/
!
! process atom selections;
! specify atom group
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      if (associated(iselct)) then ; nslct=size(iselct) ; else ; nslct=0 ; endif
!
      if (nslct.lt.4) then ! require at least four atoms, otherwise, can never define frame uniquely
       write(0,*) 'WARNING FROM: ',whoami,': ',' FEWER THAN FOUR ATOMS SELECTED.'
      endif
!
      call int_vector_init(frame_list)
      do i=1,size(iselct)
       j=int_vector_add(frame_list,iselct(i))
       if (j.le.0) then
        write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD ATOM INDEX TO FRAME LIST'
       endif
      enddo
      if (associated(iselct)) deallocate(iselct)
!
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       write(msg___, 665) whoami ; write(0,'(A)') msg___
      endif
!
  665 format(/A,' WILL ADD NEW REFERENCE FRAME')
! now attempt to add frame
      ok=(frames_add(frame_list).gt.0) !
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD FRAME'
      endif
! deallocate lists
      call int_vector_done(frame_list)
! done!
      end subroutine smcv_frame_add
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_quaternion_add(comlyn, comlen)
      use cv_qcomp
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
      implicit none
!
 character(len=80) :: msg___
!
      character(len=*) :: comlyn
      integer :: comlen
! locals
      character(len=21) :: whoami
      character(len=10) :: cv_name
      logical :: ok
      real*8 :: k, gam, w
      integer :: f1, f2
!
      data whoami/' SMCV_QUATERNION_ADD>'/
      data cv_name/'QUATERNION'/
!
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction
      f1=atoi(get_remove_parameter(comlyn, 'FRA1', comlen), 0) ! coordinate frame index 1
      f2=atoi(get_remove_parameter(comlyn, 'FRA2', comlen), 0) ! coordinate frame index 2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
!
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
       write(0,'(A)') msg___
!
       if (f1.ge.1) then
        write(msg___,'(A,I3)') whoami//' FRAME1: LOCAL FRAME #',f1
       else
        write(msg___,'(A)') whoami//' FRAME1: ABSOLUTE FRAME'
       endif
       write(0,'(A)') msg___
!
       if (f2.ge.1) then
        write(msg___,'(A,I3)') whoami//' FRAME2: LOCAL FRAME #',f2
       else
        write(msg___,'(A)') whoami//' FRAME2: ABSOLUTE FRAME'
       endif
       write(0,'(A)') msg___
!
      endif
!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' AND GAMMA =',F7.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F7.3, &
     & ' WEIGHT =',F7.3,' AND GAMMA =',F7.3)
!
! now attempt to add quaternion components, one by one:
      ok= cv_qcomp_add(qcomp_1, f1, f2, k, gam, w) &
     & .and.cv_qcomp_add(qcomp_2, f1, f2, k, gam, w) &
     & .and.cv_qcomp_add(qcomp_3, f1, f2, k, gam, w) &
     & .and.cv_qcomp_add(qcomp_4, f1, f2, k, gam, w)
!
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD QUATERNION CV'
      endif
! done!
      end subroutine smcv_quaternion_add
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_cvrms_add(comlyn, comlen)
! this CV is an RMS combination of existing CV; experimental and not fully implemented (intented for steered dynamics)
      use cv_cvrms
      use parselist
      use ivector
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use parser
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      implicit none
      character(len=80) :: msg___
!
      character(len=*) :: comlyn
      integer :: comlen
! locals
      character(len=16) :: whoami
      character(len=5) :: cv_name
      logical :: ok
      real*8 :: k, gam, w
      type (int_vector) :: cv_list ! for storing cv indices used to calculate RMS
!
      data whoami/' SMCV_CVRMS_ADD>'/
      data cv_name/'CVRMS'/
!
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction coeff.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! specify cv to use in the RMS definition
      call parse_list(cv_list, comlyn) ! will return allocated cv_list with the indices
!
! print short summary
      if (MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0) then
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
      endif
      write(0,'(A)') msg___
!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F8.3, &
     & ' AND GAMMA =',F8.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F8.3, &
     & ' WEIGHT =',F8.3,' AND GAMMA =',F8.3)
!
! now attempt to add CV
      ok=cv_cvrms_add(cv_list,k,gam,w) ! no mass weighting
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD CV'
      endif
! deallocate list
      call int_vector_done(cv_list)
! done!
      end subroutine smcv_cvrms_add
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smcv_rmsd_add(comlyn, comlen, cvtype)
! much of the code taken from RTMD sources in ../misc
      use cv_rmsd
      use cv_drmsd
      use cv_proj, only : cv_proj_add
      use cv_types
      use output,only:message,warning,plainmessage,output_init,output_done,fatal_warning,fout
      use psf
      use system, only : r, rcomp, m, bfactor, occupancy
      use parser
      use constants
      use multicom_aux !!**CHARMM_ONLY**!##MULTICOM
      use mpi
      use system, only : system_getind
!
      implicit none
      character(len=*) :: comlyn
      integer :: comlen
      integer :: cvtype
! locals
      character(len=15) :: whoami
      character(len=5) :: cv_name
      logical :: ok, oset
      real*8 :: k, gam, w, rtarget_com(3)
!
      real*8, pointer :: rtarget_o(:,:), rtarget_f(:,:), &
     & rtarget1_o(:,:), rtarget1_f(:,:), &
     & orientWeights(:), forcedWeights(:)
      integer, pointer :: iatom_o(:), iatom_f(:)
      integer :: norient, nforced
!
 integer :: isele, i__, iend; integer, pointer::iselct(:);character(LEN=20)::word__
 integer, pointer :: jselct(:), kselct(:)
 character(len=80) :: msg___
!
      integer :: i, j, n
      real*8 :: a, b
!
      logical :: use_main, use_comp, qroot, qmass, qtwo ! qtwo: true if using two target structures
!
      data whoami/' SMCV_RMSD_ADD>'/
!
      select case (cvtype)
       case(rmsd ); cv_name='RMSD '; qtwo=.false.
       case(drmsd); cv_name='DRMSD'; qtwo=.true.
       case(proj ); cv_name='PROJ '; qtwo=.true.
       case default
        write(0,*) 'WARNING FROM: ',whoami,': ','UNKNOWN CV REQUESTED. NOTHING DONE.';
        return
      end select
!
      qroot=(MPI_COMM_STRNG.ne.MPI_COMM_NULL.and.ME_STRNG.eq.0)
!
! first check for CV options
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      k=atof(get_remove_parameter(comlyn, 'FORC', comlen), 0.0d0) ! can specify force constant manually
      w=atof(get_remove_parameter(comlyn, 'WEIG', comlen), -1.0d0) ! can specify weight manually
      gam=atof(get_remove_parameter(comlyn, 'GAMM', comlen), 1.0d0) ! friction coeff.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! process atom selections
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nullify(iselct, jselct)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! look for orientation atom selection
      oset=.false. ! is there a separate orientation specification?
      i=remove_tag(comlyn,'ORIE',comlen) ! find the position of `ORIE'
      oset=(i.gt.0)
! first determine whether a selection keyword follows orie
      if (oset) then
       n=comlen-i+1
       j=find_tag(comlyn(i:comlen), 'SELE', n)
       if (j.le.0) then ! only if the ORIE directive exists
         write(0,*) 'WARNING FROM: ',whoami,': ','ATOM SELECTION MUST BE SPECIFIED AFTER ORIE.'
         return
        endif
      endif
!*************************************************************************************
! look for the first selection
!*************************************************************************************
      j=find_tag(comlyn, 'SELE', comlen)
      if (j.eq.0) then ! sele keyword does not exist
       write(0,*) 'WARNING FROM: ',whoami,': ','RMSD ATOMS NOT SPECIFIED. NOTHING DONE.';
       return
      elseif (j.le.i.or..not.oset) then ! no 'orie', or first occurrence of sele before orie (deleted by __INDXa above)
! assume that this is the forcing set selection
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      jselct=>iselct ; nullify(iselct)
! orientation selection
       if (oset) then
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
       endif
!
      else ! first selection after orie
! first selection (orientation atoms)
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      kselct=>iselct; nullify(iselct) ! have to do this because macro uses iselct
! second selection (forcing atoms)
          nullify(iselct)
          isele=find_tag(comlyn, 'SELE', comlen)
          msg___=comlyn(isele:comlen) ! part of string that begins with the selection
          i__=comlen-isele+1
          iend=find_tag(msg___, 'END', i__) ! location of selection termination
          iend=iend-1+isele ! index into comlyn starting from 1
          msg___=comlyn(isele:iend) ! part of string that begins with the selection and ends before ' END'
          word__=pop_string(msg___) ! remove first word (which we know is 'SELE*') from msg___
! process selection:
          nullify(iselct)
          iselct=>system_getind(msg___)
! remove selection string from command line:
          msg___=comlyn(iend:comlen) ! command line starting with 'END' (see above)
          word__=pop_string(msg___) ! remove 'END*' e.g. 'ENDING' is ok too
          comlyn(isele:isele)=' ';
          comlyn(isele+1:)=msg___ ! selection has been removed from command
          comlen=len_trim(comlyn)
      jselct=>iselct
      iselct=>kselct ; nullify(kselct)
!
      if (associated(iselct)) then ; norient=size(iselct) ; else ; norient=0 ; endif
      if (associated(jselct)) then ; nforced=size(jselct) ; else ; nforced=0 ; endif

!
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (nforced.eq.0) then
       write(0,*) 'WARNING FROM: ',whoami,': ','RMSD ATOMS NOT SPECIFIED. NOTHING DONE.'
       return
      endif
!
      if (.not.oset) then ! use forced atoms for orientation too
       norient=nforced



       iselct=>jselct

      endif
!
! currently we require at least three atoms for orientation
!
      if (norient.lt.3) then
        write(0,*) 'WARNING FROM: ',whoami,': ',' FEWER THAN THREE ATOMS SELECTED FOR ORIENTATION. ABORT.'
        return
      endif
!
      qmass=(remove_tag(comlyn,'MASS',comlen).gt.0) ! mass-weighting flag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! done with selections
      allocate(iatom_o(norient), &
     & iatom_f(nforced), &
     & rtarget_o(norient,3), &
     & rtarget_f(nforced,3), &
     & orientWeights(norient), &
     & forcedWeights(nforced))
      if (qtwo) allocate(rtarget1_o(norient,3),rtarget1_f(nforced,3))
!
! initialize arrays
      iatom_o=0; iatom_f=0;
      rtarget_o=0.0d0; rtarget_f=0.0d0;
      if (qtwo) then ; rtarget1_o=0.0d0; rtarget1_f=0.0d0; endif ! second structure
      orientWeights=1d0; forcedWeights=1d0
!
! build index arrays
! NOTE that in DMOL, iselct has a different meaning: it stores the atom indices, rather than flags
      if (associated(iselct)) then ; iatom_o=iselct ; deallocate(iselct) ; endif
      if (associated(jselct)) then ; iatom_f=jselct ; deallocate(jselct) ; endif
!!!!!!!!! mass-weighting
      if (qmass) then
        do i=1,norient
         orientWeights(i)= &
     & m(iatom_o(i))*orientWeights(i)
        enddo
!
        do i=1, nforced
         forcedWeights(i)= &
     & m(iatom_f(i))*forcedWeights(i)
        enddo
!
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load and save coordinates of the target structure
      use_main=((remove_tag(comlyn,'MAIN',comlen)).gt.0)
      use_comp=((remove_tag(comlyn,'COMP',comlen)).gt.0)
      if (use_comp) then
        if (use_main) then
         write(0,*) 'WARNING FROM: ',whoami,': ','MAIN AND COMP CANNOT BOTH BE SPECIFIED. USING MAIN.'
         use_comp=.false.
        endif
      else
        use_main=.true. ! default
      endif
!
      if (use_comp) then
! orient
        do i=1,norient
         rtarget_o(i,1)=rcomp(1,iatom_o(i))
         rtarget_o(i,2)=rcomp(2,iatom_o(i))
         rtarget_o(i,3)=rcomp(3,iatom_o(i))
        enddo
! forced
        do i=1,nforced
         rtarget_f(i,1)=rcomp(1,iatom_f(i))
         rtarget_f(i,2)=rcomp(2,iatom_f(i))
         rtarget_f(i,3)=rcomp(3,iatom_f(i))
        enddo
! second reference structure:
        if (qtwo) then
         do i=1,norient
          rtarget1_o(i,1)=r(1,iatom_o(i))
          rtarget1_o(i,2)=r(2,iatom_o(i))
          rtarget1_o(i,3)=r(3,iatom_o(i))
         enddo
!
         do i=1,nforced
          rtarget1_f(i,1)=r(1,iatom_f(i))
          rtarget1_f(i,2)=r(2,iatom_f(i))
          rtarget1_f(i,3)=r(3,iatom_f(i))
         enddo
        endif ! qtwo
!
      else ! use main coordinates
! orient
        do i=1,norient
         rtarget_o(i,1)=r(1,iatom_o(i))
         rtarget_o(i,2)=r(2,iatom_o(i))
         rtarget_o(i,3)=r(3,iatom_o(i))
        enddo
! forced
        do i=1,nforced
         rtarget_f(i,1)=r(1,iatom_f(i))
         rtarget_f(i,2)=r(2,iatom_f(i))
         rtarget_f(i,3)=r(3,iatom_f(i))
        enddo
! second reference structure:
        if (qtwo) then
         do i=1,norient
          rtarget1_o(i,1)=rcomp(1,iatom_o(i))
          rtarget1_o(i,2)=rcomp(2,iatom_o(i))
          rtarget1_o(i,3)=rcomp(3,iatom_o(i))
         enddo
!
         do i=1,nforced
          rtarget1_f(i,1)=rcomp(1,iatom_f(i))
          rtarget1_f(i,2)=rcomp(2,iatom_f(i))
          rtarget1_f(i,3)=rcomp(3,iatom_f(i))
         enddo
        endif ! qtwo
!
      endif ! use_comp
!
! check for undefined values
      if (any(rtarget_o.eq.unknownf).or.any(rtarget_f.eq.unknownf)) then
        write(0,*) 'WARNING FROM: ',whoami,': ','FIRST TARGET STRUCTURE HAS UNDEFINED COORDINATES. QUITTING.'
        return
      elseif (qtwo) then
       if (any(rtarget1_o.eq.unknownf).or.any(rtarget1_f.eq.unknownf)) then
        write(0,*) 'WARNING FROM: ',whoami,': ','SECOND TARGET STRUCTURE HAS UNDEFINED COORDINATES. QUITTING.'
        return
       endif
      endif
!
! normalize weighting coefficients
! note: routines in rtmd_aux do _not_ perform normalization; result affects FP precision, analytically, there is no difference
!
      a=sum(orientWeights)
      b=sum(forcedWeights)
      if (abs(a).gt.errtol()) then
        a=1d0/a
        orientWeights=a*orientWeights
      endif
!
      if (abs(b).gt.errtol()) then
        b=1d0/b
        forcedWeights=b*forcedWeights
      endif
!
! translate target structure so that its centroid is at the origin
      rtarget_com=matmul(transpose(rtarget_o),orientWeights)
!
      rtarget_o(:,1)=rtarget_o(:,1)-rtarget_com(1)
      rtarget_o(:,2)=rtarget_o(:,2)-rtarget_com(2)
      rtarget_o(:,3)=rtarget_o(:,3)-rtarget_com(3)
!
      rtarget_f(:,1)=rtarget_f(:,1)-rtarget_com(1)
      rtarget_f(:,2)=rtarget_f(:,2)-rtarget_com(2)
      rtarget_f(:,3)=rtarget_f(:,3)-rtarget_com(3)
!
      if (qtwo) then ! repeat for second structure
       rtarget_com=matmul(transpose(rtarget1_o),orientWeights)
!
       rtarget1_o(:,1)=rtarget1_o(:,1)-rtarget_com(1)
       rtarget1_o(:,2)=rtarget1_o(:,2)-rtarget_com(2)
       rtarget1_o(:,3)=rtarget1_o(:,3)-rtarget_com(3)
!
       rtarget1_f(:,1)=rtarget1_f(:,1)-rtarget_com(1)
       rtarget1_f(:,2)=rtarget1_f(:,2)-rtarget_com(2)
       rtarget1_f(:,3)=rtarget1_f(:,3)-rtarget_com(3)
      endif
!
! print summary
      if (qroot) then
!
       if (w.lt.0d0) then
        write(msg___, 664) whoami,cv_name,k,gam
       else
        write(msg___, 665) whoami,cv_name,k,w,gam
       endif
       write(0,'(A)') msg___
!
       write(msg___,103) whoami, nforced ; write(0,'(A)') msg___
       write(msg___,100) whoami, norient ; write(0,'(A)') msg___
       if (qmass) then ; write(msg___,102) whoami ; write(0,'(A)') msg___ ; endif
       if (qtwo) then
        if (use_comp) then
         write(msg___,105) whoami, 'COMP'
        else
         write(msg___,105) whoami, 'MAIN'
        endif
       else
        if (use_comp) then
         write(msg___,104) whoami, 'COMP'
        else
         write(msg___,104) whoami, 'MAIN'
        endif
       endif
       write(0,'(A)') msg___!
  664 format(/A,' WILL ADD ',A,' CV WITH K =',F8.3, &
     & ' AND GAMMA =',F8.3,'.')
  665 format(/A,' WILL ADD ',A,' CV WITH K =',F8.3, &
     & ' WEIGHT =',F8.3,' AND GAMMA =',F8.3)
!
 100 format(A,' WILL ORIENT TARGET STRUCTURE(S) BASED ON ',i5, &
     & ' ATOMS')
 102 format(A,' WILL USE MASS-WEIGHTING.')
 103 format(A,' ',i5,' FORCED ATOMS FOUND.')
 104 format(A,' TARGET STRUCTURE TAKEN FROM ',A,' SET.')
 105 format(A,' FIRST TARGET STRUCTURE TAKEN FROM ',A,' SET.')
!
      endif ! qroot
! now attempt to add CV
      ok=.false.
      select case(cvtype)
       case(rmsd);
        ok=cv_rmsd_add(iatom_o, iatom_f, rtarget_o, rtarget_f, &
     & orientWeights, forcedWeights, &
     & k, gam, w)
       case(drmsd);
        ok=cv_drmsd_add(iatom_o, iatom_f, rtarget_o, rtarget_f, &
     & rtarget1_o, rtarget1_f, &
     & orientWeights, forcedWeights, &
     & k, gam, w)
       case(proj);
        ok=cv_proj_add(iatom_o, iatom_f, rtarget_o, rtarget_f, &
     & rtarget1_o, rtarget1_f, &
     & orientWeights, forcedWeights, &
     & k, gam, w)
      end select
!
      if (.not.ok) then
       write(0,*) 'WARNING FROM: ',whoami,': ','COULD NOT ADD CV'
      endif
! deallocate variables
      deallocate(iatom_o, iatom_f, rtarget_o, rtarget_f, &
     & orientWeights, forcedWeights)
      if (qtwo) deallocate(rtarget1_o, rtarget1_f)
!
      end subroutine smcv_rmsd_add
!
!**CHARMM_ONLY**!##ENDIF
