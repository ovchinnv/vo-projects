! ********************************************************************************!
! Preprocessed by GNU CPP !
! ********************************************************************************!
      interface
       subroutine frames_align_string(x,y,z,mass,min_rmsd,ind)
        use stream
        implicit none
        real(chm_real) :: x(:), y(:), z(:), mass(:)
        logical, optional :: min_rmsd
        integer, optional :: ind ! frame index
       end subroutine frames_align_string
!
       subroutine frame_align_rmsd(x,y,z,mass,ind)
        use stream
        implicit none
        real(chm_real) :: x(:), y(:), z(:), mass(:)
        integer, optional :: ind ! frame index
       end subroutine frame_align_rmsd
!
       subroutine frame_align_voro(x,y,z,mass,ind)
        use stream
        implicit none
        real(chm_real) :: x(:), y(:), z(:), mass(:)
        integer, optional :: ind ! frame index
       end subroutine frame_align_voro
!
       subroutine smcv_init(maxcv)
        implicit none
        integer, optional :: maxcv
       end subroutine smcv_init
!
       function sm_get_column(cmd_, l, qcoltag, missing) result(C)
        implicit none
        character(len=*) :: cmd_
        integer :: l, missing
        logical :: qcoltag
        integer :: C
       end function sm_get_column
!
      end interface
