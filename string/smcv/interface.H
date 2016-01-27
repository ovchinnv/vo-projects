      interface
       subroutine frames_align_string(x,y,z,mass,min_rmsd,ind)
        __DEP_OUTPUT
        __IMPNONE
        float :: x(:), y(:), z(:), mass(:)
        bool, optional :: min_rmsd
        int, optional :: ind ! frame index
       end subroutine frames_align_string
!
       subroutine frame_align_rmsd(x,y,z,mass,ind)
        __DEP_OUTPUT
        __IMPNONE
        float :: x(:), y(:), z(:), mass(:)
        int, optional :: ind ! frame index
       end subroutine frame_align_rmsd
!
       subroutine frame_align_voro(x,y,z,mass,ind)
        __DEP_OUTPUT
        __IMPNONE
        float :: x(:), y(:), z(:), mass(:)
        int, optional :: ind ! frame index
       end subroutine frame_align_voro
!
       subroutine smcv_init(maxcv)
        __IMPNONE
        int, optional :: maxcv
       end subroutine smcv_init
!
       function sm_get_column(cmd_, l, qcoltag, missing, coltag_) result(C)
        __IMPNONE
        character(len=*) :: cmd_
        int :: l, missing
        bool :: qcoltag
        int :: C
        character(len=*), optional :: coltag_
       end function sm_get_column
!
      end interface
