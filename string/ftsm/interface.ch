! ********************************************************************************!
! Preprocessed by GNU CPP !
! ********************************************************************************!
      interface
!
       function sm_get_column(cmd_, l, qcoltag, missing, coltag_) result(C)
        implicit none
        character(len=*) :: cmd_
        integer :: l, missing
        logical :: qcoltag
        integer :: C
        character(len=*), optional :: coltag_
       end function sm_get_column
!
      end interface
