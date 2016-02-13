      interface
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
