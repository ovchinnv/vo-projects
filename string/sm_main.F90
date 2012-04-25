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
!
!**CHARMM_ONLY**!##IF STRINGM
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sm_main(COMLYN,COMLEN)
      use sm0k, only: sm0k_main
      use ftsm, only: ftsm_parse
      implicit none
!-----------------------------------------------
! calls string method parsers
!-----------------------------------------------
      character(len=*) :: comlyn
      integer :: comlen
!
! local variables
!
      character(len=8) :: keyword
      character(len=9) :: whoami
 character(len=80) :: msg__
!
      data whoami /' SM_MAIN>'/
!
      keyword=pop_string(comlyn,comlen) ; comlen=len_trim(comlyn)
! there are two parsers, 0K , CV , and FTSM
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (( keyword(1:4).eq.'ZERO'(1:4) )) then
        call sm0k_main(comlyn, comlen)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'COLV'(1:4) )) then
        call smcv(comlyn, comlen) ! SMCV parser
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'FTSM'(1:4) )) then
        call ftsm_parse(comlyn, comlen) ! FTSM parser
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'OPEN'(1:4) )) then
        call sm_open(comlyn, comlen, .true.) ! open a file on each replica
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      elseif (( keyword(1:4).eq.'CLOS'(1:4) )) then
        call sm_close(comlyn, comlen) ! close a file on each replica
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
            write(msg__,*)'UNRECOGNIZED SUBCOMMAND: ',keyword;write(0,*) 'WARNING FROM: ',whoami,': ',msg__
      endif
!
      end subroutine sm_main
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       SUBROUTINE SM_OPEN(COMLYN,COMLEN, loud)
!-----------------------------------------------------------------------
! Opens 1 file per replica process
!-----------------------------------------------------------------------
! much of this code is duplicated from the main charmm open routine (vopen)
       use smcv_var, only: smcv_initialized
       use sm0k, only: sm0k_initialized
       use ftsm_var, only: ftsm_initialized
!
       use files; use parser
       use mpi
       use multicom_aux
!
       ! implicit none
!
       character(len=*) :: COMLYN
       integer :: COMLEN
       character(len=132) :: filex, formt, facc
       integer :: unum, flen
!
       logical :: loud
       character(len=(9) :: whoami
       data whoami /' SM_OPEN>'/
!
       if (.not.(smcv_initialized &
     & .or.sm0k_initialized &
     & .or.ftsm_initialized)) then
        if (ME_GLOBAL.eq.0) &
        write(0,*) 'WARNING FROM: ',whoami,': ',' STRING METHOD NOT INITIALIZED. NOTHING DONE.'
        return
       endif
!
       if (ME_STRNG.ne.MPI_UNDEFINED) then ! only roots work
! unum
        UNUM=atoi(get_remove_parameter(COMLYN, 'UNIT', COMLEN), -1)
        IF (UNUM.LT.0) THEN
         write(0,*) 'WARNING FROM: ',whoami,': ',' NO UNIT NUMBER SPECIFIED'
        ENDIF
! filename
        FILEX=get_remove_parameter(COMLYN, COMLEN, 'NAME'); FLEN=len_trim(FILEX)
        IF (FLEN.LE.0) THEN
         write(0,*) 'WARNING FROM: ',whoami,': ','NO FILE NAME GIVEN'
         return
        ENDIF
! format
        IF (remove_tag(COMLYN,'UNFO',COMLEN).GT.0) THEN
         FORMT='UNFORMATTED'
        ELSE IF (remove_tag(COMLYN,'FILE',COMLEN).GT.0) THEN
         FORMT='UNFORMATTED'
        ELSE IF (remove_tag(COMLYN,'FORM',COMLEN).GT.0) THEN
         FORMT='FORMATTED'
        ELSE IF (remove_tag(COMLYN,'CARD',COMLEN).GT.0) THEN
         FORMT='FORMATTED'
        ELSE
 write(0,*) 'WARNING FROM: ',whoami,': ','NO FORMAT SPECIFIED, WILL USE "UNFORMATTED"'
         FORMT='UNFORMATTED'
        ENDIF
!
! access
        IF (remove_tag(COMLYN,'APPE',COMLEN).GT.0) THEN
         FACC='APPEND'
        ELSE IF (remove_tag(COMLYN,'READ',COMLEN).GT.0) THEN
         FACC='READ'
        ELSE IF (remove_tag(COMLYN,'WRIT',COMLEN).GT.0) THEN
         FACC='WRITE'
        ELSE
          FACC='READ'
        ENDIF
       call files_open(unum,filex,formt,facc)
       comlen=min(max(0,comlen),len(comlyn));comlyn(comlen+1:)='';call adjustleft(comlyn,(/' ',tab/));comlen=len_trim(comlyn)
       call mpi_bcast(comlyn,len(comlyn),MPI_CHARACTER,0,MPI_COMM_STRNG,ierror)
!
       END SUBROUTINE SM_OPEN
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine sm_close(comlyn, comlen)
!
       use sm_var, only: smcv_initialized
       use sm0k, only: sm0k_initialized
       use ftsm_var, only: ftsm_initialized
!
       use files; use parser
       use mpi
       use multicom_aux
!
       implicit none
!
       CHARACTER(len=*) :: COMLYN
       integer :: COMLEN
!
       integer :: UNUM
!
       character(len=10) :: whoami
       data whoami /' SM_CLOSE>'/
!
       if (.not.(smcv_initialized &
     & .or.sm0k_initialized &
     & .or.ftsm_initialized)) then
        if (ME_GLOBAL.eq.0) &
        write(0,*) 'WARNING FROM: ',whoami,': ',' STRING METHOD NOT INITIALIZED. NOTHING DONE.'
        return
       endif
!
       UNUM=atoi(get_remove_parameter(COMLYN, 'UNIT', COMLEN), -1)
!
       if(ME_STRNG.ne.MPI_UNDEFINED) then ! only ensemble heads work
!
        call files_close(unum)
       endif
!
       end subroutine sm_close
!**********************************************************************
!**CHARMM_ONLY**!##ENDIF
