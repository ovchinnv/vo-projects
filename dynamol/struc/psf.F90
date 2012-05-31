/*COORDINATES AND MASSES:*/
/*#define __INDX(__STR, __STRLEN, __TEST, __TESTLEN)  index(__STR(1:min(__STRLEN,len(__STR))),__TEST(1:min(__TESTLEN,len(__TEST))))*/
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
      module psf
      use parser, only: adjustleft, numword
      use tlist
      use psfatom
!
      private
! read, parse & store CHARMM parameter/topology file;
!
      character, parameter :: tab=char(9)
      character, parameter :: comment(1) = (/'*'/)
!
      logical, save :: psf_initialized=.false.
      type (toplist), public, save :: blist, alist, dlist, ilist, clist
      type (atomlist), public, save :: atoms
! missing : hb donor/acceptor lists, nb exclusion lists
!
      public psf_read ! read input file and store all parameters
      private psf_init
      public psf_done
      public psf_info ! list parameters
      public psf_print
      public psf_natom ! return number of atoms
!
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine psf_init()
       implicit none
! initialize params structure
       call psf_done()
       call atomlist_init(atoms) ! atom list
       call toplist_init(blist, 3) ! bond list
       call toplist_init(alist, 4) ! angle list
       call toplist_init(dlist, 5) ! dihedral list
       call toplist_init(ilist, 5) ! improper list
! call toplist_init(clist, 9) ! CMAP list
       psf_initialized=.true.
       end subroutine psf_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine psf_done()
       implicit none
       call atomlist_done(atoms) ! atom list
       call toplist_done(blist) ! bond list
       call toplist_done(alist) ! angle list
       call toplist_done(dlist) ! dihedral list
       call toplist_done(ilist) ! improper list
       call toplist_done(clist) ! CMAP list
       psf_initialized=.false.
       end subroutine psf_done
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       subroutine psf_read(fid, &
     & xpl)
!
       use output
       implicit none
!
!
       integer :: fid ! input file handle
       logical, optional :: xpl
       logical :: xplor
!
       integer :: i, l, j, k
       character(len=8), parameter :: whoami = 'PSF_READ'
       character(len=200) :: keyword, cmdline
!
       integer, parameter :: qhead=1, qtitl=2, qatom=3, qbond=4, qangl=5, qdihe=6, qimpr=7, qdon=8,&
& qacc=9, qnbex=10, qgrp=11, qmolnt=12, qlp=13, qcmap=14
       logical, dimension(14), save :: qread=.false. ! this array
       integer :: which
       integer :: natom=-1, nbond=-1, nangl=-1, ndihe=-1, nimpr=-1, ndon=-1, nacc=-1,&
& nnbx=-1, ngrp=-1, nst2=-1, nmolnt=-1, nlp=-1, nlph=-1, ncmap=-1
       integer, allocatable :: idummy(:)
       integer :: ioerr
! for atomlist:
       character(len=8) :: segid, resid, resname, type, aname
       integer :: typeid
       real*8 :: charge, mass
!
       if (present(xpl)) then ; xplor=xpl; else ; xplor=.false.; endif
!
       which=qhead
       if (.not. psf_initialized) call psf_init()
! call message(whoami, 'Reading CHARMM Structure file (PSF).')
       do while (.true.)




        read(fid,'(A)',IOSTAT=ioerr) cmdline

        if (ioerr.eq.0) then



         call adjustleft(cmdline)
         l=len_trim(cmdline)
! ignore lines that begin with comment symbols
         if (any(comment.eq.cmdline(1:1))) cycle
! skip empty lines (they indicate the end of the current PSF section and the beginning of the next section)
         if (l.eq.0) then
          which=min(which, -which) ! set to negative so that we can't read until which is defined properly
          cycle
         endif
!
         if (l.ge.3) then
          keyword=cmdline(1:3)
          if (keyword.eq.'PSF') then
           if (.not.(any(qread))) then ! can only be the first section
            qread(qhead)=.true.
            which=qhead+1
            cycle
           else
            call error(whoami,'MISPLACED "PSF" TAG IN STRUCTURE FILE. ABORTING.',0)
            call error(whoami,cmdline(1:l),-1)
            return
           endif
          endif
         endif
! see if command string has a comment that identifies the next section:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         i=scan(cmdline(1:l), '!')
         if (i.gt.0.and.i+1.le.l) then
           keyword=cmdline(i+1:l)
           call adjustleft(keyword)
           if (keyword(1:6).eq.'NTITLE') then
            if (.not.qread(qtitl)) then
             qread(qtitl)=.true.
             which=qtitl+1
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:5).eq.'NATOM') then
            if (.not.qread(qatom)) then
             qread(qatom)=.true.
             read(cmdline,*) natom
             which=qatom
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:5).eq.'NBOND') then
            if (.not.qread(qbond)) then
             qread(qbond)=.true.
             read(cmdline,*) nbond
             which=qbond
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:6).eq.'NTHETA') then
            if (.not.qread(qangl)) then
             qread(qangl)=.true.
             read(cmdline,*) nangl
             which=qangl
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:4).eq.'NPHI') then
            if (.not.qread(qdihe)) then
             qread(qdihe)=.true.
             read(cmdline,*) ndihe
             which=qdihe
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:6).eq.'NIMPHI') then
            if (.not.qread(qimpr)) then
             qread(qimpr)=.true.
             read(cmdline,*) nimpr
             which=qimpr
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:4).eq.'NDON') then
            if (.not.qread(qdon)) then
             qread(qdon)=.true.
             read(cmdline,*) ndon
             which=qdon
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:4).eq.'NACC') then
            if (.not.qread(qacc)) then
             qread(qacc)=.true.
             read(cmdline,*) nacc
             which=qacc
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:4).eq.'NNB') then
            if (.not.qread(qnbex)) then
             qread(qnbex)=.true.
             read(cmdline,*) nnbx
             which=qnbex
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:4).eq.'NGRP') then
            if (.not.qread(qgrp)) then
             qread(qgrp)=.true.
             read(cmdline,*) ngrp, nst2
             which=qgrp
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:5).eq.'MOLNT') then
            if (.not.qread(qmolnt)) then
             qread(qmolnt)=.true.
             read(cmdline,*) nmolnt
             which=qmolnt
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:5).eq.'NUMLP') then
            if (.not.qread(qlp)) then
             qread(qlp)=.true.
             read(cmdline,*) nlp, nlph
             which=qlp
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (keyword(1:7).eq.'NCRTERM') then
            if (.not.qread(qcmap)) then
             qread(qcmap)=.true.
             read(cmdline,*) ncmap
             which=qcmap
             cycle
            else
             call error(whoami,'DUPLICATE ENTRY IN STRUCTURE FILE. ABORTING.',0)
             call error(whoami,cmdline(1:l),-1)
             return
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           endif
         endif ! i.gt.0
!
         if (which.lt.0) then
          which=max(which,-which)+1
          if (which.eq.qatom) then
            read(cmdline,*) natom
            cycle
          elseif (which.eq.qbond) then
            read(cmdline,*) nbond
            cycle
          elseif (which.eq.qangl) then
            read(cmdline,*) nangl
            cycle
          elseif (which.eq.qdihe) then
            read(cmdline,*) ndihe
            cycle
          elseif (which.eq.qimpr) then
            read(cmdline,*) nimpr
            cycle
          elseif (which.eq.qdon) then
            read(cmdline,*) ndon
            cycle
          elseif (which.eq.qacc) then
            read(cmdline,*) nacc
            cycle
          elseif (which.eq.qnbex) then
            read(cmdline,*) nnbx
            cycle
          elseif (which.eq.qgrp) then
            read(cmdline,*) ngrp, nst2
            cycle
          elseif (which.eq.qmolnt) then
            read(cmdline,*) nmolnt
            cycle
          elseif (which.eq.qlp) then
            read(cmdline,*) nlp, nlph
            cycle
          elseif (which.eq.qcmap) then
            read(cmdline,*) ncmap
            cycle
          endif
         endif
!

! many psf entries are skipped for now
         select case(which)
          case(qatom);
           if (.not.xplor) then
            read(cmdline,*) i, segid, resid, resname, aname, typeid, charge, mass ! skip the rest: imove, ECH, EHA
            k=atomlist_uadd_ch(atoms, i, segid, resid, resname, aname, typeid, charge, mass)
           else
            read(cmdline,*) i, segid, resid, resname, aname, type, charge, mass ! skip the rest: imove, ECH, EHA
            k=atomlist_uadd_xp(atoms, i, segid, resid, resname, aname, type, charge, mass)
           endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          case(qbond);
           i=numword(cmdline) ! get the number of blank-separated words
           if (mod(i,2).ne.0) then ! not divisible by 2 ==> some bond specification invalid
            call error(whoami, 'INCOMPLETE BOND SPECIFICATION. ABORTING.',0)
            call error(whoami, cmdline(1:l),-1)
           else
            if (allocated(idummy)) deallocate(idummy)
            allocate(idummy(i))
            read(cmdline,*) (idummy(j), j=1,i)
            do j=1,i,2
             k=toplist_uadd(blist,3,(/ idummy(j), idummy(j+1), -1 /))
            enddo
           endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          case(qangl);
           i=numword(cmdline) ! get the number of blank-separated words
           if (mod(i,3).ne.0) then ! not divisible by 3 ==> some angle specification invalid
            call error(whoami, 'INCOMPLETE ANGLE SPECIFICATION. ABORTING.',0)
            call error(whoami, cmdline(1:l),-1)
           else
            if (allocated(idummy)) deallocate(idummy)
            allocate(idummy(i))
            read(cmdline,*) (idummy(j), j=1,i)
            do j=1,i,3
             k=toplist_uadd(alist,4,(/ idummy(j), idummy(j+1), idummy(j+2), -1 /))
            enddo
           endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          case(qdihe, qimpr);
           i=numword(cmdline) ! get the number of blank-separated words
           if (mod(i,4).ne.0) then ! not divisible by 3 ==> some angle specification invalid
            call error(whoami, 'IMCOMPLETE DIHEDRAL ANGLE SPECIFICATION. ABORTING.',0)
            call error(whoami, cmdline(1:l),-1)
           else
            if (allocated(idummy)) deallocate(idummy)
            allocate(idummy(i))
            read(cmdline,*) (idummy(j), j=1,i)
            do j=1,i,4
             if (which.eq.qdihe) then
              k=toplist_uadd(dlist,5,(/ idummy(j), idummy(j+1), idummy(j+2), idummy(j+3), -1 /)) ! add regular dihedral entry
             else
              k=toplist_uadd(ilist,5,(/ idummy(j), idummy(j+1), idummy(j+2), idummy(j+3), -1 /)) ! add improper entry
             endif
            enddo
           endif ! mod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          case(qdon); cycle
          case(qacc); cycle
          case(qnbex); cycle
          case(qgrp); cycle
          case(qmolnt); cycle
          case(qlp); cycle
          case(qcmap); cycle
          case default
           call error(whoami, 'UNRECOGNIZED LINE IN STRUCTURE FILE. ABORTING.',0)
           call error(whoami, cmdline(1:l),-1)
         end select
!
        else ! end of file
         exit
        endif ! ioerr
       enddo ! over all lines in the file
!
       if (allocated(idummy)) deallocate(idummy)
! compare number of entries read with the corresponding number indicated
       if (natom.ne.atoms%last) write(fout,'(2(A,I5))') 'WARNING('//whoami//') : EXPECTED ',natom,' ATOMS, BUT FOUND ',atoms%last
       if (nbond.ne.blist%last) write(fout,'(2(A,I5))') 'WARNING('//whoami//') : EXPECTED ',nbond,' ATOMS, BUT FOUND ',blist%last
       if (nangl.ne.alist%last) write(fout,'(2(A,I5))') 'WARNING('//whoami//') : EXPECTED ',nangl,' ATOMS, BUT FOUND ',alist%last
       if (ndihe.ne.dlist%last) write(fout,'(2(A,I5))') 'WARNING('//whoami//') : EXPECTED ',ndihe,' ATOMS, BUT FOUND ',dlist%last
       if (nimpr.ne.ilist%last) write(fout,'(2(A,I5))') 'WARNING('//whoami//') : EXPECTED ',nimpr,' ATOMS, BUT FOUND ',ilist%last
! other PSF entries are omitted at present
       call message(whoami, 'Structure file read.')
!
       end subroutine psf_read
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine psf_print( &



     & )
       use output
!



!
       implicit none
!





!
       character(len=9), parameter :: whoami='PSF_PRINT'
       integer :: i
!



!
       if (.not.psf_initialized) call psf_init()
!



        call message(whoami,'THE FOLLOWING ATOMS ARE PRESENT')
        call message(whoami,'===============================')
!
        do i=1,atoms%last
         write(fout, '(I10,5(1X,A8), I4, 2G14.6)') &
& atoms%atomid(i),atoms%segid(i),atoms%resid(i),atoms%resname(i),atoms%aname(i),&
& atoms%type(i),atoms%typeid(i),atoms%charge(i),atoms%mass(i)
        enddo
        call message(whoami,'===============================')
        call message(whoami,'THE FOLLOWING BONDS ARE PRESENT')
!
        do i=1,blist%last
         write(fout, '(I6," -- ",I6,"; ")', advance='no') blist%ind(1,i),blist%ind(2,i)
         if (mod(i,6).eq.0) write(fout,*)
        enddo
        write(fout,*)
        call message(whoami,'================================')
        call message(whoami,'THE FOLLOWING ANGLES ARE PRESENT')
!
        do i=1,alist%last
         write(fout, '(I6," -- ",I6," -- ",I6,"; ")', advance='no') alist%ind(1,i),alist%ind(2,i),alist%ind(3,i)
         if (mod(i,5).eq.0) write(fout,*)
        enddo
        write(fout,*)
        call message(whoami,'=========================================')
        call message(whoami,'THE FOLLOWING DIHEDRAL ANGLES ARE PRESENT')
!
        do i=1,dlist%last
         write(fout, '(3(I6," -- "),I6,"; ")', advance='no') dlist%ind(1,i),dlist%ind(2,i),dlist%ind(3,i),dlist%ind(4,i)
         if (mod(i,3).eq.0) write(fout,*)
        enddo
        write(fout,*)
        call message(whoami,'==================================================')
        call message(whoami,'THE FOLLOWING IMPROPER DIHEDRAL ANGLES ARE PRESENT')
!
        do i=1,ilist%last
         write(fout, '(3(I6," -- "),I6,"; ")', advance='no') ilist%ind(1,i),ilist%ind(2,i),ilist%ind(3,i),ilist%ind(4,i)
         if (mod(i,3).eq.0) write(fout,*)
        enddo
        write(fout,*)
! there are other parts of a psf that are omitted at present



       end subroutine psf_print
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       subroutine psf_info( &



     & )
       use output
       use parser, only: tab
!



!
       implicit none
!





!
       character(len=8), parameter :: whoami='PSF_INFO'
       integer :: i
!



!
       if (.not.psf_initialized) call psf_init()
!



        call message(whoami,'Structure summary:')
        write(fout,'(A)') tab//'================='
        write(fout,'(A,I6,A)' ) tab,atoms%last,' ATOMS'
        write(fout,'(A,I6,A)' ) tab,blist%last,' BONDS'
        write(fout,'(A,I6,A)' ) tab,alist%last,' ANGLES'
        write(fout,'(A,I6,A)' ) tab,dlist%last,' DIHEDRALS'
        write(fout,'(A,I6,A)' ) tab,ilist%last,' IMPROPERS'
        write(fout,'(A)') tab//'================='
        write(fout,'(A,G14.6)' ) tab//'Total charge: ',sum(atoms%charge(1:atoms%last))
        write(fout,'(A,G16.8)' ) tab//'Total mass  : ',sum(atoms%mass(1:atoms%last))
        write(fout,'(A)') tab//'================='



! other interactions are currently ignored
       end subroutine psf_info
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function psf_natom result(natom)
       implicit none
       integer :: natom
       if (psf_initialized) then
        natom=atoms%last
       else
        natom=0
       endif
       end function psf_natom
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end module psf
