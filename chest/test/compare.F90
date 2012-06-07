
 implicit none
 float, dimension(258,258)  :: uanl, uches
 float :: grid(258,258,2)
 float :: x1(258), x2(258), y1(258), y2(258)
 float :: griderror, error
 integer :: where(2)
!
! compare grids
!
 open(1,file='xg_test.dat',form='formatted', status='old')
 read(1,*) x1
 close(1)
 open(1,file='yg_test.dat',form='formatted', status='old')
 read(1,*) y1
 close(1)
! bring grids to center
 x2(2:258)=0.5d0*(x1(1:258-1)+x1(2:258)); x2(1)=x1(1)-0.5d0*(x1(2)-x1(1)); x1=x2;
 y2(2:258)=0.5d0*(y1(1:258-1)+y1(2:258)); y2(1)=y1(1)-0.5d0*(y1(2)-y1(1)); y1=y2;
!
! grid output from CHEST is at center
!
 open(1,file='xy_test.xyz',form='formatted', status='old')
 read(1,*)
 read(1,*) grid
 close(1)
!
 x2=grid(:,1,1); y2=grid(1,:,2)
 griderror=max(maxval(abs(x1-x2)), maxval(abs(y1-y2)))
! write(0,*) x1-x2
! write(0,*) y1-y2
! write(0,*) y2
!
! compare solutions
!
 open(1,file='uexact_test.dat',form='formatted', status='old')
 read(1,*)
 read(1,*) uanl
 close(1)
!
 open(1,file='solution.dat',form='formatted', status='old')
 read(1,*)
 read(1,*)
 read(1,*) uches
 close(1)
!
 uches=abs(uanl-uches);
 error=maxval(uches(2:258-1,2:258-1))
 where=maxloc(uches(2:258-1,2:258-1))
 write(6,'(2F25.12,2I10)') griderror, error, where
 end 
