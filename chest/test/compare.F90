
 implicit none
 float, dimension(66,66,66)  :: uanl, uches
 float :: grid(66,66,66,3)
 float :: x1(66), x2(66), y1(66), y2(66), z1(66), z2(66)
 float :: griderror, error
 integer :: where(3)
!
! compare grids
!
 open(1,file='xg_test.dat',form='formatted', status='old')
 read(1,*) x1
 close(1)
 open(1,file='yg_test.dat',form='formatted', status='old')
 read(1,*) y1
 close(1)
 open(1,file='zg_test.dat',form='formatted', status='old')
 read(1,*) z1
 close(1)
! bring grids to center
 x2(2:66)=0.5d0*(x1(1:66-1)+x1(2:66)); x2(1)=x1(1)-0.5d0*(x1(2)-x1(1)); x1=x2;
 y2(2:66)=0.5d0*(y1(1:66-1)+y1(2:66)); y2(1)=y1(1)-0.5d0*(y1(2)-y1(1)); y1=y2;
 z2(2:66)=0.5d0*(z1(1:66-1)+z1(2:66)); z2(1)=z1(1)-0.5d0*(z1(2)-z1(1)); z1=z2;
!
! grid output from CHEST is at center
!
 open(1,file='xyz_test.xyz',form='formatted', status='old')
 read(1,*)
 read(1,*) grid
 close(1)
!
 x2=grid(:,1,1,1); y2=grid(1,:,1,2); z2=grid(1,1,:,3)
 griderror=max(maxval(abs(x1-x2)), maxval(abs(y1-y2)), maxval(abs(z1-z2)))
! write(0,*) x1-x2
! write(0,*) y1-y2
! write(0,*) z1-z2
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
 error=maxval(uches(2:66-1,2:66-1,2:66-1))
 where=maxloc(uches(2:66-1,2:66-1,2:66-1))
 write(6,'(2F25.12,3I10)') griderror, error, where
 end 
