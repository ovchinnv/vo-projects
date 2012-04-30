 program bfit

 use bestfit
 use constants

 float :: A(3,3)
 
 float :: eval(3), evec(3,3), b0(3,3)
 
! ok
 A=reshape ( (/ 4.1d0, ERRTOL(), 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 3d0/),(/3,3/) )
! ok
 A=reshape ( (/ 4.1d0, ERRTOL(), ERRTOL(), 0d0, ERRTOL(), ERRTOL(), ERRTOL(), 0d0, ERRTOL()/),(/3,3/) )
! ok
 A=reshape ( (/  ERRTOL(), 1d0, 0d0, & 
                 1d0,      0d0, 3d0, & 
                 0d0,      3d0, ERRTOL()*1000/),(/3,3/) )
! A=reshape ( (/  ERRTOL(), 0., 0., 1., 0., 0., 0., 3., ERRTOL()*1000/),(/3,3/) )
! A=reshape ( (/  ERRTOL(), 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, ERRTOL()*1000/),(/3,3/) )


! this failed:
! in this case the eigenvalues were very close to zero (but computed as 1e-8), which lead to two 'duplicate' evectors; fixed
! A=reshape( (/ &
! 0.32753820909748578016352d-01, 0.75997504774541049799019d-01, 0.76284169485594674475770d-01,&
! 0.75997504774541049799019d-01, 0.17633425876849012636782d+00, 0.17699939651858714873001d+00,&
! 0.76284169485594674475770d-01, 0.17699939651858714873001d+00, 0.17766704318685860863880d+00 /), (/3,3/) )
!
! this also failed:

! A=reshape( (/ &
! 0.27427009764448101458711d+00, 0.10798602617951455484580d+00,-0.16209868811080687251457d-06,&
! 0.10798602617951455484580d+00, 0.42516417029020051043808d-01,-0.63821733861371916309836d-07,&
!-0.16209868811080687251457d-06,-0.63821733861371916309836d-07, 0.95803315464978367580174d-13 /), (/3,3/) )
!
! fixed by using the fact that eigenvectors coresponding to different EV of a symmetric matrix are orthogonal
 
 
 call eig3s(A, eval, evec)
 
 write(0,*) eval

 write(0,'(3G22.14)') evec
 
 
! check diagonalization directly
        b0(:,1)=eval(1)*evec(:,1);b0(:,2)=eval(2)*evec(:,2);b0(:,3)=eval(3)*evec(:,3)! aa
        write(0,*)' A - U x M x U^T (should be 0):'
        write(0,'(3E22.15)') matmul(b0,transpose(evec))-A ! aa

 
 
 end