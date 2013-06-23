* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*            A D A M S  -  B A S H F O R T H   S T E P            *
*                                                                 *
*    All necessary routines to perform an A-B step are included   *
*                                                                 *
*           Predictor step  :   USTARAB                           *
*           Divergence      :   DIVAB                             *
*           Velocity corect.:   CORAB                             *
*           Store prev. val.:   STOREAB                           *
*                                                                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c V.O 8/11/03.  Implicit mode solves for velocity _correction_ (delta form) 

      SUBROUTINE predictor(uo,vo,wo,p,uab1,vab1,wab1,
     &                     uab2,vab2,wab2,ucn,vcn,wcn,tv,
     &                     us,vs,ws,xf,dpdx,dpdy,dpdz,dtm,a1,a2,botbc,topbc,
     &                     im,jm,km                           )
c
      implicit none
      include 'headers/common.h'
      integer im,jm,km,i,j,k,ii
      real tvip,tvim,delx
c
      REAL   uo(im,jm,km),vo(im,jm,km),wo(im,jm,km),
     &       p(im,jm,km),
     &       uab1(im,jm,km),vab1(im,jm,km),wab1(im,jm,km),
     &       uab2(im,jm,km),vab2(im,jm,km),wab2(im,jm,km),
     &       ucn(im,jm,km),vcn(im,jm,km),wcn(im,jm,km),
     &       tv(im,jm,km),
     &       us(im,jm,km),vs(im,jm,km),ws(im,jm,km)
      REAL   xf(im)
c
      REAL   a(mmx-2),b(mmx-2),c(mmx-2)
c
      REAL   dtm,dpdx,dpdy,dpdz,a1,a2
      REAL   coef
c
      REAL   botbc(3,jm,km),topbc(3,jm,km)
c
**** Compute provisional values for each momentum eq.
*
c      for implicit part
      coef=   0.5*dtm

c------------------------------------------------------------------
c                                                    compute RHS
c------------------------------------------------------------------
c     uab1,vab1,wab1 store the rhs of the mom. equation since they will 
c     not be used after this

      do k=kz1+1,kz2
         do j=jy1+1,jy2
            do i=ix1+1,ix2
               uab1(i,j,k)=
     &              (a2*uab2(i,j,k)-a1*uab1(i,j,k)  &  !
     &              -au(i)*(p(i+1,j,k)-p(i,j,k))/dx
     &              -dpdx      
     &              +ucn(i,j,k))*dtm
               vab1(i,j,k)=
     &              (a2*vab2(i,j,k)-a1*vab1(i,j,k)
     &              -bv(j)*(p(i,j+1,k)-p(i,j,k))/dy
     &              -dpdy
     &              +vcn(i,j,k))*dtm
               wab1(i,j,k)=
     &              (a2*wab2(i,j,k)-a1*wab1(i,j,k)
     &              -cw(k)*(p(i,j,k+1)-p(i,j,k))/dz
     &              -dpdz
     &              +wcn(i,j,k))*dtm
            enddo
         enddo
      enddo
      if (impfla.eq.1) then 
c     x-mom
      do k=kz1+1,kz2
       do j=jy1+1,jy2
        do i=ibu,ieu
          ii=i-1
	  tvip=ru1+2.*tv(i+1,j,k)
	  tvim=ru1+2.*tv(i  ,j,k)
	  a(ii)=-coef*au(i)*app(i)*tvim
	  c(ii)=-coef*au(i)*app(i+1)*tvip
	  b(ii)=1.-a(ii)-c(ii)
        enddo !i
c         botbc
        delx=xf(2)-xf(1)
        b(1)= b(1)*ibc(1,1)           !dirichlet
     &      +(a(1)+b(1))*(1.-ibc(1,1))!neumann
        uab1(2,j,k)=(uab1(2,j,k)-a(1)*     botbc(1,j,k))*ibc(1,1) !dirichlet
     &             +(uab1(2,j,k)+a(1)*delx*botbc(1,j,k))*(1.-ibc(1,1)) !neumann 1st o.
c         topbc
        delx=xf(ix2)-xf(ix2-1)	  
        b(ieu-1)= b(ieu-1)*ibc(1,2)               !dirichlet
     &          +(c(ieu-1)+b(ieu-1))*(1.-ibc(1,2))!neumann
        uab1(ieu,j,k)=(uab1(ieu,j,k)-c(ieu-1)*     topbc(1,j,k))*ibc(1,2) !dirichlet
     &               +(uab1(ieu,j,k)-c(ieu-1)*delx*topbc(1,j,k))*(1.-ibc(1,2)) !neumann 1st o.
        call tridag(a,b,c,uab1(2,j,k),uab1(2,j,k),ieu-1,mmx-2)
       enddo !j
      enddo !k	
c     y-mom        
      do k=kz1+1,kz2
       do j=jbv,jev
        do i=ix1+1,ix2
          ii=i-1
	  tvip=ru1+0.25*(tv(i+1,j,k)+tv(i+1,j+1,k)
     & 	               + tv(i  ,j,k)+tv(i  ,j+1,k))
	  tvim=ru1+0.25*(tv(i  ,j,k)+tv(i  ,j+1,k)
     & 	               + tv(i-1,j,k)+tv(i-1,j+1,k))
	  a(ii)=-coef*av(i)*au(i-1)*tvim
	  c(ii)=-coef*av(i)*au(i  )*tvip
	  b(ii)=1.-a(ii)-c(ii)
c         botbc
        enddo !i
        delx=xf(2)-xf(1)
        b(1)=(b(1)-a(1))*ibc(2,1)           !dirichlet
     &      +(a(1)+b(1))*(1.-ibc(2,1))      !neumann
        vab1(2,j,k)=(vab1(2,j,k)-2.*a(1)*     botbc(2,j,k))*ibc(2,1) !dirichlet
     &             +(vab1(2,j,k)+   a(1)*delx*botbc(2,j,k))*(1.-ibc(2,1)) !neumann 1st o.
c         topbc
        delx=xf(ix2+1)-xf(ix2)	  
        b(ix2-1)= (b(ix2-1)-c(ix2-1))*ibc(2,2)     !dirichlet
     &          + (c(ix2-1)+b(ix2-1))*(1.-ibc(2,2))!neumann
        vab1(ix2,j,k)=(vab1(ix2,j,k)-2.*c(ix2-1)*  topbc(2,j,k))*ibc(2,2) !dirichlet
     &               +(vab1(ix2,j,k)-c(ix2-1)*delx*topbc(2,j,k))*(1.-ibc(2,2)) !neumann 1st o.
        call tridag(a,b,c,vab1(2,j,k),vab1(2,j,k),ix2-1,mmx-2)
       enddo !j
      enddo !k	
c     z-mom        
      do k=kbw,kew
       do j=jy1+1,jy2
        do i=ix1+1,ix2
          ii=i-1
	  tvip=ru1+0.25*(tv(i+1,j,k)+tv(i+1,j,k+1)
     & 	               + tv(i  ,j,k)+tv(i  ,j,k+1))
	  tvim=ru1+0.25*(tv(i  ,j,k)+tv(i  ,j,k+1)
     & 	               + tv(i-1,j,k)+tv(i-1,j,k+1))
	  a(ii)=-coef*aw(i)*au(i-1)*tvim
	  c(ii)=-coef*aw(i)*au(i  )*tvip
	  b(ii)=1.-a(ii)-c(ii)
c         botbc
        enddo !i
        delx=xf(2)-xf(1)
        b(1)=(b(1)-a(1))*ibc(3,1)           !dirichlet
     &      +(a(1)+b(1))*(1.-ibc(3,1))       !neumann
        wab1(2,j,k)=(wab1(2,j,k)-2.*a(1)*     botbc(3,j,k))*ibc(3,1) !dirichlet
     &             +(wab1(2,j,k)+   a(1)*delx*botbc(3,j,k))*(1.-ibc(3,1)) !neumann 1st o.
c         topbc
        delx=xf(ix2+1)-xf(ix2)	  
        b(ix2-1)= (b(ix2-1)-c(ix2-1))*ibc(3,2)               !dirichlet
     &          + (c(ix2-1)+b(ix2-1))*(1.-ibc(3,2))!neumann
        wab1(ix2,j,k)=(wab1(ix2,j,k)-2.*c(ix2-1)*  topbc(3,j,k))*ibc(3,2) !dirichlet
     &               +(wab1(ix2,j,k)-c(ix2-1)*delx*topbc(3,j,k))*(1.-ibc(3,2)) !neumann 1st o.
        call tridag(a,b,c,wab1(2,j,k),wab1(2,j,k),ix2-1,mmx-2)
       enddo !j
      enddo !k	
      endif
c
c     now compute the velocities
c------------------------------------------------------------------
c                                                    compute ustar
c------------------------------------------------------------------

      do k=kz1+1,kz2
         do j=jy1+1,jy2
            do i=ibu,ieu
	      us(i,j,k)=uo(i,j,k)+uab1(i,j,k)
            enddo
         enddo
      enddo
c------------------------------------------------------------------
c                                                    compute vstar
c------------------------------------------------------------------

      do k=kz1+1,kz2
         do j=jbv,jev
            do i=ix1+1,ix2
	      vs(i,j,k)=vo(i,j,k)+vab1(i,j,k)
            enddo
         enddo
      enddo
c------------------------------------------------------------------
c                                                    compute wstar
c------------------------------------------------------------------
      do k=kbw,kew
         do j=jy1+1,jy2
            do i=ix1+1,ix2
	      ws(i,j,k)=wo(i,j,k)+wab1(i,j,k)
            enddo
         enddo
      enddo

  999 continue 
      RETURN
      END    

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE DIVERGENCE(US,VS,WS,IM,JM,KM,DTM,DIV)
*                                                                 *
*          Computes the divergence for an AB step                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      DIMENSION US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),
     &          DIV(IM,JM,KM)
      real rflag,coef
c
**** Compute divergence
*
c
      rflag=real(impfla)
      coef = 1.+ rflag

      do 10 i=ix1+1,ix2
       do 20 j=jy1+1,jy2
        do 30 k=kz1+1,kz2
         div(i,j,k)=coef*(ap(i)*(us(i,j,k)-us(i-1,j,k))/dx
     &              +bp(j)*(vs(i,j,k)-vs(i,j-1,k))/dy
     &              +cp(k)*(ws(i,j,k)-ws(i,j,k-1))/dz)
     &               /(dtm)
30      continue
20     continue
10    continue


      RETURN
      END

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE CORRECT(US,VS,WS,DPP,DTM,IM,JM,KM,UUO,VVO,WWO,
     &                 IT,RSU,RSV,RSW)
*                                                                 *
*                Performs an AB corector step                     *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      DIMENSION UUO(IM,JM,KM),VVO(IM,JM,KM),WWO(IM,JM,KM),
     &          US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),
     &          DPP(IM,JM,KM)
      real rflag,coef 

      rflag=real(impfla)
      coef = 1./(1.+ rflag)

*
**** Update u values
*
        do  k=kz1+1,kz2
           do  j=jy1+1,jy2
              do  i=ibu,ieu


                 uuo(i,j,k)=us(i,j,k)-coef*au(i)*( dpp(i+1,j,k)-dpp(i,j,k) )
     &                             *(dtm/dx )
              enddo
           enddo
        enddo
*
**** Update v values
*
                do  k=kz1+1,kz2
                   do  j=jbv,jev
                      do  i=ix1+1,ix2


                         vvo(i,j,k)=vs(i,j,k)-coef*bv(j)*
     &                        ( dpp(i,j+1,k)-dpp(i,j,k) )
     &                        *(dtm/dy )
                      enddo
                   enddo
                enddo

*
**** Update w values
*

        do k=kbw,kew
           do  j=jy1+1,jy2
              do i=ix1+1,ix2


                 wwo(i,j,k)=ws(i,j,k)-coef*cw(k)*( dpp(i,j,k+1)-dpp(i,j,k) )
     &                             *(dtm/dz )
         
      enddo
      enddo
      enddo

      RETURN
      END


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE STOREN(UUAB1,VVAB1,WWAB1,UUAB2,VVAB2,WWAB2,
     &                   tmeth,IM,JM,KM                      )
*                                                                 *
*              Stores RHS + Presure gradient                      *                
*            to be used in the predictor step                     *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      DIMENSION UUAB1(IM,JM,KM),VVAB1(IM,JM,KM),WWAB1(IM,JM,KM),
     &          UUAB2(IM,JM,KM),VVAB2(IM,JM,KM),WWAB2(IM,JM,KM)
c
       REAL     tmeth
*
**** Store RHS + presure gradient for 3 momentum equations
*
        do  k=kz1+1,kz2
           do  j=jy1+1,jy2
              do  i=ibu,ieu


                uuab1(i,j,k)= (1.-tmeth)*uuab1(i,j,k)+tmeth*uuab2(i,j,k)
              enddo
           enddo
        enddo

*


        do  k=kz1+1,kz2
           do  j=jbv,jev
              do  i=ix1+1,ix2
                 vvab1(i,j,k)=(1.-tmeth)*vvab1(i,j,k)+tmeth*vvab2(i,j,k)
              enddo
           enddo
        enddo
*
        do  k=kbw,kew
           do  j=jy1+1,jy2
               do  i=ix1+1,ix2
                 wwab1(i,j,k)=(1.-tmeth)*wwab1(i,j,k)+tmeth*wwab2(i,j,k)
               enddo
            enddo
         enddo
      
      RETURN
      END




