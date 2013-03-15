C
C-----------------------------------------------------------------------
C                 ***************************                         
C                 *         info.f          *                        
C                 ***************************                       
C----------------------------------------------------------------------- 
C
C	- SCRINFO   : write screen infos
C	- CALCME    : comp. mass error
C	- CALCTKE   : comp. turb. kinetic energy
C	- CALCDIV   : comp. max. divergence
C	- CALCTWALL : comp. max. twall
C	- CALCD1    : comp. displ. thickness
C	- PRINTVAR  : "left-over" from former times
C       - SCRVEL
C	- SCRRHS
C	- SCRUVWP
C	- SCRCS
C
C----------------------------------------------------------------------- 
C
C 
C-----SUBROUTINE-SRCINFO------------------------P. FLOHR--23/01/1994----
C
      SUBROUTINE SCRINFO(uo,vo,wo,uav,vav,wav,DisplHist,TwallHist,
     &                   qold, dpdx, tke, divmax,divmed,
     &                   imx,jmx,kmx,erm,amass,yplusmin,
     &                   xf, yf, zf,icycle, dtm,Nstep, Miter,
     &                   nx, ny, nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - print infos to the screen
C                                                                       
C-----------------------------------------------------------------------
C 
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
C.....common variables
C
      real uo(nx,ny,*),vo(nx,ny,*),wo(nx,ny,*),
     &          uav(nx,ny,*), vav(nx,ny,*),wav(nx,ny,*),
     &          xf(*), yf(*), zf(*),
     &          DisplHist(Nstep,Miter), TwallHist(Nstep,Miter)
      INTEGER   nx, ny, nz, kmx, imx, jmx
      REAL      qold, dpdx, divmax, tke, tlevel, dtm, divmed, yplusmin
      integer Nstep, Miter, icycle
      real amass, erm
C
C.....local variables
C
      CHARACTER*(10) c(3)
      
      INTEGER i,date(8)
C
C-----------------------------------------------------------------------
C     GET TIME AND DATE
      IF(MYRANK.EQ.0) THEN
c         CALL DATE_AND_TIME(C(1),C(2),C(3),DATE)
C
         WRITE(6,*) '*...................................'
         WRITE(6,*) '*..iteration nr:                       ', ICYCLE
         WRITE(6,*) '*..time step      = ',DTM
         WRITE(6,*) '*..total time     = ',Tlevel
c         WRITE(6,100) DATE(1),DATE(2),DATE(3),C(3)
c         WRITE(6,101) DATE(5),DATE(6),DATE(7),DATE(4)
         WRITE(6,*) '*..nx ny nz       = ',nx, ny, nz
         WRITE(6,*) '*..x  y  z        = ',xf(ix2), yf(jy2), zf(kz2)
         WRITE(6,*) '*..mass error(max)= ',erm
         WRITE(6,*) '*..mass error(L2) = ',amass
         WRITE(6,*) '*..max. div.      = ',divmax,"  (",imx,jmx,kmx,")"
         WRITE(6,*) '*..div (L2)       = ',divmed
         WRITE(6,*) '*..tke            = ',tke
         WRITE(6,*) '*..mass flow      = ',qold
         WRITE(6,*) '*..dpdx           = ',dpdx
         WRITE(6,*) '*..z+ (min.)      = ',yplusmin
         WRITE(6,*) '*.........................'
C     
         DO 10 i = 1, Nstep
            WRITE(6,*) '*..delta1(x)   = ', DisplHist(i,icycle), Xpos(i)
 10      CONTINUE
C     
         DO 20 i = 1, Nstep
            WRITE(6,*) '*..tw (low.)(x)= ', TwallHist(i,icycle), Xpos(i)
 20      CONTINUE
C
      ENDIF
 100  FORMAT(1x,'*..machine day',2x,i4,'-',i2,'-',i2,
     &     ' (YYYY-MM-DD)',a10)
 101  FORMAT(1x,'*..machine time',2x,i2,':',i2,':',i2,' (HH:MM:SS)')

      RETURN
      END
C
C 
C-----SUBROUTINE-CALCME-------------------------P. FLOHR--23/02/1994----
C
      SUBROUTINE CALCME(UUO,VVO,WWO,XF,YF,ZF,dtm,IM,JM,KM)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - mass flow error
C
C-----------------------------------------------------------------------
C 
      INCLUDE "mpif.h"
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real UUO(IM,JM,KM),VVO(IM,JM,KM),WWO(IM,JM,KM),
     &          XF(*),YF(*),ZF(*)
      REAL      erm,er,amass,dtm
      integer   iim,jjm,kkm
      INTEGER   nx, ny, nz
      real in(2),out(2), dyy,dxx,dzz, div
C
      INTEGER   i, j, k,mytype1,mytype2
      integer im,jm,km,index, mrm1,ier,istop
C
C-----------------------------------------------------------------------
C
#ifdef SINGLE
      mytype1=MPI_2REAL
      mytype2=MPI_REAL
#endif
#ifdef DOUBLE
      mytype1=MPI_2DOUBLE_PRECISION
      mytype2=MPI_DOUBLE_PRECISION
#endif

      er =0.0
      erm=0.0
      amass = 0.0
C
 
      istop=ix2

       do 10 k=kz1+1,kz2
         dzz=zf(k)-zf(k-1)
         do 10 j=jy1+1,jy2
            dyy=yf(j)-yf(j-1)
         
            do 10 i=ix1+1,istop
               dxx=xf(i)-xf(i-1)
               div=(uuo(i,j,k)-uuo(i-1,j,k))*dyy*dzz+
     &              (vvo(i,j,k)-vvo(i,j-1,k))*dxx*dzz+
     &              (wwo(i,j,k)-wwo(i,j,k-1))*dxx*dyy
               
               er=abs(div*dtm)/(dxx*dyy*dzz)
               amass = amass + er**2
               erm   = amax1(erm,er)
         if(erm.eq.er) then
            iim=i
            jjm=j
            kkm=k
         endif
 10   CONTINUE

      
      in(1)=erm
      in(2)=iim+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MAXLOC,0,MPI_COMM_HYDRO,ier)

      erm=out(1)
      index=out(2)
      
      mrm1=index/(im*jm*km)
      kkm=(index-(im*jm*km)*mrm1)/(im*jm)
      jjm=(index-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(index-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jjm=jjm+1

      kkm=kkm+1+mrm1*(km-2)
      IF(MYRANK.EQ.0) then
         write(6,1955) erm,iim,jjm,kkm
      endif
      
      CALL MPI_REDUCE(AMASS,ERM,1,MYTYPE2,MPI_SUM,0,MPI_COMM_HYDRO,ier)

      IF(MYRANK.EQ.0) then
         amass = sqrt(erm/(ix2*jy2*kz2*mysize))
         write(6,1956) amass
      endif
 1955 format(' *...max mass error  = ',1pe15.5,' (',3i4,')')
 1956 format(' *...avg mass error  = ',1pe15.5)
      RETURN
      END
C
C 
C-----SUBROUTINE-CALCTKE------------------------P. FLOHR--23/01/1994----
C
      SUBROUTINE CALCTKE(tke,uo,vo,wo,uav,vav,wav,zf,x,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - calc. averaged turbulent kinetic energy
C                                                                       
C-----------------------------------------------------------------------
C 
      include 'headers/common.h'
C
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     &          uav(nx,nz), vav(nx,nz), wav(nx,nz),
     &          x(nx)
      REAL      zf(nz)
      INTEGER   nx,ny,nz
      REAL      tke
C
      INTEGER   i,j,k, imin, imax
      REAL      ufl, vfl, wfl
      REAL      ubar(mmz), vbar(mmz), wbar(mmz)
      REAL      ucord(mmz)
      REAL      coef1, coef2, coef3, coef4, 
     *          dudzb, fijp, h1, h2, h3, h4, utaub, utaut
C 
C-----------------------------------------------------------------------
C
C
C
      fijp  = 1.0/real((ix2 - (ix1+1) + 1)*(jy2 - (jy1+1) + 1))
C
c     ucord(1) = zf(1)
c     DO 3 k = 2, kz2
C       Initial stuff.
c       ucord(k) = 0.5*(zf(k)+zf(k-1))
c       ubar(k)  = 0.0
c       vbar(k)  = 0.0
c       wbar(k)  = 0.0
C       Mean profile.
c       DO 5 j=2, jy2
c         DO 7 i=2, ix2
c           ubar(k) = ubar(k) + uo(i,j,k)
c           vbar(k) = vbar(k) + vo(i,j,k)
c           wbar(k) = wbar(k) + wo(i,j,k)
c   7     CONTINUE
c   5   CONTINUE
c       ubar(k) = ubar(k)*fijp
c       vbar(k) = vbar(k)*fijp
c       wbar(k) = wbar(k)*fijp
c   3 CONTINUE
c
c...<u_t>
c
c      h1 =ucord(2)
c      h2 =ucord(3)
c      coef1 = h2/(h1*(h2-h1))
c      coef2 = h1/(h2*(h2-h1))
c      dudzb = coef1*ubar(2)-coef2*ubar(3)
c      utaub = sqrt(ru1*dudzb)
c
c      h3 =  ucord(kz2)   - zf(kz2)
c      h4 =  ucord(kz2-1) - zf(kz2)
c      print*,'h',h1,h2,h3,h4
c      coef3 = h4/(h3*(h4-h3))
c      coef4 = h3/(h4*(h4-h3))

c      dudzt = - coef3*ubar(kz2) + coef4*ubar(kz2-1)
c      utaut = sqrt(ru1*dudzt)
c      print*,' '
c      print*,'ubar_b',ubar(2),ubar(3),ubar(4),ubar(5)
c      print*,'ubar_t',ubar(kz2),ubar(kz2-1),ubar(kz2-2),ubar(kz2-3)
c      print*,' '
C
      tke  = 0.
      DO 10 i=2,ix2
        DO 10 j=2,jy2
          DO 10 k=2,kz2
            ufl = .5*(uo(i,j,k)+uo(i-1,j,k))-uav(i,k)
            vfl = .5*(vo(i,j,k)+vo(i,j-1,k))-vav(i,k)
            wfl = .5*(wo(i,j,k)+wo(i,j,k-1))-wav(i,k)
c           ufl = .5*(uo(i,j,k)+uo(i-1,j,k))-ubar(k)
c           vfl = .5*(vo(i,j,k)+vo(i,j-1,k))-vbar(k)
c           wfl = .5*(wo(i,j,k)+wo(i,j,k-1))-wbar(k)
            tke = tke + ufl*ufl+vfl*vfl+wfl*wfl
 10   CONTINUE
      tke = 0.5 * SQRT(tke)   
C

      write(6,*) '*...volume average tke = ',tke
c     write(*,*) '*...shear stress at the walls =', utaub, utaut

      RETURN
      END
C
C 
C-----SUBROUTINE-CALCDIV------------------------P. FLOHR--23/01/1994----
C
      SUBROUTINE CALCDIV(u,v,w,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - calc. max. divergence
C                                                                       
C-----------------------------------------------------------------------
C 
      INCLUDE "mpif.h"
      include 'headers/common.h'
C
      real u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      REAL      divmax,amass,erm
      INTEGER   nx, ny, nz
      
C
      INTEGER   i, j, k,mytype1,mytype2
      REAL      divm,in(2),out(2), divmed
      integer im,jm,km,mrm1,iim,jjm,kkm,index, ier
      
C 
C-----------------------------------------------------------------------
C
#ifdef SINGLE
      mytype1=MPI_2REAL
      mytype2=MPI_REAL
#endif
#ifdef DOUBLE
      mytype1=MPI_2DOUBLE_PRECISION
      mytype2=MPI_DOUBLE_PRECISION
#endif

      im=nx
      jm=ny
      km=nz
      
      iim=0
      jjm=0
      kkm=0

      divmax = 0.
      divmed = 0.
C
      DO 11 i=ix1+1,ix2
        DO 11 j=jy1+1,jy2
          DO 11 k=kz1+1,kz2
            divm = ABS( ap(i)*(u(i,j,k)-u(i-1,j,k))/dx
     &                + bp(j)*(v(i,j,k)-v(i,j-1,k))/dy
     &                + cp(k)*(w(i,j,k)-w(i,j,k-1))/dz )
            divmed = divmed + divm**2
            if(divm .gt. divmax) then
             divmax = divm
             iim=i
             jjm=j
             kkm=k
            endif
 11   CONTINUE
C
      in(1)=divmax
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MAXLOC,0,MPI_COMM_HYDRO,ier)

      erm=out(1)
      index=out(2)
      
      mrm1=index/(im*jm*km)
      kkm=(index-(im*jm*km)*mrm1)/(im*jm)
      jjm=(index-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(index-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jjm=jjm+1
      iim=iim+1
      kkm=kkm+1+mrm1*(km-2)



      IF(MYRANK.EQ.0) then
         write(6,100) erm,iim,jjm,kkm
      endif

      amass=divmed
      CALL MPI_REDUCE(AMASS,ERM,1,MYTYPE2,MPI_SUM,0,MPI_COMM_HYDRO,ier)

      IF(MYRANK.EQ.0) then
         divmed = sqrt(erm/(ix2*jy2*kz2*mysize))
         write(6,110) divmed
      endif

 100  format(' *...max divergence  = ',1pe15.5,' (',3i4,')')
 110  format(' *...avg divergence  = ',1pe15.5)
C
      RETURN
      END
C
C 
C-----SUBROUTINE-CALCTWALL----------------------P. FLOHR--23/02/1994----
C
      SUBROUTINE CALCTWALL(twist,twav,x,Twallhist,Nstep,Miter,
     &                     icycle,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - calc. inst. twall (low.) at diff. positions
C                                                                       
C-----------------------------------------------------------------------
C 
      include 'headers/common.h'
C
      real twist(nx,ny,nz,3), twav(nx,ny,nz),
     &          TwallHist(Nstep, Miter),x(nx)
      INTEGER   nx, ny, nz, Nstep, Miter
C
      INTEGER   i, j, ii,icount_tmp, icycle
      REAL      Twlomax, twloc
C 
C-----------------------------------------------------------------------
C
      DO 10 i = 1, Nstep
        ii = INT(ix2 * Xpos(i)/x(ix2))
        Twlomax = 0.
        twloc = 0.
        DO 101 j=2,jy2
          IF (Twlomax .LT. twist(ii,j,kz1+1,1)) THEN
            Twlomax = twist(ii,j,kz1+1,1)
c          IF (twloc .LT. twav(ii,j,kz1+1)) THEN
c            twloc = twav(ii,j,kz1+1)
          END IF
 101    CONTINUE
        TwallHist(i,icycle) = Twlomax

        icount_tmp = MOD(ICYCLE,INFNUM)
        IF(icount_tmp.EQ.0) THEN
          j=(jy2+jy1)/2.
          WRITE(6,*) '*..twav(x)        = ', twav(ii,j,kz1+1), Xpos(i)
        END IF
 10   CONTINUE
C
      RETURN
      END
C
C 
C-----SUBROUTINE-CALCD1-------------------------P. FLOHR--01/02/1994----
C
      SUBROUTINE CALCD1(DisplHist,u,twav,x,z,Miter,Nstep,
     &                  icycle,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - COMPUTE DISPLACEMENT THICKNESS FROM
C                   AVERAGED PROFILE AT DIFFERENT LOCATIONS
C                                                       
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real u(nx,ny,nz), x(nx), z(nz), twav(nx,ny,nz),
     &          DisplHist(Nstep,Miter)
      INTEGER   nx, ny, nz, miter, nstep
C
      REAL      xloc(25),yloc(25),y2(25)
      REAL      yp1,ypn, dzloc, yspl
      REAL      dzblloc, zblloc1, zblloc2
      REAL      intgrl, intgrl2, intgrl3, twall, xx
      INTEGER   i, j, k, ii, icycle
C 
C-----------------------------------------------------------------------
C
      j     = (jy2-jy1)/2
C
      DO 10 i = 1, Nstep
        ii = INT(ix2 * Xpos(i)/x(ix2))
C
C.....integrate computed velocity-field (trapezoid rule)
C
        intgrl = 0.
        DO 101 k = 1, kz2-10
          intgrl = intgrl + (z(k+1) - z(k)) * 
     &                      (1 - u(ii,j,k+1) )
C        print*,k,z(k+1),u(ii,j,k+1),z(k+1) - z(k)
 101    CONTINUE
        DisplHist(i,icycle) = intgrl
C
C.....cubic spline (Euler- MacLaurin), linear at the boundary
C
c        xloc(1) = 0.
c        yloc(1) = 0.
c        DO 102 k = 2, 25
c          xloc(k) = 0.5 * (z(k) + z(k-1))
c          yloc(k) = u(ii,j,k)
c 102    CONTINUE
c        yp1 = 2 * u(ii,j,2) / (z(2) - z(1))
c        ypn = (u(ii,j,25)-u(ii,j,24))/(z(25)-z(24))
c        CALL SPLINE(xloc,yloc,25,yp1,ypn,y2)
C
c        intgrl2 = 0.
c        dzloc = 0.5 * (z(3) - z(2))
c        intgrl2 = intgrl2 + dzloc * 0.5 * (2 - yloc(1) - yloc(25))
c        DO 103 k  = 2, 47
c          xspl = dzloc * (k - 1)
c          CALL SPLINT(xloc,yloc,y2,25,xspl,yspl)
c          intgrl2 = intgrl2 + dzloc * (1 - yspl)
c 103    CONTINUE
c        intgrl2 = intgrl2 - dzloc*dzloc/12*(ypn-yp1)
C
C.....cubic spline (Euler- MacLaurin), slope from t_wall
C.....watch it: this routine can be Re-dependent!!!!!!!!!!!!
C.....use increased number of points
C
c        xloc(1) = 0.
c        yloc(1) = 0.
c        DO 104 k = 2, 25
c          xloc(k) = 0.5 * (z(k) + z(k-1))
c          yloc(k) = u(ii,j,k)
c 104    CONTINUE
cc        xx = 340. + Xpos(i)
cc        twall = 0.0225 * 0.37**(-0.25) * (xx / ru1)**(-0.2)
cc        write(6,*)'twall, predicted : ', twall
c        twall = twav(ii,j,kz1+1)
cc        write(6,*)'twall, computed  : ', twall
c        yp1 = (1/ru1) * twall
c        ypn = (u(ii,j,25)-u(ii,j,24))/(z(25)-z(24))
c        CALL SPLINE(xloc,yloc,25,yp1,ypn,y2)
C
c        intgrl3 = 0.
c        dzloc = 0.5 * (z(3) - z(2))
c        intgrl3 = intgrl3 + dzloc*0.5*(yloc(1)+yloc(25))
c        DO 105 k  = 2, 47
c          xspl = dzloc * (k - 1)
c          CALL SPLINT(xloc,yloc,y2,25,xspl,yspl)
c          intgrl3 = intgrl3 + dzloc * yspl
c 105    CONTINUE
c        intgrl3 = intgrl3 - dzloc*dzloc/12*(ypn-yp1)
c        intgrl3 = xloc(25) - intgrl3
C
C.....boundary layer thickness (99%)
C
c        zbpos1 = 2
c        zbpos2 = 30
c 106    CONTINUE
c        IF (u(ii,j,zbpos1).GE.0.985) GOTO 107
c        zbpos1 = zbpos1+1
c        GOTO 106 
c 107    CONTINUE
c        IF (u(ii,j,zbpos2) .LE. 0.995) GOTO 108
c        zbpos2 = zbpos2-1
c        GOTO 107
c 108    CONTINUE
c        IF (zbpos1.EQ.zbpos2) THEN
c          dzblloc = 1
c        ELSE
c          dzblloc = (0.99-u(ii,j,zbpos1)) /
c     &              (u(ii,j,zbpos2)-u(ii,j,zbpos1))
c        END IF
c        zblloc1  = (1-dzloc)*z(zbpos1) + dzloc*z(zbpos2)
c         zblloc1 = 10.
C
C.....b.l. thickness from spline
C
c        xloc(1) = 0.
c        yloc(1) = 0.
c        DO 109 k = 2, 25
c          xloc(k) = u(ii,j,k)
c          yloc(k) = 0.5 * (z(k) + z(k-1))
c 109    CONTINUE
c        yp1 = (z(2) - z(1)) / (2 * u(ii,j,2))
c        ypn = (z(25)-z(24)) / (u(ii,j,25)-u(ii,j,24))
c        CALL SPLINE(xloc,yloc,25,yp1,ypn,y2)
C
c        xspl = 0.99
c        CALL SPLINT(xloc,yloc,y2,25,xspl,yspl)
c        zblloc2 = yspl
C
C.....
C
c      WRITE(6,*) '*...............iter....... ',icycle
c      WRITE(6,*) '*..x position             = ',xpos(i)
c      WRITE(6,*) '*..d1(trap/spllin/spllog) = ',intgrl, intgrl2, intgrl3
c      WRITE(6,*) '*..d(simple/spl.)         = ',zblloc1, zblloc2
c      DisplHist(i,icycle) = intgrl3
 10   CONTINUE
C
c      WRITE(6,*) '*...................................'
      RETURN 
      END 
C
C 
C-----SUBROUTINE-CalcShape----------------------P. FLOHR--08/03/1994----
C
      SUBROUTINE CalcShape(u,x,z,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - some parameters to check the boundary layer
C                   development
C                                                       
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real u(nx,ny,nz), x(nx), z(nz)
      INTEGER   nx, ny, nz
C
      REAL      integral
      INTEGER   i, j, k, ii
C 
C-----------------------------------------------------------------------
C
C.....integrate computed velocity-field (trapeziod rule)
C
      j     = (jy2-jy1)/2
C
      DO 10 i = 1, 3
        ii = INT(ix2 * Xpos(i)/x(ix2))
        integral = 0.
        DO 101 k = 1, kz2-5
          integral = integral + (z(k+1) - z(k)) * 
     &                        (1 - (u(ii,j,k+1)/uinfty) )
 101    CONTINUE
        WRITE(6,*) '*..d1             = ',integral
 10   CONTINUE
C
      RETURN 
      END 
C
C
C
C
C 
C-----SUBROUTINE-CalcShape----------------------P. FLOHR--08/03/1994----
C
      SUBROUTINE CalcMoment(u,x,z,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - some parameters to check the boundary layer
C                   development
C                                                       
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real u(nx,ny,nz), x(nx), z(nz)
      INTEGER   nx, ny, nz
C
      REAL      integral,ucenter
      INTEGER   i, j, k, ii
C 
C-----------------------------------------------------------------------
C
C.....integrate computed velocity-field (trapeziod rule)
C
      j      = (jy2-jy1)/2
      ucenter= 1.5
C
      DO 10 i = 1, 3
        ii = INT(ix2 * Xpos(i)/x(ix2))
        integral = 0.
        DO 101 k = 2, kz2
        integral = integral + (z(k) - z(k-1))* 
     &               (0.25 - (u(ii,j,k)-ucenter)**2.)
 101    CONTINUE
        WRITE(6,*) '*...momentum thickness at x=',xpos(i),' d=',integral
 10   CONTINUE
C
      RETURN 
      END 







*//////////////////////////////////////////////////////////////////
*//                                                              //
      SUBROUTINE PRINTVAR(NU,NAME,V1,V2,V3,IM,JM,KM)
*//                                                              //
*//       Prints out any block of 3 variables                    //
*//                                                              //
*//////////////////////////////////////////////////////////////////

      real V1(IM,JM,KM),V2(IM,JM,KM),V3(IM,JM,KM)
      integer i,j,k,im,jm,km, nu
   
      CHARACTER*5 NAME
      open(NU,file=NAME,form='formatted',status='new')

      WRITE(NU,*) 'u-component'
      DO i=1,12
       WRITE(NU,*) i
       WRITE(NU,100) ((V1(i,j,k),j=1,12),k=12,1,-1)
      ENDDO

      WRITE(NU,*) 'v-component'
      DO i=1,12
       WRITE(NU,*) i
       WRITE(NU,100) ((V2(i,j,k),j=1,12),k=12,1,-1)
      ENDDO

      WRITE(NU,*) 'w-component'
      DO i=1,12
       WRITE(NU,*) i
       WRITE(NU,100) ((V3(i,j,k),j=1,12),k=12,1,-1)
      ENDDO

100   FORMAT( 12(2X,F8.2) )

       close (nu)
c
      STOP 
      END
C 
C-----SUBROUTINE-SCRVEL-------------------------P. FLOHR--01/02/1994----
C
      SUBROUTINE SCRVEL(u,v,w,p,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - display inst. velocity field
C                                                                       
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), p(nx,ny,nz)
      INTEGER nx, ny, nz
C
      CHARACTER*1 string(100)
      REAL f, h, sum, utmp
      REAL maxu, minu, maxv, minv, maxw, minw, maxp, minp
      INTEGER i, j, k, ii
C 
C-----------------------------------------------------------------------
C
C 
C.....check for maximum and minimum velocity
C
      maxu  = 0.
      minu  = 0.
      maxv  = 0.
      minv  = 0.
      maxw  = 0.
      minw  = 0.
      minp  = 0.
      maxp  = 0.
      DO i = 1, 100
        string(i) = ' '
      END DO
C
      DO 10 i = ix1+1, ix2
        DO 10 j = jy1+1, jy2
          DO 10 k = kz1+1, kz2
            IF (maxu .LT. u(i,j,k)) maxu = u(i,j,k)
            IF (minu .GT. u(i,j,k)) minu = u(i,j,k)
            IF (maxv .LT. v(i,j,k)) maxv = v(i,j,k)
            IF (minv .GT. v(i,j,k)) minv = v(i,j,k)
            IF (maxw .LT. w(i,j,k)) maxw = w(i,j,k)
            IF (minw .GT. w(i,j,k)) minw = w(i,j,k)
            IF (maxp .LT. p(i,j,k)) maxp = p(i,j,k)
            IF (minp .GT. p(i,j,k)) minp = p(i,j,k)
 10   CONTINUE
      IF (maxu .EQ. minu) maxu = 1
      IF (maxv .EQ. minv) maxv = 1
      IF (maxw .EQ. minw) maxw = 1
      IF (maxp .EQ. minp) maxp = 1
C
C.....write velocity field on the screen
C
      j = (jy2 - jy1)/2
      WRITE(6,*)'*.................................U..................'
      DO 20 k = kz2, kz1+1, -1
        ii = 0
        DO 201 i = ix1+1, ix2
          ii = ii + 1
          string(ii) = CHAR(48 
     &               + INT(9.9 * (u(i,j,k)-minu)/ (maxu-minu)))
c          string(ii) = CHAR(97 + INT(26. * u(i,j,k) / maxu))
 201    CONTINUE
        WRITE(6,1000) (string(i),i=1,100)
 20   CONTINUE
C
      WRITE(6,*)'*.................................V..................'
      DO 30 k = kz2, kz1+1, -1
        ii = 0
        DO 301 i = ix1+1, ix2
          ii = ii + 1
          string(ii) = CHAR(48 
     &               + INT(9.9 * (v(i,j,k)-minv) / (maxv-minv)))
c          string(ii) = CHAR(97 + INT(26. * v(i,j,k) / maxv))
 301    CONTINUE
        WRITE(6,1000) (string(i),i=1,100)
 30   CONTINUE
C
C
      WRITE(6,*)'*.................................W..................'
      DO 40 k = kz2, kz1+1, -1
        ii = 0
        DO 401 i = ix1+1, ix2
          ii = ii + 1
          string(ii) = CHAR(48
     &               + INT(9.9 * (w(i,j,k)-minw) / (maxw-minw)))
c          string(ii) = CHAR(97 + INT(26. * w(i,j,k) / maxw))
 401    CONTINUE
        WRITE(6,1000) (string(i),i=1,100)
 40   CONTINUE
C
      WRITE(6,*)'*.................................P..................'
      DO 50 k = kz2, kz1+1, -1
        ii = 0
        DO 501 i = ix1+1, ix2
          ii = ii + 1
          string(ii) = CHAR(48
     &               + INT(9.9 * (p(i,j,k)-minp) / (maxp-minp)))
c          string(ii) = CHAR(97 + INT(26. * p(i,j,k) / maxp))
 501    CONTINUE
        WRITE(6,1000) (string(i),i=1,100)
 50   CONTINUE

 1000 FORMAT (1X,100(A1))
C
      RETURN 
      END 
C 
C-----SUBROUTINE-SCRRHS-------------------------P. FLOHR--01/02/1994----
C
      SUBROUTINE SCRRHS(rhsu,rhsv,rhsw,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - display RHS
C                                                                       
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real rhsu(nx,ny,nz), rhsv(nx,ny,nz), rhsw(nx,ny,nz)
      INTEGER nx, ny, nz
C
      REAL maxu, minu, maxv, minv, maxw, minw
      INTEGER i, j, k
      INTEGER iumin,iumax,ivmin,ivmax,iwmin,iwmax,
     &        jumin,jumax,jvmin,jvmax,jwmin,jwmax,
     &        kumin,kumax,kvmin,kvmax,kwmin,kwmax

C 
C-----------------------------------------------------------------------
C
C 
C.....check for maximum and minimum values
C
      maxu  = 0.
      minu  = 0.
      maxv  = 0.
      minv  = 0.
      maxw  = 0.
      minw  = 0.
C
      DO 10 i = ix1+1, ix2
        DO 10 j = jy1+1, jy2
          DO 10 k = kz1+1, kz2
            IF (maxu .LT. rhsu(i,j,k)) THEN
              maxu  = rhsu(i,j,k)
              iumax = i
              jumax = j
              kumax = k
            END IF
            IF (minu .GT. rhsu(i,j,k)) THEN
              minu = rhsu(i,j,k)
              iumin = i
              jumin = j
              kumin = k
            END IF
            IF (maxv .LT. rhsv(i,j,k)) THEN
              maxv = rhsv(i,j,k)
              ivmax = i
              jvmax = j
              kvmax = k
            END IF
            IF (minv .GT. rhsv(i,j,k)) THEN
              minv = rhsv(i,j,k)
              ivmin = i
              jvmin = j
              kvmin = k
            END IF
            IF (maxw .LT. rhsw(i,j,k)) THEN
              maxw = rhsw(i,j,k)
              iwmax = i
              jwmax = j
              kwmax = k
            END IF
            IF (minw .GT. rhsw(i,j,k)) THEN
              minw = rhsw(i,j,k)
              iwmin = i
              jwmin = j
              kwmin = k
            END IF
 10   CONTINUE
      WRITE(6,*)'*..RHS U min      = ',minu, iumin,jumin,kumin
      WRITE(6,*)'*..RHS U max      = ',maxu, iumax,jumax,kumax
      WRITE(6,*)'*..RHS V min      = ',minv, ivmin,jvmin,kvmin
      WRITE(6,*)'*..RHS V max      = ',maxv, ivmax,jvmax,kvmax
      WRITE(6,*)'*..RHS W min      = ',minw, iwmin,jwmin,kwmin
      WRITE(6,*)'*..RHS W max      = ',maxw, iwmax,jwmax,kwmax
C
      RETURN 
      END 
C 
C-----SUBROUTINE-SCRUVWP------------------------P. FLOHR--01/02/1994----
C
      SUBROUTINE SCRUVWP(u,v,w,p,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - CHECK FOR MAX. AND MIN. VELOCITY IN THE FLOW FIELD
C                                                       
C-----------------------------------------------------------------------
C
      INCLUDE "mpif.h"
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),
     &          p(nx,ny,nz)
      INTEGER   nx, ny, nz
C
      REAL      f, h, sum, utmp,in(2),out(2)
      REAL      maxu, minu, maxv, minv, maxw, minw, minp, maxp      
      INTEGER   i, j, k,mytype1,mytype2,indu(3),indv(3),indw(3),indp(3),
     &                                  induM(3),indvM(3),indwM(3),indpM(3)
      integer mrm1, kminv, kmaxv, kmaxu, kminp, kmaxw, kminw, kminu, kmaxp,
     &        jminv, jminw, jminu, jmaxv, jmaxw, jmaxu, jminp, jmaxp,
     &        iminv, iminw, iminu, imaxv, imaxw, imaxu, iminp, imaxp
      
      integer iim, jjm, kkm, im, jm, km, indx, ier
     
#ifdef SINGLE
      mytype1=MPI_2REAL
      mytype2=MPI_REAL
#endif
#ifdef DOUBLE
      mytype1=MPI_2DOUBLE_PRECISION
      mytype2=MPI_DOUBLE_PRECISION
#endif

      
C 
C-----------------------------------------------------------------------
C 
C.....check for maximum and minimum velocity
C*******************MIN********************
C
      im=nx
      jm=ny
      km=nz
      minu=0.
      maxu=0.
      minv=0.
      maxv=0.
      minw=0.
      maxw=0.
      minp=0.
      maxp=0.
      do k=1+(mysize-1+myrank)/mysize,nz-1+(myrank+1)/mysize
       do j=1,ny
        do i=1,nx
	 if (u(i,j,k).lt.minu) then
	  minu=u(i,j,k)
	  indu(1)=i
	  indu(2)=j
	  indu(3)=k
	 endif 
	 if (v(i,j,k).lt.minv) then
	  minv=v(i,j,k)
	  indv(1)=i
	  indv(2)=j
	  indv(3)=k
	 endif 
	 if (w(i,j,k).lt.minw) then
	  minw=w(i,j,k)
	  indw(1)=i
	  indw(2)=j
	  indw(3)=k
	 endif 
	 if (p(i,j,k).lt.minp) then
	  minp=p(i,j,k)
	  indp(1)=i
	  indp(2)=j
	  indp(3)=k
	 endif 
	 if (u(i,j,k).gt.maxu) then
	  maxu=u(i,j,k)
	  induM(1)=i
	  induM(2)=j
	  induM(3)=k
	 endif 
	 if (v(i,j,k).gt.maxv) then
	  maxv=v(i,j,k)
	  indvM(1)=i
	  indvM(2)=j
	  indvM(3)=k
	 endif 
	 if (w(i,j,k).gt.maxw) then
          maxw=w(i,j,k)
	  indwM(1)=i
	  indwM(2)=j
	  indwM(3)=k
	 endif 
	 if (p(i,j,k).gt.maxp) then
          maxp=p(i,j,k)
	  indpM(1)=i
	  indpM(2)=j
	  indpM(3)=k
	 endif 
        enddo
       enddo
      enddo      		 

      iim=indu(1)
      jjm=indu(2)
      kkm=indu(3)

      in(1)=minu
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)

      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MINLOC,0,MPI_COMM_HYDRO,ier)
      
      minu=out(1)
      indx=out(2)

      mrm1=indx/(im*jm*km)               !myrank
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm) 
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jminu=jjm+1
      iminu=iim+1
      kminu=kkm+1+mrm1*(km-2)

c***************************MAX      

      iim=induM(1)
      jjm=induM(2)
      kkm=induM(3)

      in(1)=maxu
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MAXLOC,0,MPI_COMM_HYDRO,ier)
      
      maxu=out(1)
      indx=out(2)
      
      mrm1=indx/(im*jm*km)
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm)
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jmaxu=jjm+1
      imaxu=iim+1
      kmaxu=kkm+1+mrm1*(km-2)
      
      

      IF(MYRANK.EQ.0) then
      
         WRITE(*,100) '*..u=', minu, maxu,'(',iminu,jminu,kminu,')',
     &        '(',imaxu,jmaxu,kmaxu,')'
      ENDIF

C.....check for maximum and minimum velocity
C*******************MIN********************
C
      iim=indv(1)
      jjm=indv(2)
      kkm=indv(3)

      in(1)=minv
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MINLOC,0,MPI_COMM_HYDRO,ier)
      
      minv=out(1)
      indx=out(2)
      
      mrm1=indx/(im*jm*km)
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm)
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jminv=jjm+1
      iminv=iim+1
      kminv=kkm+1+mrm1*(km-2)

c***************************MAX      
      iim=indvM(1)
      jjm=indvM(2)
      kkm=indvM(3)

      in(1)=maxv
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MAXLOC,0,MPI_COMM_HYDRO,ier)
      
      maxv=out(1)
      indx=out(2)
      
      mrm1=indx/(im*jm*km)
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm)
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jmaxv=jjm+1
      imaxv=iim+1
      kmaxv=kkm+1+mrm1*(km-2)
      
      

      IF(MYRANK.EQ.0) then
      
         WRITE(*,100) '*..v=', minv, maxv,'(',iminv,jminv,kminv,')',
     &        '(',imaxv,jmaxv,kmaxv,')'
      ENDIF

C.....check for maximum and minimum velocity
C*******************MIN********************
C
      iim=indw(1)
      jjm=indw(2)
      kkm=indw(3)

      in(1)=minw
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MINLOC,0,MPI_COMM_HYDRO,ier)
      
      minw=out(1)
      indx=out(2)
      
      mrm1=indx/(im*jm*km)
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm)
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jminw=jjm+1
      iminw=iim+1
      kminw=kkm+1+mrm1*(km-2)

c***************************MAX      
      iim=indwM(1)
      jjm=indwM(2)
      kkm=indwM(3)

      in(1)=maxw
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MAXLOC,0,MPI_COMM_HYDRO,ier)
      
      maxw=out(1)
      indx=out(2)
      
      mrm1=indx/(im*jm*km)
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm)
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jmaxw=jjm+1
      imaxw=iim+1
      kmaxw=kkm+1+mrm1*(km-2)
      
      

      IF(MYRANK.EQ.0) then
      
         WRITE(*,100) '*..w=', minw, maxw,'(',iminw,jminw,kminw,')',
     &        '(',imaxw,jmaxw,kmaxw,')'
      ENDIF

C.....check for maximum and minimum velocity
C*******************MIN********************
C
      iim=indp(1)
      jjm=indp(2)
      kkm=indp(3)

      in(1)=minp
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MINLOC,0,MPI_COMM_HYDRO,ier)
      
      minp=out(1)
      indx=out(2)
      
      mrm1=indx/(im*jm*km)
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm)
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jminp=jjm+1
      iminp=iim+1
      kminp=kkm+1+mrm1*(km-2)

c***************************MAX      
      iim=indpM(1)
      jjm=indpM(2)
      kkm=indpM(3)

      in(1)=maxp
      in(2)=iim-1+im*(jjm-1)+im*jm*(kkm-1)+im*jm*km*(myrank)
      CALL MPI_REDUCE(in,out,1,mytype1,MPI_MAXLOC,0,MPI_COMM_HYDRO,ier)
      
      maxp=out(1)
      indx=out(2)
      
      mrm1=indx/(im*jm*km)
      kkm=(indx-(im*jm*km)*mrm1)/(im*jm)
      jjm=(indx-(im*jm*km)*mrm1-kkm*jm*im)/im
      iim=(indx-(im*jm*km)*mrm1-kkm*jm*im-jjm*im)
      
      jmaxp=jjm+1
      imaxp=iim+1
      kmaxp=kkm+1+mrm1*(km-2)
      
      

      IF(MYRANK.EQ.0) then
      
         WRITE(*,100) '*..p=', minp, maxp,'(',iminp,jminp,kminp,')',
     &        '(',imaxp,jmaxp,kmaxp,')'
      ENDIF



C
 100  FORMAT(1X, A, 1p2e10.3, A, 3I5, 2A, 3I5, A)
      RETURN 
      END 
