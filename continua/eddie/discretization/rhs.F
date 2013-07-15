C
C 
C-----SUBROUTINE-Rhs---------------------------E. Balaras  4/12/98------
C
      SUBROUTINE Rhs(uo,vo,wo,upr,vpr,wpr,tv,wfrg,ru,rv,rw,nx,ny,nz )
C
C      Driver routine
C      Calls different RHS routines
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
      float uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),upr(nx,ny,nz),
     &      vpr(nx,ny,nz),wpr(nx,ny,nz),tv(nx,ny,nz),wfrg(nx,ny,nz),
     &      ru(nx,ny,nz),rv(nx,ny,nz),rw(nx,ny,nz)
      int   nx,ny,nz
C 
C-----------------------------------------------------------------------
C
c...rhs in conservative formulation
      CALL Rhscent(uo,vo,wo,upr,vpr,wpr,tv,ru,rv,rw,nx,ny,nz)

c...modify rhs for fringe method
c      IF (Lfrg) CALL Rhsfrg(uo,vo,wo,ru,rv,rw,wfrg,nx,ny,nz)
c
      END subroutine Rhs
C
C
C-----SUBROUTINE-Rhscent----------------------E. Balaras  7/1/99------
C 
      SUBROUTINE Rhscent(UNI,VNI,WNI,upr,vpr,wpr,
     &                   TV,RU,RV,RW,IMAX,JMAX,KMAX)
C
C      CREATES THE RIGHT SIDE OF THE MOM. EQ.                         
C      (ADVECTION AND VISCOUS TERMS)                                  
C         EXPLICIT VERSION                                         
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      integer imax, jmax, kmax
      REAL   UNI(IMAX,JMAX,KMAX),VNI(IMAX,JMAX,KMAX),
     &       WNI(IMAX,JMAX,KMAX),TV(IMAX,JMAX,KMAX),
     &       RU(IMAX,JMAX,KMAX),RV(IMAX,JMAX,KMAX),
     &       RW(IMAX,JMAX,KMAX)
      REAL   upr( IMAX,JMAX,KMAX),vpr( IMAX,JMAX,KMAX),
     &       wpr( IMAX,JMAX,KMAX)
      REAL   rflag, dx1, dy1, dz1
      
      integer i,j,k
      real wzminus, wzplus, wyminus, wyplus, wxminus, wxplus,
     &     vzminus, vzplus, vyminus, vyplus, vxminus, vxplus,
     &     uzminus, uzplus, uyminus, uyplus, uxminus, uxplus,
     &     txzp, txzm, tzzp, tzzm, tyzm, tyyp, tyzp, tvip, txyp,
     &     tvkp, txxm, txym, txxp, tyym, tvjp, tvim, tvjm, tvkm,
     &     dwdzm, dwdzp, dwdxp, dwdym, dvdzp, dvdym, dvdxp,
     &     dudym, dvdxm, dudzp, dudzm, dvdyp, dvdzm, dwdyp,
     &     dudxp, dudxm, dudyp, dwdxm
      
c
      rflag=real(impfla)
C
CC
CCCCCCCCCCCCCCCCCCC     FOR THE U COMPONENT     CCCCCCCCCCCCCCCCCC
CC
C
c..get coefficients
       DX1=1.0/DX
       DY1=1.0/DY
       DZ1=1.0/DZ
c
      DO 10 K=KZ1+1,KZ2
       DO 20 J=JY1+1,JY2
        DO 30 I=IBU,IEU
c...get velocities at 1/2 locations
         UXPLUS =(UNI(I+1,J  ,K  )+UNI(I  ,J  ,K  ))*0.5
         UXMINUS=(UNI(I  ,J  ,K  )+UNI(I-1,J  ,K  ))*0.5
c
         VXPLUS =(VNI(I+1,J  ,K  )+VNI(I  ,J  ,K  ))*0.5
         VXMINUS=(VNI(I+1,J-1,K  )+VNI(I  ,J-1,K  ))*0.5
c
         WXPLUS =(WNI(I+1,J  ,K  )+WNI(I  ,J  ,K  ))*0.5
         WXMINUS=(WNI(I+1,J  ,K-1)+WNI(I  ,J  ,K-1))*0.5
c
         UYPLUS =(UNI(I  ,J+1,K  )+UNI(I  ,J  ,K  ))*0.5
         UYMINUS=(UNI(I  ,J  ,K  )+UNI(I  ,J-1,K  ))*0.5
c
         UZPLUS =(UNI(I  ,J  ,K+1)+UNI(I  ,J  ,K  ))*0.5
         UZMINUS=(UNI(I  ,J  ,K  )+UNI(I  ,J  ,K-1))*0.5

c...get derivatives at 1/2 locations
         dudxp= ap(i+1)*(uni(i+1,j,k)-uni(i  ,j,k))*dx1
         dudxm= ap(i  )*(uni(i  ,j,k)-uni(i-1,j,k))*dx1
         dudyp= bv(j  )*(uni(i,j+1,k)-uni(i,j  ,k))*dy1
         dudym= bv(j-1)*(uni(i,j  ,k)-uni(i,j-1,k))*dy1
         dudzp= cw(k  )*(uni(i,j,k+1)-uni(i,j,k  ))*dz1
         dudzm= cw(k-1)*(uni(i,j,k  )-uni(i,j,k-1))*dz1
         dvdxp= au(i)*(vni(i+1,j  ,k)-vni(i,j  ,k))*dx1
         dvdxm= au(i)*(vni(i+1,j-1,k)-vni(i,j-1,k))*dx1
         dwdxp= au(i)*(wni(i+1,j,k  )-wni(i,j,k  ))*dx1
         dwdxm= au(i)*(wni(i+1,j,k-1)-wni(i,j,k-1))*dx1
c...get nu_t where needed
         tvjp=0.25*(tv(i,j  ,k)+tv(i+1,j  ,k)
     %             +tv(i,j+1,k)+tv(i+1,j+1,k))
         tvjm=0.25*(tv(i,j  ,k)+tv(i+1,j  ,k)
     %             +tv(i,j-1,k)+tv(i+1,j-1,k))
         tvkp=0.25*(tv(i,j,k  )+tv(i+1,j,k  )
     %             +tv(i,j,k+1)+tv(i+1,j,k+1))
         tvkm=0.25*(tv(i,j,k  )+tv(i+1,j,k  )
     %             +tv(i,j,k-1)+tv(i+1,j,k-1))
c...flux of normal total stresses
         txxp=(ru1+2.*tv(i+1,j,k))*dudxp
         txxm=(ru1+2.*tv(i  ,j,k))*dudxm
         tyyp=(ru1+tvjp)*dudyp
         tyym=(ru1+tvjm)*dudym
         tzzp=(ru1+tvkp)*dudzp
         tzzm=(ru1+tvkm)*dudzm
c...flux of cross sgs stresses
         txyp=tvjp*dvdxp
         txym=tvjm*dvdxm
         txzp=tvkp*dwdxp
         txzm=tvkm*dwdxm
c
c......calculate RHS for u-momentum
c
         RU(i,j,k)=
c..advective term in conservative formulation
     %   -au(i)*(uxplus*uxplus-uxminus*uxminus)*DX1
     %   -bu(j)*(vxplus*uyplus-vxminus*uyminus)*DY1
     %   -cu(k)*(wxplus*uzplus-wxminus*uzminus)*DZ1
c..viscous+part of sgs diffusion
     %   +au(i)*(txxp-txxm)*dx1*(1.-rflag)
     %   +bu(j)*(tyyp-tyym)*dy1
     %   +cu(k)*(tzzp-tzzm)*dz1
c..rest of sgs diffusion
     %   +bu(j)*(txyp-txym)*dy1
     %   +cu(k)*(txzp-txzm)*dz1

c..store explicit part of wall normal diffusio for C.N.
         upr(i,j,k) = au(i)*(txxp-txxm)*dx1*rflag
c


30       continue     
20      continue
10     continue
C
CC
CCCCCCCCCCCCCCCCCCC     FOR THE V COMPONENT     CCCCCCCCCCCCCCCCCC
CC
C
      do 40 k=kz1+1,kz2
       do 50 J=JBV,JEV
        do 60 i=ix1+1,ix2
c...get velocities at 1/2 locations
         vxplus= (VNI(i+1,j  ,k  )+VNI(i  ,j  ,k  ))*0.5
         vxminus=(VNI(i  ,j  ,k  )+VNI(i-1,j  ,k  ))*0.5
c
         vyplus= (VNI(i  ,j+1,k  )+VNI(i  ,j  ,k  ))*0.5
         vyminus=(VNI(i  ,j  ,k  )+VNI(i  ,j-1,k  ))*0.5
c
         vzplus= (VNI(i  ,j  ,k+1)+VNI(i  ,j  ,k  ))*0.5
         vzminus=(VNI(i  ,j  ,k  )+VNI(i  ,j  ,k-1))*0.5
c
         uyplus= (UNI(i  ,j+1,k  )+UNI(i  ,j  ,k  ))*0.5
         uyminus=(UNI(i-1,j+1,k  )+UNI(i-1,j  ,k  ))*0.5
c
         wyplus= (WNI(i  ,j+1,k  )+WNI(i  ,j  ,k  ))*0.5   
         wyminus=(WNI(i  ,j+1,k-1)+WNI(i  ,j  ,k-1))*0.5

c...get derivatives at 1/2 locations
         dvdxp= au(i  )*(vni(i+1,j,k)-vni(i  ,j,k))*dx1
         dvdxm= au(i-1)*(vni(i  ,j,k)-vni(i-1,j,k))*dx1
         dvdyp= bp(j+1)*(vni(i,j+1,k)-vni(i,j  ,k))*dy1
         dvdym= bp(j  )*(vni(i,j  ,k)-vni(i,j-1,k))*dy1
         dvdzp= cw(k  )*(vni(i,j,k+1)-vni(i,j,  k))*dz1
         dvdzm= cw(k-1)*(vni(i,j,k  )-vni(i,j,k-1))*dz1
         dudyp= bv(j)*(uni(i  ,j+1,k)-uni(i  ,j,k))*dy1
         dudym= bv(j)*(uni(i-1,j+1,k)-uni(i-1,j,k))*dy1
         dwdyp= bv(j)*(wni(i,j+1,k  )-wni(i,j,k  ))*dy1
         dwdym= bv(j)*(wni(i,j+1,k-1)-wni(i,j,k-1))*dy1
c...get nu_t where needed
         tvip=0.25*(tv(i,j  ,k)+tv(i+1,j  ,k)
     %             +tv(i,j+1,k)+tv(i+1,j+1,k))
         tvim=0.25*(tv(i,j  ,k)+tv(i-1,j  ,k)
     %             +tv(i,j+1,k)+tv(i-1,j+1,k))
         tvkp=0.25*(tv(i,j  ,k)+tv(i,j  ,k+1)
     %             +tv(i,j+1,k)+tv(i,j+1,k+1))
         tvkm=0.25*(tv(i,j  ,k)+tv(i,j  ,k-1)
     %             +tv(i,j+1,k)+tv(i,j+1,k-1))
c...flux of normal total stresses
         txxp=(ru1+tvip)*dvdxp
         txxm=(ru1+tvim)*dvdxm
         tyyp=(ru1+2.*tv(i,j+1,k))*dvdyp
         tyym=(ru1+2.*tv(i,j  ,k))*dvdym
         tzzp=(ru1+tvkp)*dvdzp
         tzzm=(ru1+tvkm)*dvdzm
c...flux of cross sgs stresses
         txyp=tvip*dudyp
         txym=tvim*dudym
         tyzp=tvkp*dwdyp
         tyzm=tvkm*dwdym
c
c......calculate RHS for v-momentum
c
           RV(i,j,k)=
c..advective term in conservative formulation
     %      -av(i)*(uyplus*vxplus-uyminus*vxminus)*DX1
     %      -bv(j)*(vyplus*vyplus-vyminus*vyminus)*DY1
     %      -cv(k)*(wyplus*vzplus-wyminus*vzminus)*DZ1
c..viscous+part of sgs diffusion
     %      +av(i)*(txxp-txxm)*dx1*(1.-rflag)
     %      +bv(j)*(tyyp-tyym)*dy1
     %      +cv(k)*(tzzp-tzzm)*dz1
c..rest of sgs diffusion
     %      +av(i)*(txyp-txym)*dx1
     %      +cv(k)*(tyzp-tyzm)*dz1

c..store explicit part of wall normal diffusio for C.N.
           vpr(i,j,k) = av(i)*(txxp-txxm)*dx1*rflag
c
 60      continue
 50     continue
 40    continue
C
CC
CCCCCCCCCCCCCCCCCCC     FOR THE W COMPONENT     CCCCCCCCCCCCCCCCCC
CC
C
      do 70 K=KBW,KEW
       do 80 j=jy1+1,jy2
        do 90 i=ix1+1,ix2
c...get velcities at 1/2 locations
          wxplus =(WNI(i+1,j  ,k  )+WNI(i  ,j  ,k  ))*0.5
          wxminus=(WNI(i  ,j  ,k  )+WNI(i-1,j  ,k  ))*0.5
c
          wyplus =(WNI(i  ,j+1,k  )+WNI(i  ,j  ,k  ))*0.5   
          wyminus=(WNI(i  ,j  ,k  )+WNI(i  ,j-1,k  ))*0.5
c
          wzplus =(WNI(i  ,j  ,k+1)+WNI(i  ,j  ,k  ))*0.5
          wzminus=(WNI(i  ,j  ,k  )+WNI(i  ,j  ,k-1))*0.5
c
          uzplus =(UNI(i  ,j  ,k+1)+UNI(i  ,j  ,k  ))*0.5
          uzminus=(UNI(i-1,j  ,k+1)+UNI(i-1,j  ,k  ))*0.5
c
          vzplus =(VNI(i  ,j  ,k+1)+VNI(i  ,j  ,k  ))*0.5
          vzminus=(VNI(i  ,j-1,k+1)+VNI(i  ,j-1,k  ))*0.5

c...get derivatives at 1/2 locations
          dwdxp= au(i  )*(wni(i+1,j,k)-wni(i  ,j,k))*dx1
          dwdxm= au(i-1)*(wni(i  ,j,k)-wni(i-1,j,k))*dx1
          dwdyp= bv(j  )*(wni(i,j+1,k)-wni(i,j  ,k))*dy1
          dwdym= bv(j-1)*(wni(i,j  ,k)-wni(i,j-1,k))*dy1
          dwdzp= cp(k+1)*(wni(i,j,k+1)-wni(i,j,k  ))*dz1
          dwdzm= cp(k  )*(wni(i,j,k  )-wni(i,j,k-1))*dz1
          dudzp= cw(k)*(uni(i  ,j,k+1)-uni(i  ,j,k))*dz1
          dudzm= cw(k)*(uni(i-1,j,k+1)-uni(i-1,j,k))*dz1
          dvdzp= cw(k)*(vni(i,j  ,k+1)-vni(i,j  ,k))*dz1
          dvdzm= cw(k)*(vni(i,j-1,k+1)-vni(i,j-1,k))*dz1
c...get nu_t where needed
         tvip=0.25*(tv(i,j  ,k)+tv(i+1,j  ,k)
     %             +tv(i,j,k+1)+tv(i+1,j,k+1))
         tvim=0.25*(tv(i,j  ,k)+tv(i-1,j  ,k)
     %             +tv(i,j,k+1)+tv(i-1,j,k+1))
         tvjp=0.25*(tv(i,j  ,k)+tv(i,j+1,k  )
     %             +tv(i,j,k+1)+tv(i,j+1,k+1))
         tvjm=0.25*(tv(i,j  ,k)+tv(i,j-1,k  )
     %             +tv(i,j,k+1)+tv(i,j-1,k+1))
c...flux of normal total stresses
         txxp=(ru1+tvip)*dwdxp
         txxm=(ru1+tvim)*dwdxm
         tyyp=(ru1+tvjp)*dwdyp
         tyym=(ru1+tvjm)*dwdym
         tzzp=(ru1+2.*tv(i,j,k+1))*dwdzp
         tzzm=(ru1+2.*tv(i,j,k  ))*dwdzm
c...flux of cross sgs stresses
         txzp=tvip*dudzp
         txzm=tvim*dudzm
         tyzp=tvjp*dvdzp
         tyzm=tvjm*dvdzm
c
c......calculate RHS for w-momentum
c
        RW(i,j,k)=
c..advective term in conservative formulation
     %   -aw(i)*(uzplus * wxplus - uzminus * wxminus)*dx1
     %   -bw(j)*(vzplus * wyplus - vzminus * wyminus)*dy1
     %   -cw(k)*(wzplus * wzplus - wzminus * wzminus)*dz1
c..viscous+part of sgs diffusion
     %   +aw(i)*(txxp-txxm)*dx1*(1.-rflag)
     %   +bw(j)*(tyyp-tyym)*dy1
     %   +cw(k)*(tzzp-tzzm)*dz1
c..rest of sgs diffusion
     %   +aw(i)*(txzp-txzm)*dx1
     %   +bw(j)*(tyzp-tyzm)*dy1

c..store explicit part of wall normal diffusio for C.N.
           wpr(i,j,k) = aw(i)*(txxp-txxm)*dx1*rflag
c
 90     continue
 80    continue
 70   continue
c
      return
      end

