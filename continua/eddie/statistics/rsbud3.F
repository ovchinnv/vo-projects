      module rsbud
c     reynolds stress budget computation in the boundary layer
c     v.3 -- same as v.1 but with subgrid scale terms
c    DEBUG VERSION





!ccccccccVARIABLES FOR REYNOLDS STRESS BUDGETS -- SEE RSBUD1.F routine
!        NOTE:  TO SAVE MEMORY rsblk_1 overwrites uab2 (1st 4*7 nx*nz planes) This works if ny>=28!!!
!                if LES is on rsblk_3 overwrites vab2  (1st 4*4 nx*nz planes)
      real rsblk_1(nx,nz,ny)
      real rsblk_2(nx,0:nz+1,33)
      real rsblkg_1(nx,nzg,7) 
#ifdef LES
      real rsblk_3 (nx,nz,ny)     !temporary LES block
      real rsblk_4 (nx,0:nz+1,25) !permanent LES block
      real rsblkg_3(nx,nzg,4)     !temporary global LES block
#endif
!

!
      character*4 rsnames(0:3)
      
      
      
      


      
      
      subroutine rsbud(uu_conv,uu_prod,uu_diss,uu_ttrans,uu_ptrans,uu_pstr,uu_vdiff,
     &     vv_conv,vv_prod,vv_diss,vv_ttrans,vv_ptrans,vv_pstr,vv_vdiff,
     &     ww_conv,ww_prod,ww_diss,ww_ttrans,ww_ptrans,ww_pstr,ww_vdiff,
     &     uw_conv,uw_prod,uw_diss,uw_ttrans,uw_ptrans,uw_pstr,uw_vdiff,
#ifdef LES
     &     uu_prod_sgs,uu_diss_sgs,uu_tran_sgs,
     &     vv_prod_sgs,vv_diss_sgs,vv_tran_sgs,
     &     ww_prod_sgs,ww_diss_sgs,ww_tran_sgs,
     &     uw_prod_sgs,uw_diss_sgs,uw_tran_sgs, !12
#endif                     
     &     u_,v_,w_,p_,uu_,uw_,pu_,uuu_,uuw_,pdudx_,dudx2_,dudy2_,dudz2_,
     &     vv_,vvu_,vvw_,pdvdy_,dvdx2_,dvdy2_,dvdz2_,ww_,wwu_,www_,pw_,
     &     pdwdz_,dwdx2_,dwdy2_,dwdz2_,duwdxx_,duwdyy_,duwdzz_,pdudz_,pdwdx_,
#ifdef LES
     &     nududx_,nuduwdzx_,nududx2_,nuduwdzxdudz_,nuduvdyxdudy_,nudu2dx_,nuduwdzxu_,
     &	   nudvdy2_,nuduvdyxdvdx_,nudvwdzydvdz_,nudvudxyv_,nudvwdzyv_,nudwdz_,nudwdz2_,
     &     nuduwdzxdwdx_,nudvwdzydwdy_,nudw2dz_,nuduwdzxw_,nududxdwdx_,nudwdzdudz_,
     &     nuduwdzxduwdxz_,nuduvdyxdwdy_,nudwvdyzdudy_,nududxw_,nudwdzu_, !25
#endif
     &     uo,vo,wo,po,
#ifdef LES
     &     tv,
#endif          
     &     samples)
c	   
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
c     variable declarations
      integer samples,i,j,k,kbg,knd
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),po(nx,ny,nz)
      real fac
c      
      real uu_conv(nx,nz),uu_prod(nx,nz),uu_diss(nx,nz),uu_ttrans(nx,nz),uu_ptrans(nx,nz),uu_pstr(nx,nz),uu_vdiff(nx,nz),
     &     vv_conv(nx,nz),vv_prod(nx,nz),vv_diss(nx,nz),vv_ttrans(nx,nz),vv_ptrans(nx,nz),vv_pstr(nx,nz),vv_vdiff(nx,nz),
     &     ww_conv(nx,nz),ww_prod(nx,nz),ww_diss(nx,nz),ww_ttrans(nx,nz),ww_ptrans(nx,nz),ww_pstr(nx,nz),ww_vdiff(nx,nz),
     &     uw_conv(nx,nz),uw_prod(nx,nz),uw_diss(nx,nz),uw_ttrans(nx,nz),uw_ptrans(nx,nz),uw_pstr(nx,nz),uw_vdiff(nx,nz)
#ifdef LES
      real uu_prod_sgs(nx,nz),vv_prod_sgs(nx,nz),ww_prod_sgs(nx,nz),uw_prod_sgs(nx,nz), 
     &     uu_diss_sgs(nx,nz),vv_diss_sgs(nx,nz),ww_diss_sgs(nx,nz),uw_diss_sgs(nx,nz), 
     &     uu_tran_sgs(nx,nz),vv_tran_sgs(nx,nz),ww_tran_sgs(nx,nz),uw_tran_sgs(nx,nz),
     &     tv(nx,ny,nz) 
#endif                     
      real u_(nx,0:nz+1),v_(nx,0:nz+1),w_(nx,0:nz+1),p_(nx,0:nz+1),uu_(nx,0:nz+1),uw_(nx,0:nz+1)
     &     ,pu_(nx,0:nz+1),uuu_(nx,0:nz+1),uuw_(nx,0:nz+1),pdudx_(nx,0:nz+1),dudx2_(nx,0:nz+1)
     &     ,dudy2_(nx,0:nz+1),dudz2_(nx,0:nz+1),vv_(nx,0:nz+1),vvu_(nx,0:nz+1),vvw_(nx,0:nz+1)
     &     ,pdvdy_(nx,0:nz+1),dvdx2_(nx,0:nz+1),dvdy2_(nx,0:nz+1),dvdz2_(nx,0:nz+1),ww_(nx,0:nz+1)
     &     ,wwu_(nx,0:nz+1),www_(nx,0:nz+1),pw_(nx,0:nz+1),pdwdz_(nx,0:nz+1),dwdx2_(nx,0:nz+1)
     &     ,dwdy2_(nx,0:nz+1),dwdz2_(nx,0:nz+1),duwdxx_(nx,0:nz+1),duwdyy_(nx,0:nz+1),duwdzz_(nx,0:nz+1)
     &     ,pdudz_(nx,0:nz+1),pdwdx_(nx,0:nz+1)
#ifdef LES
      real nududx_(nx,0:nz+1),nuduwdzx_(nx,0:nz+1),nududx2_(nx,0:nz+1),nuduwdzxdudz_(nx,0:nz+1)
     &    ,nuduvdyxdudy_(nx,0:nz+1),nudu2dx_(nx,0:nz+1),nuduwdzxu_(nx,0:nz+1),nudvdy2_(nx,0:nz+1)
     &    ,nuduvdyxdvdx_(nx,0:nz+1),nudvwdzydvdz_(nx,0:nz+1),nudvudxyv_(nx,0:nz+1),nudvwdzyv_(nx,0:nz+1)
     &    ,nudwdz_(nx,0:nz+1)
     &    ,nudwdz2_(nx,0:nz+1),nuduwdzxdwdx_(nx,0:nz+1),nudvwdzydwdy_(nx,0:nz+1),nudw2dz_(nx,0:nz+1)
     &    ,nuduwdzxw_(nx,0:nz+1),nududxdwdx_(nx,0:nz+1),nudwdzdudz_(nx,0:nz+1),nuduwdzxduwdxz_(nx,0:nz+1)
     &    ,nuduvdyxdwdy_(nx,0:nz+1),nudwvdyzdudy_(nx,0:nz+1),nududxw_(nx,0:nz+1),nudwdzu_(nx,0:nz+1) !25
#endif
c     
      real u,v,w,p,uu,uw,pu,uuu,uuw,pdudx,dudx2,dudy2,dudz2,vv,vvu,
     &     vvw,pdvdy,dvdx2,dvdy2,dvdz2,ww,wwu,www,pw,pdwdz,dwdx2,
     &     dwdy2,dwdz2,duwdxx,duwdyy,duwdzz,pdudz,pdwdx
      real dvdz,dvdx,dvdy,dudx,dwdz,dudy,dudz,dwdx,dwdy,
     &     uutr,uutl,uur,uul,uubr,uubl,ucent,ucen,ucenr,ucenl,ucenb,
     &     wcen,wcent,wcenb,vcen,vcenr,vcenl,vcent,vcenb,wcenr,wcenl,
     &     vvt,vvb,vvr,vvl,dvdxr,dvdxl,dvdzt,dvdzb,wwr,wwl,wwb,wwbr,wwt,
     &     wwbl,uwr,uwl,uwt,uwb,dudztr,dudzbr,dudzbl,dudztl,dwdxr,
     &     dwdxbr,dwdxl,dwdxbl,uuur,uuul,uuwt,uuwb,uuwr,uuwl,vvur,vvul,
     &     vvwt,vvwb,wwwt,wwwb,wwur,wwul,wwub,wwut,
     &     pur,pul,put,pub,pwt,pwb,pwr,pwl,uurr,uull,wwbb
#ifdef LES
      real nuxz,nuxy,nuyz,dudzc,dwdxc,dudyc,dvdxc,dvdzc,dwdyc,nududx,nuduwdzx,nududx2,nuduwdzxdudz,nuduvdyxdudy
     &    ,ut,ub,nudu2dx,nuduwdzxu,nudvdy2,nuduvdyxdvdx,nudvwdzydvdz,vr,vt,nudvudxyv,nudvwdzyv,nudwdz,nudwdz2,nuduwdzxdwdx
     &    ,nudvwdzydwdy
     &    ,wt,wr,wl,nudw2dz,nuduwdzxw,nududxdwdx,nudwdzdudz,nuduwdzxduwdxz,nuduvdyxdwdy,nudwvdyzdudy,wcor,ucor,nududxw
     &    ,nudwdzu !25
     &    ,nuduwdzx_c,nuduwdzxdudz_c,nuduwdzxdwdx_c,nuduvdyxdudy_c,nuduvdyxdvdx_c,nudvwdzydwdy_c,nudvwdzydvdz_c
     &    ,nuduwdzxu_t,nuduwdzxu_b,nuduwdzx_t,nuduwdzx_b,nuduwdzxw_r,nuduwdzxw_l,nuduwdzx_r,nuduwdzx_l,nuduwdzxu_r
     &    ,nuduwdzxu_l,nuduwdzxw_t,nuduwdzxw_b
#endif

c     samples=samples+1 ! increment samples before calling this routine
c     compute means and stresses
c     
c      write(0,*) "samples: ",samples
      fac=1./(jy2-1)
c     
      do i=2,ix2
         do k=2,kz2

            u=0.
            v=0.
            w=0.
            p=0.
            uu=0.
            uw=0.
            pu=0.
            uuu=0.
            uuw=0.
            pdudx=0.
            dudx2=0.
            dudy2=0.
            dudz2=0.
            	    
c     additional terms for vv
            vv=0.
            vvu=0.
            vvw=0.
            pdvdy=0.
            dvdx2=0.
            dvdy2=0.
            dvdz2=0.
c     additional terms for ww
            ww=0.
            wwu=0.
            www=0.
            pw=0.
            pdwdz=0.
            dwdx2=0.
            dwdy2=0.
            dwdz2=0.
c     additional terms for uw
            duwdxx=0.
            duwdyy=0.
            duwdzz=0.
            pdudz=0.
            pdwdx=0.
#ifdef LES
            nududx=0.
	    nuduwdzx=0.
	    nududx2=0.
	    nuduwdzxdudz=0.
	    nuduvdyxdudy=0.
            nudu2dx=0.
            nuduwdzxu=0.
	    nudvdy2=0.
	    nuduvdyxdvdx=0.
	    nudvwdzydvdz=0.
	    nudvudxyv=0.
	    nudvwdzyv=0.
	    nudwdz=0.
	    nudwdz2=0.
	    nuduwdzxdwdx=0.
	    nudvwdzydwdy=0.
            nudw2dz=0.
	    nuduwdzxw=0.
	    nududxdwdx=0.
	    nudwdzdudz=0.
	    nuduwdzxduwdxz=0.
	    nuduvdyxdwdy=0.
	    nudwvdyzdudy=0.
	    nududxw=0.
            nudwdzu=0.
#endif
            do j=2,jy2

               dvdy=(vo(i,j,k)-vo(i,j-1,k))*bp(j)
               dudx=(uo(i,j,k)-uo(i-1,j,k))*ap(i)
               dwdz=(wo(i,j,k)-wo(i,j,k-1))*cp(k)

               u=u+uo(i,j,k)    !cor
               v=v+vo(i,j,k)    !cen
               w=w+wo(i,j,k)    !cor
               p=p+po(i,j,k)    !cen
               uu=uu+uo(i,j,k)**2 !cor
               uw=uw+0.25*(uo(i,j,k)+uo(i-1,j,k))*(wo(i,j,k)+wo(i,j,k-1)) !cen
               pu=pu+0.5*po(i,j,k)*(uo(i,j,k)+uo(i-1,j,k)) !cen
               uuu=uuu+uo(i,j,k)**3 !cor
               uuw=uuw+0.125*(uo(i,j,k)+uo(i-1,j,k))**2*(wo(i,j,k)+wo(i,j,k-1)) !cen
               pdudx=pdudx+po(i,j,k)*dudx !cen

               dudx2=dudx2+dudx**2 !cen
               dudy2=dudy2+((uo(i,j+1,k)-uo(i,j,k))*bv(j))**2 !cor
               dudz2=dudz2+((uo(i,j,k+1)-uo(i,j,k))*cw(k))**2 !cor-cor

c     additional terms for vv
               vv=vv+vo(i,j,k)**2 !cen
               vvu=vvu+0.125*(uo(i,j,k)+uo(i-1,j,k))*(vo(i,j,k)+vo(i,j-1,k))**2 !cen
               vvw=vvw+0.125*(vo(i,j,k)+vo(i,j-1,k))**2*(wo(i,j,k)+wo(i,j,k-1)) !cen
               pdvdy=pdvdy+po(i,j,k)*dvdy !cen
               dvdy2=dvdy2+dvdy**2 !cen
               dvdx2=dvdx2+((vo(i+1,j,k)-vo(i,j,k))*au(i))**2 !cor
               dvdz2=dvdz2+((vo(i,j,k+1)-vo(i,j,k))*cw(k))**2 !cor               
c     additional terms for ww                
               ww=ww+wo(i,j,k)**2 !cor
               wwu=wwu+0.125*(wo(i,j,k)+wo(i,j,k-1))**2*(uo(i,j,k)+uo(i-1,j,k)) !cen
               www=www+wo(i,j,k)**3 !cor
               pw=pw+0.5*po(i,j,k)*(wo(i,j,k)+wo(i,j,k-1)) !cen
               pdwdz=pdwdz+po(i,j,k)*dwdz
               dwdz2=dwdz2+dwdz**2 !cen
               dwdx2=dwdx2+((wo(i+1,j,k)-wo(i,j,k))*au(i))**2 !cor
               dwdy2=dwdy2+((wo(i,j+1,k)-wo(i,j,k))*bv(j))**2 !cor    
c     additional terms for uw
               dudy=0.25*( (uo(i,j+1,k)-uo(i,j,k)+uo(i-1,j+1,k)-uo(i-1,j,k))*bv(j)+
     &              (uo(i,j,k)-uo(i,j-1,k)+uo(i-1,j,k)-uo(i-1,j-1,k))*bv(j-1) ) !cen
               dudz=0.25*( (uo(i,j,k+1)-uo(i,j,k)+uo(i-1,j,k+1)-uo(i-1,j,k))*cw(k)+
     &              (uo(i,j,k)-uo(i,j,k-1)+uo(i-1,j,k)-uo(i-1,j,k-1))*cw(k-1) ) !cen
               dwdx=0.25*( (wo(i+1,j,k)-wo(i,j,k)+wo(i+1,j,k-1)-wo(i,j,k-1))*au(i)+
     &              (wo(i,j,k)-wo(i-1,j,k)+wo(i,j,k-1)-wo(i-1,j,k-1))*au(i-1) ) !cen
               dwdy=0.25*( (wo(i,j+1,k)-wo(i,j,k)+wo(i,j+1,k-1)-wo(i,j,k-1))*bv(j)+
     &              (wo(i,j,k)-wo(i,j-1,k)+wo(i,j,k-1)-wo(i,j-1,k-1))*bv(j-1) ) !cen
               dvdx=0.25*( (vo(i+1,j,k)-vo(i,j,k)+vo(i+1,j-1,k)-vo(i,j-1,k))*au(i)+
     &	            (vo(i,j,k)-vo(i-1,j,k)+vo(i,j-1,k)-vo(i-1,j-1,k))*au(i-1) ) !cen
               dvdz=0.25*( (vo(i,j,k+1)-vo(i,j,k)+vo(i,j-1,k+1)-vo(i,j-1,k))*cw(k)+
     &	            (vo(i,j,k)-vo(i,j,k-1)+vo(i,j-1,k)-vo(i,j-1,k-1))*cw(k-1) ) !cen
c     additional terms for uw
               duwdxx=duwdxx+dudx*dwdx
               duwdyy=duwdyy+dudy*dwdy
               duwdzz=duwdzz+dudz*dwdz
               pdudz=pdudz+po(i,j,k)*dudz
               pdwdx=pdwdx+po(i,j,k)*dwdx
c
c     subgrid stress terms
#ifdef LES
               nuxz=0.25*(tv(i,j,k)+tv(i,j,k+1)+tv(i+1,j,k)+tv(i+1,j,k+1))!cor-cor
               nuxy=0.25*(tv(i,j,k)+tv(i,j+1,k)+tv(i+1,j,k)+tv(i+1,j+1,k))!cor-cor
               nuyz=0.25*(tv(i,j,k)+tv(i,j+1,k)+tv(i,j,k+1)+tv(i,j+1,k+1))!cor-cor
c
               dudzc=(uo(i,j,k+1)-uo(i,j,k))*cw(k)!cor-cor
	       dwdxc=(wo(i+1,j,k)-wo(i,j,k))*au(i)!cor-cor
               dudyc=(uo(i,j+1,k)-uo(i,j,k))*bv(j)!cor-cor
	       dvdxc=(vo(i+1,j,k)-vo(i,j,k))*au(i)!cor-cor
               dvdzc=(vo(i,j,k+1)-vo(i,j,k))*cw(k)!cor-cor
	       dwdyc=(wo(i,j+1,k)-wo(i,j,k))*bv(j)!cor-cor
c     for uu
c     production
               nududx=nududx+tv(i,j,k)*dudx!cen
	       nuduwdzx=nuduwdzx+nuxz*(dudzc+dwdxc)!cor-cor
c     dissipation
               nududx2=nududx2+tv(i,j,k)*dudx**2
	       nuduwdzxdudz=nuduwdzxdudz+nuxz*(dudzc+dwdxc)*dudzc
	       nuduvdyxdudy=nuduvdyxdudy+nuxy*(dudyc+dvdxc)*dudyc
c     transport
               ut=0.5*(uo(i,j,k)+uo(i,j,k+1))
c	       
               nudu2dx=nudu2dx+tv(i,j,k)*(uo(i,j,k)**2-uo(i-1,j,k)**2)*ap(i)
	       nuduwdzxu=nuduwdzxu+nuxz*(dudzc+dwdxc)*ut
	       
c     for vv
c     no production
c     dissipation 	       
               nudvdy2=nudvdy2+tv(i,j,k)*dvdy**2
	       nuduvdyxdvdx=nuduvdyxdvdx+nuxy*(dvdxc+dudyc)*dvdxc
	       nudvwdzydvdz=nudvwdzydvdz+nuyz*(dvdzc+dwdyc)*dvdzc
c     transport	       
               vr=0.5*(vo(i,j,k)+vo(i+1,j,k))
               vt=0.5*(vo(i,j,k)+vo(i,j,k+1))
	       
	       
	       nudvudxyv=nudvudxyv+nuxy*(dvdxc+dudyc)*vr
	       nudvwdzyv=nudvwdzyv+nuyz*(dvdzc+dwdyc)*vt
c     for ww       	       
c     production   
               nudwdz=nudwdz+tv(i,j,k)*dwdz!cen
c     dissipation
               nudwdz2=nudwdz2+tv(i,j,k)*dwdz**2
	       nuduwdzxdwdx=nuduwdzxdwdx+nuxz*(dudzc+dwdxc)*dwdxc
	       nudvwdzydwdy=nudvwdzydwdy+nuyz*(dvdzc+dwdyc)*dwdyc
c     transport
               wt=0.5*(wo(i,j,k)+wo(i+1,j,k))
c	       
               nudw2dz=nudw2dz+tv(i,j,k)*(wo(i,j,k)**2-wo(i,j,k-1)**2)*cp(k)
      	       nuduwdzxw=nuduwdzxw+nuxz*(dudzc+dwdxc)*wt
               
c     for uw
c     production -- no additional terms
c     dissipation
               nududxdwdx=nududxdwdx+tv(i,j,k)*dudx*dwdx
	       nudwdzdudz=nudwdzdudz+tv(i,j,k)*dwdz*dudz
	       nuduwdzxduwdxz=nuduwdzxduwdxz+tv(i,j,k)*(dudz+dwdx)*(dudx+dwdz)
	       nuduvdyxdwdy=nuduvdyxdwdy+tv(i,j,k)*(dudy+dvdx)*dwdy
	       nudwvdyzdudy=nudwvdyzdudy+tv(i,j,k)*(dwdy+dvdz)*dudy
c     transport      	       
               wcen=0.5*(wo(i,j,k)+wo(i,j,k-1))
               ucen=0.5*(uo(i,j,k)+uo(i-1,j,k))
c      
               nududxw  =nududxw+tv(i,j,k)*dudx*wcen
               nudwdzu  =nudwdzu+tv(i,j,k)*dwdz*ucen
#endif
            enddo
c     divide by the number of samples:
            u=u*fac
            v=v*fac
            w=w*fac
            p=p*fac
            uu=uu*fac
            uw=uw*fac
            pu=pu*fac
            uuu=uuu*fac
            uuw=uuw*fac
            pdudx=pdudx*fac
            dudx2=dudx2*fac
            dudy2=dudy2*fac
            dudz2=dudz2*fac
c     additional terms for vv
            vv=vv*fac
            vvu=vvu*fac
            vvw=vvw*fac
            pdvdy=pdvdy*fac
            dvdx2=dvdx2*fac
            dvdy2=dvdy2*fac
            dvdz2=dvdz2*fac
c     additional terms for ww 
            ww=ww*fac
            wwu=wwu*fac
            www=www*fac
            pw=pw*fac
            pdwdz=pdwdz*fac
            dwdx2=dwdx2*fac
            dwdy2=dwdy2*fac
            dwdz2=dwdz2*fac            
c     additional terms for uw
            duwdxx=duwdxx*fac
            duwdyy=duwdyy*fac
            duwdzz=duwdzz*fac
            pdudz=pdudz*fac
            pdwdx=pdwdx*fac
#ifdef LES
            nududx=nududx*fac
	    nuduwdzx=nuduwdzx*fac
	    nududx2=nududx2*fac
	    nuduwdzxdudz=nuduwdzxdudz*fac
	    nuduvdyxdudy=nuduvdyxdudy*fac
            nudu2dx=nudu2dx*fac
            nuduwdzxu=nuduwdzxu*fac
	    nudvdy2=nudvdy2*fac
	    nuduvdyxdvdx=nuduvdyxdvdx*fac
	    nudvwdzydvdz=nudvwdzydvdz*fac
	    nudvudxyv=nudvudxyv*fac
	    nudvwdzyv=nudvwdzyv*fac
	    nudwdz=nudwdz*fac
	    nudwdz2=nudwdz2*fac
	    nuduwdzxdwdx=nuduwdzxdwdx*fac
	    nudvwdzydwdy=nudvwdzydwdy*fac
            nudw2dz=nudw2dz*fac
	    nuduwdzxw=nuduwdzxw*fac
	    nududxdwdx=nududxdwdx*fac
	    nudwdzdudz=nudwdzdudz*fac
	    nuduwdzxduwdxz=nuduwdzxduwdxz*fac
	    nuduvdyxdwdy=nuduvdyxdwdy*fac
	    nudwvdyzdudy=nudwvdyzdudy*fac
	    nududxw=nududxw*fac
            nudwdzu=nudwdzu*fac
#endif
            
c     now modify the running averages:
            u_(i,k)=u_(i,k)+(u-u_(i,k))/samples
            v_(i,k)=v_(i,k)+(v-v_(i,k))/samples
            w_(i,k)=w_(i,k)+(w-w_(i,k))/samples
            p_(i,k)=p_(i,k)+(p-p_(i,k))/samples
            uu_(i,k)=uu_(i,k)+(uu-uu_(i,k))/samples
            uw_(i,k)=uw_(i,k)+(uw-uw_(i,k))/samples
            pu_(i,k)=pu_(i,k)+(pu-pu_(i,k))/samples
            uuu_(i,k)=uuu_(i,k)+(uuu-uuu_(i,k))/samples
            uuw_(i,k)=uuw_(i,k)+(uuw-uuw_(i,k))/samples
            pdudx_(i,k)=pdudx_(i,k)+(pdudx-pdudx_(i,k))/samples
            dudx2_(i,k)=dudx2_(i,k)+(dudx2-dudx2_(i,k))/samples
            dudy2_(i,k)=dudy2_(i,k)+(dudy2-dudy2_(i,k))/samples
            dudz2_(i,k)=dudz2_(i,k)+(dudz2-dudz2_(i,k))/samples
c     additional terms for vv
            vv_(i,k)=vv_(i,k)+(vv-vv_(i,k))/samples
            vvu_(i,k)=vvu_(i,k)+(vvu-vvu_(i,k))/samples
            vvw_(i,k)=vvw_(i,k)+(vvw-vvw_(i,k))/samples
            pdvdy_(i,k)=pdvdy_(i,k)+(pdvdy-pdvdy_(i,k))/samples
            dvdx2_(i,k)=dvdx2_(i,k)+(dvdx2-dvdx2_(i,k))/samples
            dvdy2_(i,k)=dvdy2_(i,k)+(dvdy2-dvdy2_(i,k))/samples
            dvdz2_(i,k)=dvdz2_(i,k)+(dvdz2-dvdz2_(i,k))/samples            
c     additional terms for ww
            ww_(i,k)=ww_(i,k)+(ww-ww_(i,k))/samples
            wwu_(i,k)=wwu_(i,k)+(wwu-wwu_(i,k))/samples
            www_(i,k)=www_(i,k)+(www-www_(i,k))/samples
            pw_(i,k)=pw_(i,k)+(pw-pw_(i,k))/samples
            pdwdz_(i,k)=pdwdz_(i,k)+(pdwdz-pdwdz_(i,k))/samples
            dwdx2_(i,k)=dwdx2_(i,k)+(dwdx2-dwdx2_(i,k))/samples
            dwdy2_(i,k)=dwdy2_(i,k)+(dwdy2-dwdy2_(i,k))/samples
            dwdz2_(i,k)=dwdz2_(i,k)+(dwdz2-dwdz2_(i,k))/samples            
c     additional terms for uw
            duwdxx_(i,k)=duwdxx_(i,k)+(duwdxx-duwdxx_(i,k))/samples
            duwdyy_(i,k)=duwdyy_(i,k)+(duwdyy-duwdyy_(i,k))/samples
            duwdzz_(i,k)=duwdzz_(i,k)+(duwdzz-duwdzz_(i,k))/samples
            pdudz_(i,k)=pdudz_(i,k)+(pdudz-pdudz_(i,k))/samples
            pdwdx_(i,k)=pdwdx_(i,k)+(pdwdx-pdwdx_(i,k))/samples
#ifdef LES
            nududx_(i,k)      =nududx_(i,k)          +(nududx-nududx_(i,k))/samples
	    nuduwdzx_(i,k)    =nuduwdzx_(i,k)        +(nuduwdzx-nuduwdzx_(i,k))/samples
	    nududx2_(i,k)     =nududx2_(i,k)         +(nududx2-nududx2_(i,k))/samples
	    nuduwdzxdudz_(i,k)=nuduwdzxdudz_(i,k)    +(nuduwdzxdudz-nuduwdzxdudz_(i,k))/samples
	    nuduvdyxdudy_(i,k)=nuduvdyxdudy_(i,k)    +(nuduvdyxdudy-nuduvdyxdudy_(i,k))/samples
            nudu2dx_(i,k)     =nudu2dx_(i,k)         +(nudu2dx-nudu2dx_(i,k))/samples
            nuduwdzxu_(i,k)   =nuduwdzxu_(i,k)       +(nuduwdzxu-nuduwdzxu_(i,k))/samples
	    nudvdy2_(i,k)     =nudvdy2_(i,k)         +(nudvdy2-nudvdy2_(i,k))/samples
	    nuduvdyxdvdx_(i,k)=nuduvdyxdvdx_(i,k)    +(nuduvdyxdvdx-nuduvdyxdvdx_(i,k))/samples
	    nudvwdzydvdz_(i,k)=nudvwdzydvdz_(i,k)    +(nudvwdzydvdz-nudvwdzydvdz_(i,k))/samples
	    nudvudxyv_(i,k)   =nudvudxyv_(i,k)       +(nudvudxyv-nudvudxyv_(i,k))/samples
	    nudvwdzyv_(i,k)   =nudvwdzyv_(i,k)       +(nudvwdzyv-nudvwdzyv_(i,k))/samples
	    nudwdz_(i,k)      =nudwdz_(i,k)          +(nudwdz-nudwdz_(i,k))/samples
	    nudwdz2_(i,k)     =nudwdz2_(i,k)         +(nudwdz2-nudwdz2_(i,k))/samples
	    nuduwdzxdwdx_(i,k)=nuduwdzxdwdx_(i,k)    +(nuduwdzxdwdx-nuduwdzxdwdx_(i,k))/samples
	    nudvwdzydwdy_(i,k)=nudvwdzydwdy_(i,k)    +(nudvwdzydwdy-nudvwdzydwdy_(i,k))/samples
            nudw2dz_(i,k)     =nudw2dz_(i,k)         +(nudw2dz-nudw2dz_(i,k))/samples
	    nuduwdzxw_(i,k)   =nuduwdzxw_(i,k)       +(nuduwdzxw-nuduwdzxw_(i,k))/samples
	    nududxdwdx_(i,k)  =nududxdwdx_(i,k)      +(nududxdwdx-nududxdwdx_(i,k))/samples
	    nudwdzdudz_(i,k)  =nudwdzdudz_(i,k)      +(nudwdzdudz-nudwdzdudz_(i,k))/samples
	    nuduwdzxduwdxz_(i,k)=nuduwdzxduwdxz_(i,k)+(nuduwdzxduwdxz-nuduwdzxduwdxz_(i,k))/samples
	    nuduvdyxdwdy_(i,k)=nuduvdyxdwdy_(i,k)    +(nuduvdyxdwdy-nuduvdyxdwdy_(i,k))/samples
	    nudwvdyzdudy_(i,k)=nudwvdyzdudy_(i,k)    +(nudwvdyzdudy-nudwvdyzdudy_(i,k))/samples
	    nududxw_(i,k)     =nududxw_(i,k)         +(nududxw-nududxw_(i,k))/samples
            nudwdzu_(i,k)     =nudwdzu_(i,k)         +(nudwdzu-nudwdzu_(i,k))/samples
#endif

         enddo
      enddo
c     update running averages on outer boundary points for u,uu,uuu,w,ww,www
c     u first
      do i=1,ix2+1,ix2
       do k=2,kz2
        u=0.
	uu=0.
	uuu=0.
	do j=2,jy2
	  u=u+uo(i,j,k)
	  uu=uu+uo(i,j,k)**2
	  uuu=uuu+uo(i,j,k)**3
	enddo
        u=u*fac
        uu=uu*fac
        uuu=uuu*fac
        u_(i,k)=u_(i,k)+(u-u_(i,k))/samples
        uu_(i,k)=uu_(i,k)+(uu-uu_(i,k))/samples
        uuu_(i,k)=uuu_(i,k)+(uuu-uuu_(i,k))/samples
       enddo
      enddo 	 
c     now w
      if (myrank.eq.0) then 
        k=1
	do i=2,ix2
          w=0.
  	  ww=0.
  	  www=0.
	  do j=2,jy2
	    w=w+wo(i,j,k)
	    ww=ww+wo(i,j,k)**2
	    www=www+wo(i,j,k)**3
	  enddo
          w=w*fac
          ww=ww*fac
          www=www*fac
          w_(i,k)=w_(i,k)+(w-w_(i,k))/samples
          ww_(i,k)=ww_(i,k)+(ww-ww_(i,k))/samples
          www_(i,k)=www_(i,k)+(www-www_(i,k))/samples
	enddo
      endif	  
c
      if (myrank.eq.mysize-1) then 
        k=kz2+1
	do i=2,ix2
          w=0.
  	  ww=0.
  	  www=0.
	  do j=2,jy2
	    w=w+wo(i,j,k)
	    ww=ww+wo(i,j,k)**2
	    www=www+wo(i,j,k)**3
	  enddo
          w=w*fac
          ww=ww*fac
          www=www*fac
          w_(i,k)=w_(i,k)+(w-w_(i,k))/samples
          ww_(i,k)=ww_(i,k)+(ww-ww_(i,k))/samples
          www_(i,k)=www_(i,k)+(www-www_(i,k))/samples
	enddo
      endif	  
	             
c     refresh data on boundaries (nz has to ba a multiple of 2)
c     that's because two columns are sent, not the usual one, 
c     and this requires redimensioning the arrays as below 
      call REFRESHBC(u_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(v_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(w_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(p_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(uu_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(uw_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(pu_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(uuu_(1,0),   2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(uuw_(1,0),   2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(pdudx_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dudx2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dudy2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dudz2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(vv_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(vvu_(1,0),   2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(vvw_(1,0),   2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(pdvdy_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dvdx2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dvdy2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dvdz2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(ww_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(wwu_(1,0),   2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(www_(1,0),   2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(pw_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(pdwdz_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dwdx2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dwdy2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(dwdz2_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(duwdxx_(1,0),2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(duwdyy_(1,0),2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(duwdzz_(1,0),2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(pdudz_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(pdwdx_(1,0), 2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
#ifdef LES
      call REFRESHBC(nududx_(1,0),        2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduwdzx_(1,0),      2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nududx2_(1,0),       2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduwdzxdudz_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduvdyxdudy_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudu2dx_(1,0),       2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduwdzxu_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudvdy2_(1,0),       2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduvdyxdvdx_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudvwdzydvdz_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudvudxyv_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudvwdzyv_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudwdz_(1,0),        2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudwdz2_(1,0),       2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduwdzxdwdx_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudvwdzydwdy_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudw2dz_(1,0),       2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduwdzxw_(1,0),     2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nududxdwdx_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudwdzdudz_(1,0),    2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduwdzxduwdxz_(1,0),2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nuduvdyxdwdy_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudwvdyzdudy_(1,0),  2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nududxw_(1,0),       2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
      call REFRESHBC(nudwdzu_(1,0),       2*nx,(nz-2)/2,mysize,myrank,MPI_COMM_HYDRO)
c
#endif
c
c     assign starting and ending indices
      if (myrank.eq.0) then 
        kbg=3
      else
        kbg=2
      endif
c      
      if (myrank.eq.mysize-1) then 
        knd=kz2-1
      else
        knd=kz2
      endif
      	 
c     now build the terms
c 
      do i=3,ix2-1
         do k=kbg,knd

c     no y terms (spanwise homogeneity)

c     convective term    
            uutr=uu_(i,k+1)-u_(i,k+1)**2
            uutl=uu_(i-1,k+1)-u_(i-1,k+1)**2
            uur=uu_(i,k)-u_(i,k)**2
            uul=uu_(i-1,k)-u_(i-1,k)**2
            uubr=uu_(i,k-1)-u_(i,k-1)**2
            uubl=uu_(i-1,k-1)-u_(i-1,k-1)**2
c     
            ucent=0.5*(u_(i,k+1)+u_(i-1,k+1))
            ucen= 0.5*(u_(i,k)+u_(i-1,k))
            ucenr= 0.5*(u_(i+1,k)+u_(i,k))
            ucenl= 0.5*(u_(i-1,k)+u_(i-2,k))
            ucenb=0.5*(u_(i,k-1)+u_(i-1,k-1))
            wcen= 0.5*(w_(i,k)+w_(i,k-1))
            wcent= 0.5*(w_(i,k+1)+w_(i,k))
            wcenb= 0.5*(w_(i,k-1)+w_(i,k-2))

            vcen=v_(i,k)
            vcenr=v_(i+1,k)
            vcenl=v_(i-1,k)
            vcent=v_(i,k+1)
            vcenb=v_(i,k-1)

            wcenr= 0.5*(w_(i+1,k)+w_(i+1,k-1))
            wcenl= 0.5*(w_(i-1,k)+w_(i-1,k-1))
c     
            vv=vv_(i,k)-v_(i,k)**2      
            vvt=vv_(i,k+1)-v_(i,k+1)**2      
            vvb=vv_(i,k-1)-v_(i,k-1)**2      
            vvr=vv_(i+1,k)-v_(i+1,k)**2      
            vvl=vv_(i-1,k)-v_(i-1,k)**2      
c     
            ww=ww_(i,k)-w_(i,k)**2
            wwr=ww_(i+1,k)-w_(i+1,k)**2
            wwl=ww_(i-1,k)-w_(i-1,k)**2
            wwb=ww_(i,k-1)-w_(i,k-1)**2
            wwbr=ww_(i+1,k-1)-w_(i+1,k-1)**2
            wwbl=ww_(i-1,k-1)-w_(i-1,k-1)**2
c     
            uwr=uw_(i+1,k)-ucenr*wcenr
            uw=uw_(i,k)-ucen*wcen
            uwl=uw_(i-1,k)-ucenl*wcenl
            uwt=uw_(i,k+1)-ucent*wcent
            uwb=uw_(i,k-1)-ucenb*wcenb
c     
            dudx=(u_(i,k)-u_(i-1,k))*ap(i)
            dwdz=(w_(i,k)-w_(i,k-1))*cp(k)

            dudztr=(u_(i,k+1)-u_(i,k))*cw(k)
            dudzbr=(u_(i,k)-u_(i,k-1))*cw(k-1)
            dudztl=(u_(i-1,k+1)-u_(i-1,k))*cw(k)
            dudzbl=(u_(i-1,k)-u_(i-1,k-1))*cw(k-1)
c     
            dwdxr=(w_(i+1,k)-w_(i,k))*au(i)
            dwdxbr=(w_(i+1,k-1)-w_(i,k-1))*au(i)
            dwdxl=(w_(i,k)-w_(i-1,k))*au(i-1)
            dwdxbl=(w_(i,k-1)-w_(i-1,k-1))*au(i-1)

c     uu
            uu_conv(i,k)=ucen*(uur-uul)*ap(i)+
     &           0.25*(w_(i,k)*(uutr-uur+uutl-uul)*cw(k)+w_(i,k-1)*(uur-uubr+uul-uubl)*cw(k-1))
c     vv
            vv_conv(i,k)=0.5*(u_(i,k)*(vvr-vv)*au(i)+
     &           u_(i-1,k)*(vv-vvl)*au(i-1))+
     &           0.5*(w_(i,k)*(vvt-vv)*cw(k)+w_(i,k-1)*(vv-vvb)*cw(k-1))
c     ww
            ww_conv(i,k)=0.25*(u_(i,k)*(wwr-ww+wwbr-wwb)*au(i)+u_(i-1,k)*(ww-wwl+wwb-wwbl)*au(i-1))+
     &           wcen*(ww-wwb)*cp(k)
c     uw
            uw_conv(i,k)=0.5*(u_(i,k)*(uwr-uw)*au(i)+u_(i-1,k)*(uw-uwl)*au(i-1))+
     &           0.5*(w_(i,k)*(uwt-uw)*cw(k)+w_(i,k-1)*(uw-uwb)*cw(k-1))
c     production term
c     uu
            uu_prod(i,k)=2.*(0.5*(uur+uul)*dudx+
     &           0.5*uw*((ucent-ucen)*cw(k)+(ucen-ucenb)*cw(k-1)))
c     vv
            vv_prod(i,k)=0.

c     ww
            ww_prod(i,k)=2.*(0.5*uw*((wcenr-wcen)*au(i)+(wcen-wcenl)*au(i-1))+
     &           0.5*(ww+wwb)*dwdz)

            dudz=0.5*((ucent-ucen)*cw(k)+(ucen-ucenb)*cw(k-1))
            dvdz=0.5*((vcent-vcen)*cw(k)+(vcen-vcenb)*cw(k-1))
            dvdx=0.5*((vcenr-vcen)*au(i)+(vcen-vcenl)*au(i-1))
            dwdx=0.5*((wcenr-wcen)*au(i)+(wcen-wcenl)*au(i-1))

c     uw
            uw_prod(i,k)=uw*dudx+0.25*(ww*(dudztr+dudztl)+wwb*(dudzbr+dudzbl))+
     &           0.25*(uur*(dwdxr+dwdxbr)+uul*(dwdxl+dwdxbl))+uw*dwdz

c     dissipation term 
            
c     uu      
            uu_diss(i,k)=2.*ru1*((dudx2_(i,k)-dudx**2)+
     &           0.5*(dudy2_(i,k)+dudy2_(i-1,k))+ ! no mean spanwise derivative
     &           0.25*(dudz2_(i,k)+dudz2_(i,k-1)+dudz2_(i-1,k)+dudz2_(i-1,k-1)
     &           -dudztr**2-dudzbr**2-dudztl**2-dudzbl**2))
c     vv
            vv_diss(i,k)=2.*ru1*(0.5*(dvdx2_(i,k)+dvdx2_(i-1,k))+
     &           dvdy2_(i,k)+   ! no mean spanwise derivative
     &           0.5*(dvdz2_(i,k)+dvdz2_(i,k-1)))
c     ww
            ww_diss(i,k)=2.*ru1*(0.25*(dwdx2_(i,k)+dwdx2_(i-1,k)+dwdx2_(i,k-1)+dwdx2_(i-1,k-1)
     &           -dwdxr**2-dwdxbr**2-dwdxl**2-dwdxbl**2)+
     &           0.5*(dwdy2_(i,k)+dwdy2_(i,k-1))+
     &           dwdz2_(i,k)-dwdz**2) 
c     uw
            uw_diss(i,k)=2.*ru1*(duwdxx_(i,k)-dudx*dwdx+duwdyy_(i,k)+duwdzz_(i,k)-dudz*dwdz)

c     turbulent transport       
            uuur=uuu_(i,k)-3.*uu_(i,k)*u_(i,k)+2.*u_(i,k)**3
            uuul=uuu_(i-1,k)-3.*uu_(i-1,k)*u_(i-1,k)+2.*u_(i-1,k)**3
c     
            uuwt=uuw_(i,k+1)-2.*ucent*uw_(i,k+1)+(2.*ucent**2-0.5*(uu_(i,k+1)+uu_(i-1,k+1)))*wcent
            uuw=uuw_(i,k)-2.*ucen*uw_(i,k)+(2.*ucen**2-0.5*(uu_(i,k)+uu_(i-1,k)))*wcen
            uuwb=uuw_(i,k-1)-2.*ucenb*uw_(i,k-1)+(2.*ucenb**2-0.5*(uu_(i,k-1)+uu_(i-1,k-1)))*wcenb
            uuwr=uuw_(i+1,k)-2.*ucenr*uw_(i+1,k)+(2.*ucenr**2-0.5*(uu_(i+1,k)+uu_(i,k)))*wcenr
            uuwl=uuw_(i-1,k)-2.*ucenl*uw_(i-1,k)+(2.*ucenl**2-0.5*(uu_(i-1,k)+uu_(i-2,k)))*wcenl
c     
            vvur=vvu_(i+1,k)+2.*v_(i+1,k)**2*ucenr-vv_(i+1,k)*ucenr
            vvu=vvu_(i,k)+2.*v_(i,k)**2*ucen-vv_(i,k)*ucen
            vvul=vvu_(i-1,k)+2.*v_(i-1,k)**2*ucenl-vv_(i-1,k)*ucenl
            vvw=vvw_(i,k)+2.*v_(i,k)**2*wcen-vv_(i,k)*wcen
            vvwt=vvw_(i,k+1)+2.*v_(i,k+1)**2*wcent-vv_(i,k+1)*wcent
            vvwb=vvw_(i,k-1)+2.*v_(i,k-1)**2*wcenb-vv_(i,k-1)*wcenb
c     
            wwwt=www_(i,k)-3.*ww_(i,k)*w_(i,k)+2.*w_(i,k)**3
            wwwb=www_(i,k-1)-3.*ww_(i,k-1)*w_(i,k-1)+2.*w_(i,k-1)**3
            wwur=wwu_(i+1,k)-2.*wcenr*uw_(i+1,k)+(2.*wcenr**2-0.5*(ww_(i+1,k)+ww_(i+1,k-1)))*ucenr
            wwu=wwu_(i,k)-2.*wcen*uw_(i,k)+(2.*wcen**2-0.5*(ww_(i,k)+ww_(i,k-1)))*ucen
            wwul=wwu_(i-1,k)-2.*wcenl*uw_(i-1,k)+(2.*wcenl**2-0.5*(ww_(i-1,k)+ww_(i-1,k-1)))*ucenl
            wwut=wwu_(i,k+1)-2.*wcent*uw_(i,k+1)+(2.*wcent**2-0.5*(ww_(i,k+1)+ww_(i,k)))*ucent
            wwub=wwu_(i,k-1)-2.*wcenb*uw_(i,k-1)+(2.*wcenb**2-0.5*(ww_(i,k-1)+ww_(i,k-2)))*ucenb
c     
c     uu
c     
            uu_ttrans(i,k)=(uuur-uuul)*ap(i)+0.5*((uuwt-uuw)*cw(k)+(uuw-uuwb)*cw(k-1))
c     vv
            vv_ttrans(i,k)=0.5*( (vvur-vvu)*au(i)+(vvu-vvul)*au(i-1)+
     &           (vvwt-vvw)*cw(k)+(vvw-vvwb)*cw(k-1) )
c     ww
            ww_ttrans(i,k)=0.5*( (wwur-wwu)*au(i)+(wwu-wwul)*au(i-1) )+(wwwt-wwwb)*cp(k)
c     uw
            uw_ttrans(i,k)=0.5*( (uuwr-uuw)*au(i)+(uuw-uuwl)*au(i-1) +
     &           (wwut-wwu)*cw(k)+(wwu-wwub)*cw(k-1) )
c     pressure transport
c     uu
            pur=pu_(i+1,k)-p_(i+1,k)*ucenr
            pu=pu_(i,k)-p_(i,k)*ucen
            pul=pu_(i-1,k)-p_(i-1,k)*ucenl

            put=pu_(i,k+1)-p_(i,k+1)*ucent
            pub=pu_(i,k-1)-p_(i,k-1)*ucenb

            uu_ptrans(i,k)=((pur-pu)*au(i)+(pu-pul)*au(i-1))
c     vv -- no terms
            vv_ptrans(i,k)=0.
c     ww
            pwt=pw_(i,k+1)-p_(i,k+1)*wcent
            pwb=pw_(i,k-1)-p_(i,k-1)*wcenb
            pwr=pw_(i+1,k)-p_(i+1,k)*wcenr
            pwl=pw_(i-1,k)-p_(i-1,k)*wcenl
            pw=pw_(i,k)-p_(i,k)*wcen

            ww_ptrans(i,k)=(pwt-pw)*cw(k)+(pw-pwb)*cw(k-1)

c     uw
            uw_ptrans(i,k)=0.5*( (put-pu)*cw(k)+(pu-pub)*cw(k-1)+(pwr-pw)*au(i)+(pw-pwl)*au(i-1))

c     pressure strain
c     uu
            uu_pstr(i,k)=-2.*(pdudx_(i,k)-p_(i,k)*dudx)
c     vv
            vv_pstr(i,k)=-2.*pdvdy_(i,k) ! no mean spanwise derivative
c     ww
            ww_pstr(i,k)=-2.*(pdwdz_(i,k)-p_(i,k)*dwdz)
c     uw
            uw_pstr(i,k)=-(pdudz_(i,k)-p_(i,k)*dudz+pdwdx_(i,k)-p_(i,k)*dwdx)

c     viscous diffusion
            uurr=uu_(i+1,k)-u_(i+1,k)**2
            uull=uu_(i-2,k)-u_(i-2,k)**2
            wwt=ww_(i,k+1)-w_(i,k+1)**2
            wwbb=ww_(i,k-2)-w_(i,k-2)**2

c     uu    
            uu_vdiff(i,k)=-ru1*(((uurr-uur)*ap(i+1)-(uul-uull)*ap(i-1))*app(i)*0.5+
     &           0.5*((uutr+uutl-uur-uul)*cw(k)-(uur+uul-uubr-uubl)*cw(k-1))*cpp(k))
c     vv
            vv_vdiff(i,k)=-ru1*(((vvr-vv)*au(i)-(vv-vvl)*au(i-1))*app(i)+
     &           ((vvt-vv)*cw(k)-(vv-vvb)*cw(k-1))*cpp(k))
c     ww    
            ww_vdiff(i,k)=-ru1*(0.5*((wwr+wwbr-ww-wwb)*au(i)-(ww+wwb-wwl-wwbl)*au(i-1))*app(i)+
     &           ((wwt-ww)*cp(k+1)-(wwb-wwbb)*cp(k-1))*cpp(k)*0.5)
c     uw
            uw_vdiff(i,k)=-ru1*(((uwr-uw)*au(i)-(uw-uwl)*au(i-1))*app(i)+
     &           ((uwt-uw)*cw(k)-(uw-uwb)*cw(k-1))*cpp(k))
#ifdef LES
c     sgs production
c     uu
            nuduwdzx_c=0.25*(nuduwdzx_(i-1,k-1)+nuduwdzx_(i,k-1)+nuduwdzx_(i-1,k)+nuduwdzx_(i,k))

	    uu_prod_sgs(i,k)=-4*nududx_(i,k)*dudx-2*nuduwdzx_c*dudz
c     vv	    
	    vv_prod_sgs(i,k)=0.
c     ww    
            ww_prod_sgs(i,k)=-4.*nudwdz_(i,k)*dwdz-2*nuduwdzx_c*dwdx
c     uw  
            uw_prod_sgs(i,k)=-2.*nududx_(i,k)*dwdx-2*nudwdz_(i,k)*dwdz-nuduwdzx_c*(dudx+dwdz)
c     sgs dissipation
            nuduwdzxdudz_c=0.25*(nuduwdzxdudz_(i-1,k-1)+nuduwdzxdudz_(i,k-1)+nuduwdzxdudz_(i-1,k)+nuduwdzxdudz_(i,k))
            nuduwdzxdwdx_c=0.25*(nuduwdzxdwdx_(i-1,k-1)+nuduwdzxdwdx_(i,k-1)+nuduwdzxdwdx_(i-1,k)+nuduwdzxdwdx_(i,k))

            nuduvdyxdudy_c= 0.5*(nuduvdyxdudy_(i-1,k)+nuduvdyxdudy_(i,k))
            nuduvdyxdvdx_c= 0.5*(nuduvdyxdvdx_(i-1,k)+nuduvdyxdvdx_(i,k))
      
            nudvwdzydwdy_c= 0.5*(nudvwdzydwdy_(i,k-1)+nudvwdzydwdy_(i,k))
            nudvwdzydvdz_c= 0.5*(nudvwdzydvdz_(i,k-1)+nudvwdzydvdz_(i,k))
c     uu
            uu_diss_sgs(i,k)=-(4.*nududx2_(i,k)+2.*nuduwdzxdudz_c+2.*nuduvdyxdudy_c)
c     vv                
            vv_diss_sgs(i,k)=-(4.*nudvdy2_(i,k)+2.*nuduvdyxdvdx_c+2.*nudvwdzydvdz_c)
c     ww
            ww_diss_sgs(i,k)=-(4.*nudwdz2_(i,k)+2.*nuduwdzxdwdx_c+2.*nudvwdzydwdy_c)
c     uw  
            uw_diss_sgs(i,k)=-(2.*nududxdwdx_(i,k)+2.*nudwdzdudz_(i,k)+nuduwdzxduwdxz_(i,k)+
     &                          nuduvdyxdwdy_(i,k)+nudwvdyzdudy_(i,k))	    	    
c     sgs transport
            nuduwdzxu_t=0.5*(nuduwdzxu_(i-1,k)+nuduwdzxu_(i,k))
            nuduwdzxu_b=0.5*(nuduwdzxu_(i-1,k-1)+nuduwdzxu_(i,k-1))
            nuduwdzx_t=0.5*(nuduwdzx_(i-1,k)+nuduwdzx_(i,k))
            nuduwdzx_b=0.5*(nuduwdzx_(i-1,k-1)+nuduwdzx_(i,k-1))
c
            nuduwdzxw_r=0.5*(nuduwdzxw_(i,k-1)+nuduwdzxw_(i,k))
            nuduwdzxw_l=0.5*(nuduwdzxw_(i-1,k-1)+nuduwdzxw_(i-1,k))
            nuduwdzx_r=0.5*(nuduwdzx_(i,k-1)+nuduwdzx_(i,k))
            nuduwdzx_l=0.5*(nuduwdzx_(i-1,k-1)+nuduwdzx_(i-1,k))
c
            nuduwdzxu_r=0.5*(nuduwdzxu_(i,k-1)+nuduwdzxu_(i,k))
            nuduwdzxu_l=0.5*(nuduwdzxu_(i-1,k-1)+nuduwdzxu_(i-1,k))
c
            nuduwdzxw_t=0.5*(nuduwdzxw_(i-1,k)+nuduwdzxw_(i,k))
            nuduwdzxw_b=0.5*(nuduwdzxw_(i-1,k-1)+nuduwdzxw_(i,k-1))
	    
	    ut=0.5*(ucen+ucent)
	    ub=0.5*(ucen+ucenb)
c
	    wr=0.5*(wcen+wcenr)
	    wl=0.5*(wcen+wcenl)
c     uu
            uu_tran_sgs(i,k)=2.*(nudu2dx_(i+1,k)-nudu2dx_(i-1,k))*app(i)*0.5+
     &                       2.*(nuduwdzxu_t-nuduwdzxu_b)*cp(k)-
c     means
     &                       4.*(nududx_(i+1,k)*ucenr-nududx_(i-1,k)*ucenl)*app(i)*0.5-
     &                       2.*(nuduwdzx_t*ut-nuduwdzx_b*ub)*cp(k)     	    

c     vv
            vv_tran_sgs(i,k)=2.*(nudvudxyv_(i,k)-nudvudxyv_(i-1,k))*ap(i)+
     &                       2.*(nudvwdzyv_(i,k)-nudvwdzyv_(i,k-1))*cp(k)
c     ww
            ww_tran_sgs(i,k)=2.*(nuduwdzxw_r-nuduwdzxw_l)*ap(i)+
     &	                     2.*(nudw2dz_(i,k+1)-nudw2dz_(i,k-1))*cpp(k)*0.5-
c     means
     &                       2.*(nuduwdzx_r*wr-nuduwdzx_l*wl)*ap(i)-
     &                       4.*(nudwdz_(i,k+1)*wcent-nudwdz_(i,k-1)*wcenb)*cpp(k)*0.5
c     uw
            uw_tran_sgs(i,k)=2.*(nududxw_(i+1,k)-nududxw_(i-1,k))*app(i)*0.5+
     &	                        (nuduwdzxu_r-nuduwdzxu_l)*ap(i)+
     &				(nuduwdzxw_t-nuduwdzxw_b)*cp(k)+
     &                       2.*(nudwdzu_(i,k+1)-nudwdzu_(i,k-1))*cpp(k)*0.5-
c     means
     &                       2.*(nududx_(i+1,k)*wcenr-nududx_(i-1,k)*wcenl)*app(i)*0.5-
     &                          (nuduwdzx_r*u_(i,k)-nuduwdzx_l*u_(i-1,k))*ap(i)-
     &                          (nuduwdzx_t*w_(i,k)-nuduwdzx_b*w_(i,k-1))*cp(k)-
     &                       2.*(nudwdz_(i,k+1)*ucent-nudwdz_(i,k-1)*ucenb)*cpp(k)*0.5            			     	    

#endif
         enddo
      enddo 
      return
      end subroutine rsbud










ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        compute stress budget data
      subroutine rsbud_compute(rsblk_1,rsblkg_1,rsblk_2,
#ifdef LES      
     &                   rsblk_3,rsblk_4,rsblkg_3,
#endif     
     &                   uo,vo,wo,p,tv,utaug,tauwg,samples,statnum,rsnames,xcf,zcg,stat2d)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'

      integer i,j,k,n
      real rsblk_1(nx,nz,ny)
      real rsblk_2(nx,0:nz+1,33)
      real rsblkg_1(nx,nzg,7) 
#ifdef LES
      real rsblk_3 (nx,nz,ny)     !temporary LES block
      real rsblk_4 (nx,0:nz+1,25) !permanent LES block
      real rsblkg_3(nx,nzg,4)     !temporary global LES block
#endif
      integer stat2d
      real tauwg(nzg), utaug(nzg)
      character* rsnames(0:3))


         call rsbud(rsblk_1(1,1,1),rsblk_1(1,1,2),rsblk_1(1,1,3),rsblk_1(1,1,4),rsblk_1(1,1,5),rsblk_1(1,1,6),rsblk_1(1,1,7),
     &              rsblk_1(1,1,8),rsblk_1(1,1,9),rsblk_1(1,1,10),rsblk_1(1,1,11),rsblk_1(1,1,12),rsblk_1(1,1,13),rsblk_1(1,1,14),
     &              rsblk_1(1,1,15),rsblk_1(1,1,16),rsblk_1(1,1,17),rsblk_1(1,1,18),rsblk_1(1,1,19),rsblk_1(1,1,20),rsblk_1(1,1,21),
     &              rsblk_1(1,1,22),rsblk_1(1,1,23),rsblk_1(1,1,24),rsblk_1(1,1,25),rsblk_1(1,1,26),rsblk_1(1,1,27),rsblk_1(1,1,28),
#ifdef LES
     &              rsblk_3(1,1,1),rsblk_3(1,1,2),rsblk_3(1,1,3),rsblk_3(1,1,4),rsblk_3(1,1,5),rsblk_3(1,1,6),rsblk_3(1,1,7),
     &              rsblk_3(1,1,8),rsblk_3(1,1,9),rsblk_3(1,1,10),rsblk_3(1,1,11),rsblk_3(1,1,12),!12
#endif     
     &              rsblk_2(1,0,1),rsblk_2(1,0,2),rsblk_2(1,0,3),rsblk_2(1,0,4),rsblk_2(1,0,5),rsblk_2(1,0,6),rsblk_2(1,0,7),
     &              rsblk_2(1,0,8),rsblk_2(1,0,9),rsblk_2(1,0,10),rsblk_2(1,0,11),rsblk_2(1,0,12),rsblk_2(1,0,13),rsblk_2(1,0,14),    	           
     &              rsblk_2(1,0,15),rsblk_2(1,0,16),rsblk_2(1,0,17),rsblk_2(1,0,18),rsblk_2(1,0,19),rsblk_2(1,0,20),rsblk_2(1,0,21),
     &              rsblk_2(1,0,22),rsblk_2(1,0,23),rsblk_2(1,0,24),rsblk_2(1,0,25),rsblk_2(1,0,26),rsblk_2(1,0,27),rsblk_2(1,0,28),
     &              rsblk_2(1,0,29),rsblk_2(1,0,30),rsblk_2(1,0,31),rsblk_2(1,0,32),rsblk_2(1,0,33),!33
#ifdef LES
     &              rsblk_4(1,0,1),rsblk_4(1,0,2),rsblk_4(1,0,3),rsblk_4(1,0,4),rsblk_4(1,0,5),rsblk_4(1,0,6),rsblk_4(1,0,7),
     &              rsblk_4(1,0,8),rsblk_4(1,0,9),rsblk_4(1,0,10),rsblk_4(1,0,11),rsblk_4(1,0,12),rsblk_4(1,0,13),rsblk_4(1,0,14),    	           
     &              rsblk_4(1,0,15),rsblk_4(1,0,16),rsblk_4(1,0,17),rsblk_4(1,0,18),rsblk_4(1,0,19),rsblk_4(1,0,20),rsblk_4(1,0,21),
     &              rsblk_4(1,0,22),rsblk_4(1,0,23),rsblk_4(1,0,24),rsblk_4(1,0,25),
#endif
     &              uo,vo,wo,p,
#ifdef LES
     &              tv,
#endif          
     &              samples)
c
         do n=0,3     
	  do j=1,7
	    call MPI_GATHER(rsblk_1(1,2,n*7+j),1,stat2d,rsblkg_1(1,2,j),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
	  enddo
#ifdef LES
	  do j=1,3
	    call MPI_GATHER(rsblk_3(1,2,n*3+j),1,stat2d,rsblkg_3(1,2,j),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
	  enddo
#endif
c         matlab output 
            if (myrank.eq.mysize-1) then 
c	    
               open(62,file=rsnames(n)//index(statnum)//'.dat',status='unknown',form='unformatted')
               write(62) ru1
               write(62) nx-4,nzg-4 
               write(62) (xcf(i),i=3,ix2-1)		 
               write(62) (zcg(k),k=3,nzg-2)
               do k=3,nzg-2
                write(62) utaug(k),(
     &              -rsblkg_1(i,k,1)/tauwg(k)**2*ru1,-rsblkg_1(i,k,2)/tauwg(k)**2*ru1,
     &              -rsblkg_1(i,k,3)/tauwg(k)**2*ru1,-rsblkg_1(i,k,4)/tauwg(k)**2*ru1,
     &              -rsblkg_1(i,k,5)/tauwg(k)**2*ru1,-rsblkg_1(i,k,6)/tauwg(k)**2*ru1,
     &              -rsblkg_1(i,k,7)/tauwg(k)**2*ru1,
#ifdef LES
     &              -rsblkg_3(i,k,1)/tauwg(k)**2*ru1,rsblkg_3(i,k,2)/tauwg(k)**2*ru1,
     &               rsblkg_3(i,k,3)/tauwg(k)**2*ru1,
#endif
     &                 i=3,ix2-1)
               enddo
               close(62)
            endif !myrank
          enddo !n
      return
      end subroutine rsbud_compute
ccccccccccccccccccccccccccccccccREYNOLDS STRESS BUDGETS -- SEE 'RSBUD1.F cccccccccccccccccccccccccccccccccccccccccccccccccccc
