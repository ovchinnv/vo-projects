

       integer statini, subtract_meanvel
       parameter (statini=0, subtract_meanvel=0)



      integer statnum,samples,psample



      real uinft(nz),uinftg(nzg),ubulk(nz),cf(nz),dst(nz),theta(nz),winst(nz),tauwg(nzg),utaug(nzg),ubulkg(nzg),dstg(nzg)
      real thetag(nzg),cfg(nzg),kpg(nzg),winstg(nzg)
      real uavg(nx,nz),vavg(nx,nz),wavg(nx,nz),uuavg(nx,nz),vvavg(nx,nz),wwavg(nx,nz),uwavg(nx,nz),uuuavg(nx,nz),
     &     wwuavg(nx,nz) 
!      real omxavg(nx,nz),omyavg(nx,nz),omzavg(nx,nz)
      real wtemp(ny,nz),pavg(nx,nz),ppavg(nx,nz)
      real uuav(nz),vvav(nz),wwav(nz) ! real stat arrays
      real uav(nz),vav(nz),wav(nz) ! real stat arrays
      real guav(nzg),gvav(nzg),gwav(nzg) ! real stat arrays
      real guuav(nzg),gvvav(nzg),gwwav(nzg) ! real stat arrays
      real usave(nx,ny),vsave(nx,ny),wsave(nx,ny)




      call init_stat(samples,statnum,uinft,uav,vav,wav,uuav,vvav,wwav,
     &  uavg,vavg,wavg,pavg,uuavg,vvavg,wwavg,ppavg,uwavg,uuuavg,wwuavg,
#ifdef LES
     & tvavg,
#endif
     & rsblk_2, 
#ifdef LES     
     & rsblk_4,
#endif
     & rsnames)







c    added 7/2005

      subroutine init_stat(samples,statnum,uinft,uav,vav,wav,uuav,vvav,wwav,
     &  uavg,vavg,wavg,pavg,uuavg,vvavg,wwavg,ppavg,uwavg,uuuavg,wwuavg,
#ifdef LES
     & tvavg,
#endif
     & rsblk_2, 
#ifdef LES     
     & rsblk_4,
#endif
     & rsnames)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     initialize some statistics variables      
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/stat.h'
      
      integer samples, statnum
      
      real pavg(nx,nz),ppavg(nx,nz)   
      real uavg(nx,nz),vavg(nx,nz),wavg(nx,nz),uuavg(nx,nz),vvavg(nx,nz),wwavg(nx,nz),uwavg(nx,nz),
     &     uuuavg(nx,nz), wwuavg(nx,nz) 

      real uinft(nz)
      real uuav(nz),vvav(nz),wwav(nz)
      real uav(nz),vav(nz),wav(nz)   
#ifdef LES
      real tvavg(nx,nz)
#endif

      real rsblk_2(nx,0:nz+1,33)
#ifdef LES
      real rsblk_4 (nx,0:nz+1,25)     !temporary LES block
#endif
c
      character*4 rsnames(0:3)
      integer i,j,k
c
c  
      samples=0
      statnum=statini
      do k=kz1,kz2+1
        uinft(k)=0.
	uuav(k)=0.
	vvav(k)=0.
	wwav(k)=0.
        uav(k)=0.
	vav(k)=0.
	wav(k)=0.
        do i=ix1,ix2+1
          uavg(i,k)=0.
	  vavg(i,k)=0.
	  wavg(i,k)=0.
	  pavg(i,k)=0.
	  ppavg(i,k)=0.
	  uuavg(i,k)=0.
	  vvavg(i,k)=0.
	  wwavg(i,k)=0.
	  uwavg(i,k)=0.
	  uuuavg(i,k)=0.
	  wwuavg(i,k)=0.
c	  omxavg(i,k)=0.
c	  omyavg(i,k)=0.
c	  omzavg(i,k)=0.
#ifdef LES
	  tvavg(i,k)=0.
#endif
	enddo
      enddo
c	initialize reynolds stress variables
      rsnames(0)='rsuu'
      rsnames(1)='rsvv'
      rsnames(2)='rsww'
      rsnames(3)='rsuw'
      
      do j=1,33
        do k=1,nz
	  do i=1,nx
	    rsblk_2(i,k,j)=0.
	  enddo    
        enddo
      enddo	
#ifdef LES
      do j=1,25
        do k=1,nz
	  do i=1,nx
	    rsblk_4(i,k,j)=0.
	  enddo    
        enddo
      enddo	
#endif      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute uinfinity -- all processors, then gather on all CPUs
      subroutine uinfinity(wo,xu,zwg,uinft,uinftg,uactg,samples,real1d)      
      implicit none
      include 'headers/dimension.h'
      include 'headers/common.h'
      include 'mpif.h'
      
      real uact,uactg,uinft(nz),uinftg(nzg),uint,fac
      integer iuinf,type,samples,ierr,status(MPI_STATUS_SIZE),real1d
      real wo(nx,ny,nz),xu(nx),zwg(nzg)     
      integer i,j,k
      
#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

      uact=0.
      iuinf=ix2-3    
      do k=kz1,kz2
         uint=0. 
         do i=iuinf,ix2
          do j=jy1+1,jy2
             uint=uint+wo(i,j,k)/ap(i)
          enddo
         enddo
c       
         uint=uint/(xu(ix2)-xu(iuinf-1))/(jy2-jy1)
         uinft(k)=uinft(k)+(uint-uinft(k))/samples
         uact=uact+uinft(k)/cw(k) ! for integration over x
      enddo
c     gather
      call MPI_GATHER(uinft(2),1,real1d,uinftg(2),1,real1d,0,MPI_COMM_HYDRO,ierr)
c     include first point (optional)
      if (myrank.eq.0) uinftg(1)=uinft(1)
      call MPI_BCAST(uinftg(1),nzg,type,0,MPI_COMM_HYDRO,ierr)	     
      call MPI_ALLREDUCE(uact,uactg,1,TYPE,MPI_SUM,MPI_COMM_HYDRO,ierr)
c     get averaged uinfinity
      fac=1./(1.5*zwg(nzg-1)-0.5*zwg(nzg-2)-(1.5*zwg(1)-0.5*zwg(2)))
      uactg=uactg*fac
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stats_old(uo,vo,wo,zwg,uav,vav,wav,uuav,vvav,wwav,guav,gvav,gwav,
     &                      guuav,gvvav,gwwav,samples,statnum,imid,imidf,index,real1d)
c     compute statistics in the free-stream      
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'mpif.h'
      
      real u,v,w,uu,vv,ww
      integer np,imid,imidf,samples,statnum,ierr,real1d
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),zwg(nzg)
      
      real uav(nz),vav(nz),wav(nz) ! real stat arrays
      real uuav(nz),vvav(nz),wwav(nz) ! real stat arrays
      real guav(nzg),gvav(nzg),gwav(nzg) ! real stat arrays
      real guuav(nzg),gvvav(nzg),gwwav(nzg) ! real stat arrays
      
      integer i,j,k
      character*3 index(0:999)
      
      do k=kz1,kz2
        u=0.
        v=0.
        w=0.
        uu=0.
        vv=0.
        ww=0.
        do i=imid,imidf
         do j=jy1+1,jy2
c     spatial averaging	  
             uu=uu+(uo(i,j,k))**2 
             u=u+uo(i,j,k) 
             vv=vv+(vo(i,j,k))**2 
             v=v+vo(i,j,k) 
             ww=ww+(wo(i,j,k))**2 
             w=w+wo(i,j,k) 
          enddo
        enddo
        np=(jy2-jy1)*(imidf-imid+1)
        uu=uu/np
        u=u/np
        v=v/np
        vv=vv/np
        w=w/np
        ww=ww/np
c     
        uav(k)=uav(k)+(u-uav(k))/samples
        vav(k)=vav(k)+(v-vav(k))/samples
        wav(k)=wav(k)+(w-wav(k))/samples
        uuav(k)=uuav(k)+(uu-uuav(k))/samples
        vvav(k)=vvav(k)+(vv-vvav(k))/samples
        wwav(k)=wwav(k)+(ww-wwav(k))/samples
c     
      enddo 
      call MPI_GATHER(uav(1),1,real1d,guav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
      call MPI_GATHER(vav(1),1,real1d,gvav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
      call MPI_GATHER(wav(1),1,real1d,gwav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
      call MPI_GATHER(uuav(1),1,real1d,guuav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
      call MPI_GATHER(vvav(1),1,real1d,gvvav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
      call MPI_GATHER(wwav(1),1,real1d,gwwav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c	       
      if (myrank.eq.mysize-1) then 
c     include last point
         guav(nzg-1)=uav(kz2)
         gvav(nzg-1)=vav(kz2)
         gwav(nzg-1)=wav(kz2)
         guuav(nzg-1)=uuav(kz2)
         gvvav(nzg-1)=vvav(kz2)
         gwwav(nzg-1)=wwav(kz2)
c	       
         open(21,file='h'//index(statnum)//'.gnu',status='unknown',form='formatted')
         write(21,'(4 F20.15)') (zwg(k)-x_ref, guuav(k)-guav(k)**2, gvvav(k)-gvav(k)**2, 
     &                         	       gwwav(k)-gwav(k)**2,k=1,nzg-1)
         close(21)
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stats_new(uo,vo,wo,p,tv,xcf,xu,zcg,zwg,uavg,vavg,wavg,pavg,uuavg,vvavg,wwavg,uwavg,ppavg,
     &                     uuuavg, wwuavg,utaug,tauwg,
     &                     uinft,uinftg,dstg,thetag,uactg,
     &                     winstg,cfg,kpg,ubulkg,samples,statnum,real1d,stat2d,index,icycle)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'mpif.h'
      include 'headers/stat.h'
      
      real u,v,w,uu,vv,ww,pp,ppp,uw,uuu,wwu
      integer np,stat2d,real1d,samples,statnum
      
      real uinft(nz),ubulk(nz),dst(nz),theta(nz),winst(nz)
      
      real uinftg(nzg),uactg
      real tauwg(nzg),utaug(nzg),ubulkg(nzg),dstg(nzg),
     &     thetag(nzg),cfg(nzg),kpg(nzg),winstg(nzg)

      real uavg(nx,nz),vavg(nx,nz),wavg(nx,nz),uuavg(nx,nz),vvavg(nx,nz),wwavg(nx,nz),uwavg(nx,nz),
     &         pavg(nx,nz),ppavg(nx,nz), uuuavg(nx,nz), wwuavg(nx,nz)
      real uavgg(nx,nzg),vavgg(nx,nzg),wavgg(nx,0:nzg),uuavgg(nx,nzg),vvavgg(nx,nzg),
     &           wwavgg(nx,nzg),uwavgg(nx,nzg),pavgg(nx,nzg),ppavgg(nx,nzg), uuuavgg(nx,nzg), wwuavgg(nx,nzg)
      character*3 index(0:999)
      
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz),tv(nx,ny,nz)
      real xcf(nx),xu(nx),zcg(nzg),zwg(nzg)
      
      integer icycle,ierr,i,j,k

#ifdef LES
      real tvavg(nx,nz),tvavgg(nx,nzg),tvi
#endif
       
      np=jy2-jy1
      do k=kz1,kz2
        do i=ix1,ix2
           u=0.
           v=0.
           w=0.
           pp=0.
           uu=0.
           ww=0.
           vv=0.
           uw=0.
           ppp=0.
	   uuu=0.
	   wwu=0.
#ifdef LES	       
           tvi=0.
#endif
           do j=jy1+1,jy2
              u=u+uo(i,j,k)
              v=v+vo(i,j,k)
              w=w+wo(i,j,k)
	      pp=pp+p(i,j,k)
              uu=uu+uo(i,j,k)**2
              vv=vv+vo(i,j,k)**2
              ww=ww+wo(i,j,k)**2
              uw=uw+0.25*(uo(i+1,j,k+1)+uo(i,j,k+1))*(wo(i+1,j,k+1)+wo(i+1,j,k))! at k+1 point
              ppp=ppp+p(i,j,k)**2		  
              uuu=uuu+uo(i,j,k)**3
              wwu=wwu+0.125*(uo(i+1,j,k+1)+uo(i,j,k+1))*(wo(i+1,j,k+1)+wo(i+1,j,k))**2! at k+1 point
#ifdef LES
              tvi=tvi+tv(i,j,k+1)
#endif
           enddo  
           u=u/np
           v=v/np
           w=w/np
	   pp=pp/np
           uu=uu/np
           vv=vv/np
           ww=ww/np
           uw=uw/np
	   ppp=ppp/np
	   uuu=uuu/np
	   wwu=wwu/np
           uavg(i,k)=uavg(i,k)+(u-uavg(i,k))/samples
           vavg(i,k)=vavg(i,k)+(v-vavg(i,k))/samples
           wavg(i,k)=wavg(i,k)+(w-wavg(i,k))/samples
	   pavg(i,k)=pavg(i,k)+(pp-pavg(i,k))/samples
           uuavg(i,k)=uuavg(i,k)+(uu-uuavg(i,k))/samples
           vvavg(i,k)=vvavg(i,k)+(vv-vvavg(i,k))/samples
           wwavg(i,k)=wwavg(i,k)+(ww-wwavg(i,k))/samples
           uwavg(i+1,k+1)=uwavg(i+1,k+1)+(uw-uwavg(i+1,k+1))/samples
      	   ppavg(i,k)=ppavg(i,k)+(ppp-ppavg(i,k))/samples
      	   uuuavg(i,k)=uuuavg(i,k)+(uuu-uuuavg(i,k))/samples
      	   wwuavg(i+1,k+1)=wwuavg(i+1,k+1)+(wwu-wwuavg(i+1,k+1))/samples
#ifdef LES
	   tvi=tvi/np
	   tvavg(i,k+1)=tvavg(i,k+1)+(tvi-tvavg(i,k+1))/samples
#endif
        enddo
            
c     compute displacement thickness, momentum thickness, viscous units, skin friction coeff.
c     
        ubulk(k)=0.
        theta(k)=0.
        dst(k)=0.
        winst(k)=0.
        do i=ix1+1,ix2
           ubulk(k)=ubulk(k)+wavg(i,k)/ap(i)
           dst(k)=dst(k)+(uinft(k)-wavg(i,k))/ap(i)
           theta(k)=theta(k)+wavg(i,k)*(uinft(k)-wavg(i,k))/ap(i)
        enddo
        ubulk(k)=ubulk(k)/(xu(ix2)-xu(ix1)) 
        dst(k)=dst(k)/uinft(k)
        theta(k)=theta(k)/(uinft(k)**2)
        do j=jy1+1,jy2
           winst(k)=winst(k)+wo(2,j,k)
	enddo
	winst(k)=winst(k)/np	
      enddo
c     
c     SUM CONTRIBUTIONS FROM DIFFERENT PROCESSORS
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': reducing stats'
      call MPI_GATHER(ubulk(1),1,real1d,ubulkg(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(dst(1),1,real1d,dstg(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(theta(1),1,real1d,thetag(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
c        for instantaneous cf curve
      call MPI_GATHER(winst(1),1,real1d,winstg(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
c        
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': gathering stats'
c     gather data on the last processor for load balancing
      call MPI_GATHER(uavg(1,1),1,stat2d,uavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(vavg(1,1),1,stat2d,vavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(wavg(1,1),1,stat2d,wavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(pavg(1,1),1,stat2d,pavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(uuavg(1,1),1,stat2d,uuavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(vvavg(1,1),1,stat2d,vvavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(wwavg(1,1),1,stat2d,wwavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(uwavg(1,2),1,stat2d,uwavgg(1,2),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr) ! from 2nd point
      call MPI_GATHER(ppavg(1,1),1,stat2d,ppavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(uuuavg(1,1),1,stat2d,uuuavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
      call MPI_GATHER(wwuavg(1,2),1,stat2d,wwuavgg(1,2),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr) ! from 2nd point
#ifdef LES
      call MPI_GATHER(tvavg(1,2),1,stat2d,tvavgg(1,2),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
#endif
c       
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
c     
c *************     -- only if one needs to know the stats at k=kz2, otherwise comment
c     include the last point (k=kz2)
      if (myrank.eq.mysize-1) then
          ubulkg(nzg-1)=ubulk(kz2)
          dstg(nzg-1)  =dst(kz2)
          thetag(nzg-1)=theta(kz2)
	  winstg(nzg-1)=winst(kz2)
          do i=ix1+1,ix2 
            uavgg(i,nzg-1) =uavg(i,kz2)
            vavgg(i,nzg-1) =vavg(i,kz2)
            wavgg(i,nzg-1) =wavg(i,kz2)
            uuavgg(i,nzg-1)=uuavg(i,kz2)
            uuuavgg(i,nzg-1)=uuuavg(i,kz2)
            vvavgg(i,nzg-1)=vvavg(i,kz2)
            wwavgg(i,nzg-1)=wwavg(i,kz2)
            uwavgg(i,1)=0. ! because stats here cannot be computed (we start at k=2)
            wwuavgg(i,1)=0. ! because stats here cannot be computed (we start at k=2)
            ppavgg(i,nzg-1)=ppavg(i,kz2)
            wavgg(i,0)=0. ! for consistency with output wavgg(i,0) has to be defined
#ifdef LES
            tvavgg(i,1)=0. ! tv is unknown at first point
#endif
            pavgg (i,nzg-1)=pavg (i,kz2)
       	 enddo
      endif	 
c ************
      if (myrank .eq.mysize-1) then
c      open(23,file='stress'//index(statnum)//'.gnu',status='unknown',form='formatted')
       do k=kz1,nzg-1 
        tauwg(k)=ru1*wavgg(2,k)*au(1)*2. !second order
        utaug(k)=sqrt(tauwg(k))
        cfg(k)=tauwg(k)*2./uactg**2
ccc 	for instantaneous cf curve
        winstg(k)=ru1*winstg(k)*au(1)*2. !second order tau_wall
        winstg(k)=winstg(k)*2./uactg**2  !uactg is NOT instant, but should be ok
ccc	        
        kpg(k)=-ru1/(0.5*(uinftg(k)+uinftg(k+1)))**2*(uinftg(k+1)-uinftg(k))/(zwg(k+1)-zwg(k))
c     
      goto 999   ! obsolescent output
          write(23,*) char(10),char(10)
          write(23,*) k, zwg(k)
          do i=ix1+1,ix2
                 write(23,'(19 G20.12)') xu(i)*utaug(k)/ru1,xu(i)/dstg(k),0.5*(xu(i)+xu(i-1))*utaug(k)/ru1,
     &                 uavgg(i,k)/utaug(k),vavgg(i,k)/utaug(k),wavgg(i,k)/utaug(k),
     &                 (uuavgg(i,k)-uavgg(i,k)**2)/tauwg(k),(vvavgg(i,k)-vavgg(i,k)**2)/tauwg(k),
     &                 (wwavgg(i,k)-wavgg(i,k)**2)/tauwg(k),
     &                 (uwavgg(i,k)-uavgg(i,k)*wavgg(i,k))/tauwg(k),  
c                          normalize by macro scales
     &                 (uavgg(i,k)/uactg),(vavgg(i,k)/uactg),(wavgg(i,k)/uactg),
     &                 (uuavgg(i,k)-uavgg(i,k)**2)/uactg**2,(vvavgg(i,k)-vavgg(i,k)**2)/uactg**2,
     &                 (wwavgg(i,k)-wavgg(i,k)**2)/uactg**2,
     &                 (uwavgg(i,k)-0.25*(uavgg(i,k)+uavgg(i-1,k))*(wavgg(i,k)+wavgg(i,k-1)))/uactg**2 !data invalid at 1st point
#ifdef LES 
     &                 ,tvavgg(i,k)/ru1
#endif     
          enddo
 999  continue
       enddo !k
c	        	      
c            close(23)
c     write averaged fields (matlab)
       open(62,file='avg'//index(statnum)//'.mat',status='unknown',form='unformatted')
c              SAME OUTPUT -- ROTATE GRAPH	       
       write(62) nx-2,nzg-2 
       write(62) (xcf(i),i=2,ix2)		 
       write(62) (zcg(k),k=2,nzg-1)
       do k=2,nzg-1
        write(62) (
     &            0.5*(uavgg(i,k)+uavgg(i-1,k)),
     &            vavgg(i,k),
     &            0.5*(wavgg(i,k)+wavgg(i,k-1)),
     &            pavgg(i,k),
     &                     !char(10),
     &            i=2,ix2)
       enddo
       close(62)
c     write averaged _stress_ field
       open(62,file='avg_str'//index(statnum)//'.mat',status='unknown',form='unformatted')
c              SAME OUTPUT -- ROTATE GRAPH	       
       write(62) nx-2,nzg-2 
       write(62) (xcf(i),i=2,ix2)		 
       write(62) (zcg(k),k=2,nzg-1)
       do k=2,nzg-1
        write(62) (
     &            0.5*(uuavgg(i,k)+uuavgg(i-1,k)- (uavgg(i,k)**2+uavgg(i-1,k)**2)*subtract_meanvel ),
     &            vvavgg(i,k)- (vavgg(i,k)**2)*subtract_meanvel,
     &            0.5*(wwavgg(i,k)+wwavgg(i,k-1)- (wavgg(i,k)**2+wavgg(i,k-1)**2)*subtract_meanvel ),
     &            uwavgg(i,k)-0.25*(uavgg(i,k)+uavgg(i-1,k))*(wavgg(i,k)+wavgg(i,k-1))*subtract_meanvel,
     &            ppavgg(i,k)-(pavgg(i,k)**2)*subtract_meanvel,
     &                    !char(10),
     &            i=2,ix2)
       enddo
       close(62)
c     write _averaged_ vorticity field
       open(62,file='vortavg'//index(statnum)//'.mat',status='unknown',form='unformatted')
c              SAME OUTPUT -- ROTATE GRAPH	       
       write(62) nx-4,nzg-3 !x-dir lacks 1st and last pt; z-dir lacks last pt 
       write(62) (xcf(i),i=3,ix2-1)		 
       write(62) (zcg(k),k=2,nzg-2)
       do k=2,nzg-2
        write(62) (
     &            0.5*((vavgg(i,k)-vavgg(i,k-1))*cwg(k-1)
     &                +(vavgg(i,k+1)-vavgg(i,k))*cwg(k)),  !omx
     &          -0.25*((uavgg(i,k)-uavgg(i,k-1)+uavgg(i-1,k)-uavgg(i-1,k-1))*cwg(k-1)+
     &                 (uavgg(i,k+1)-uavgg(i,k)+uavgg(i-1,k+1)-uavgg(i-1,k))*cwg(k)-
     &	               (wavgg(i,k)-wavgg(i-1,k)+wavgg(i,k-1)-wavgg(i-1,k-1))*au(i-1)-
     &		       (wavgg(i+1,k)-wavgg(i,k)+wavgg(i+1,k-1)-wavgg(i,k-1))*au(i)),  !omy
     &           -0.5*((vavgg(i,k)-vavgg(i-1,k))*au(i-1)+
     &		       (vavgg(i+1,k)-vavgg(i,k))*au(i)),
#ifdef LES
     &                  tvavgg(i,k),
#endif
     &             i=3,ix2-1)
       enddo
       close(62)
c     write averaged triple correlation field (v^3+u^2v)
       open(62,file='avg_3_'//index(statnum)//'.mat',status='unknown',form='unformatted')
c              SAME OUTPUT -- ROTATE GRAPH	       
       write(62) nx-2,nzg-2 
       write(62) (xcf(i),i=2,ix2)		 
       write(62) (zcg(k),k=2,nzg-1)
       do k=2,nzg-1
        write(62) (
     &            0.5*(uuuavgg(i,k)+uuuavgg(i-1,k)
     &                -3.*(uuavgg(i,k)*uavgg(i,k)+uuavgg(i-1,k)*uavgg(i-1,k))
     &                +2.*(uavgg(i,k)**3+uavgg(i-1,k)**3)),
     &                 wwuavgg(i,k)-uwavgg(i,k)*(wavgg(i,k)+wavgg(i,k-1))
     &           -0.25*(uavgg(i,k)+uavgg(i-1,k))*(wwavgg(i,k)+wwavgg(i,k-1))
     &           +0.5*(uavgg(i,k)+uavgg(i-1,k))*(wavgg(i,k)**2+wavgg(i,k-1)**2),
     &            i=2,ix2)
       enddo
       close(62)
      endif	!myrank
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        compute stress budget data
      subroutine rstress(rsblk_1,rsblkg_1,rsblk_2,
#ifdef LES      
     &                   rsblk_3,rsblk_4,rsblkg_3,
#endif     
     &                   uo,vo,wo,p,tv,utaug,tauwg,samples,statnum,rsnames,xcf,zcg,stat2d,index)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'mpif.h'

      integer i,j,k,n,ierr
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz),tv(nx,ny,ny)
      real xcf(nx),zcg(nzg)
      real rsblk_1(nx,nz,ny)
      real rsblk_2(nx,0:nz+1,33)
      real rsblkg_1(nx,nzg,7) 
#ifdef LES
      real rsblk_3 (nx,nz,ny)     !temporary LES block
      real rsblk_4 (nx,0:nz+1,25) !permanent LES block
      real rsblkg_3(nx,nzg,4)     !temporary global LES block
#endif
      integer stat2d, statnum, samples
      real tauwg(nzg), utaug(nzg)
      character*4 rsnames(0:3)
      character*3 index(0:999)


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
      end
