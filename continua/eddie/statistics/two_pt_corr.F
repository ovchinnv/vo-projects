!     file for two-point correlations

      integer nxl,nzl
      integer zfreq,xfreq
      parameter (xfreq=5,zfreq=10)
      parameter (nxl=(nx-2+xfreq-1)/xfreq, nzl=(nz-2+zfreq-1)/zfreq) !   
 
      character*80 corrname
      data corrname/'corr/fstle2h.two_pt'/








!ccccccccc variables for two-pt correlations ccccccccccc           
      integer nz_act,nx_act
      real xn(nxl), zn((nzg-2+zfreq-1)/zfreq)        
      integer ks(0:ncpu-1),ke(0:ncpu-1),npt(0:ncpu-1) ! these arrays are defined on all processors
      real um2pt(nxl,nzl), vm2pt(nxl,nzl), wm2pt(nxl,nzl)
      real u2pt(nxl,ny,nzl),v2pt(nxl,ny,nzl),w2pt(nxl,ny,nzl)
      real u2tmp(nxl,ny,nzl),v2tmp(nxl,ny,nzl),w2tmp(nxl,ny,nzl) ! temporary arrays







ccccccccccccccccc
      subroutine two_pt_corr_setup(nx_act,nz_act,xn,zn,ks,ke,npt,xcf,zcg)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/two_pt.h'
c
      integer nz_act,nx_act
      real xn(nxl), zn((nzg-2+zfreq-1)/zfreq)        
      integer ks(0:ncpu-1),ke(0:ncpu-1),npt(0:ncpu-1) ! these arrays are defined on all processors
      
      real zcg(nzg), xcf(nx)
      
      integer ii,i,kk,k,j
c      compute 2-pt correlation grids (center)      
       ii=0
       do i=2,nx-1,xfreq
        ii=ii+1
        xn(ii)=xcf(i)
       enddo
       nx_act=ii
c
       kk=0
       do k=2,nzg-1,zfreq
        kk=kk+1
        zn(kk)=zcg(k)
       enddo
       
c  compute local z-indices
       nz_act=0
       do j=0,ncpu-1
        ks(j)=(nz-2)*j+2
	ke(j)=ks(j)+nz-3
	ks(j)=ks(j)+mod(zfreq-mod(ks(j)-2,zfreq),zfreq)! still global at this point
	ke(j)=ke(j)-mod(ke(j)-2,zfreq)      ! still global
	npt(j)=(ke(j)-ks(j))/zfreq+1        ! number of entries
        nz_act=nz_act+npt(j)
cccccc debug
c      if (myrank.eq.0) then 
c       write(35,'(5I)') j,ks(j),ke(j),ks(j)-(nz-2)*j,ke(j)-(nz-2)*j
c      endif
c      
cccccc debug
c       convert to local form
        ks(j)=ks(j)-(nz-2)*j	
	ke(j)=ke(j)-(nz-2)*j
       enddo
c      close(35)
   
c      do j=1,ncpu-1
c       write(36,'(6I)') j,nxl,nx_act,nzl,nz_act,npt(j)
c      enddo
c      close(36) 
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine span_corr(u,v,w,u2,v2,w2,um,vm,wm,ks,ke,samples)
c     computes spanwise correlations at every zfreq'th z-point and every xfreq'th x-point
c     VO 5/04 Please do not edit, it's more complicated than seems
      implicit none
      include 'headers/dimension.h'
      include 'headers/two_pt.h'

      integer i,j,k,ii,jj,kk,ks,ke,samples
      logical init
      integer ifaxy(13)
      real trigsy(3*(ny-2)/2+1), work(nxl*ny*nzl)
      
      real u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
      real u2(nxl,ny,nzl), v2(nxl,ny,nzl), w2(nxl,ny,nzl)
      real um(nxl,nzl),vm(nxl,nzl),wm(nxl,nzl)
      real umean,vmean,wmean
      
c     dummy array      
      real uh(ny,nxl,nzl), vh(ny,nxl,nzl), wh(ny,nxl,nzl)
      
      save init,trigsy,ifaxy
      data init /.true./
      
      if (init) then
c      initialize arrays
       do k=1,nzl
         do i=1,nxl
	  um(i,k)=0.
	  vm(i,k)=0.
	  wm(i,k)=0.
	  do j=1,ny
	   u2(i,j,k)=0.
	   v2(i,j,k)=0.
	   w2(i,j,k)=0.
	  enddo
         enddo	   
       enddo
c      set up FFT routine 1st time set up     
       call fftfax(ny-2,ifaxy,trigsy)
       init=.false.
      endif 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      COMPUTE CORRELATIONS 
c      1 -- copy into temporary array
      kk=0
      do k=ks,ke,zfreq ! outer index
        kk=kk+1
	ii=0
        do i=2,nx-1,xfreq ! inner index
	  ii=ii+1 
          umean=0
	  vmean=0
	  wmean=0
	  do j=2,ny-1
	    uh(j,ii,kk)=0.5*(u(i,j,k)+u(i-1,j,k))
	    vh(j,ii,kk)=v(i,j,k)
	    wh(j,ii,kk)=0.5*(w(i,j,k)+w(i,j,k-1))
	    umean=umean+uh(j,ii,kk) ! compute mean
	    vmean=vmean+vh(j,ii,kk) ! compute mean
	    wmean=wmean+wh(j,ii,kk) ! compute mean
	  enddo
	  umean=umean/(ny-2)
	  vmean=vmean/(ny-2)
	  wmean=wmean/(ny-2)
c
	  uh(ny,ii,kk)=uh(2,ii,kk)
	  vh(ny,ii,kk)=vh(2,ii,kk)
	  wh(ny,ii,kk)=wh(2,ii,kk)
	  uh(1,ii,kk)=uh(ny-1,ii,kk)
	  vh(1,ii,kk)=vh(ny-1,ii,kk)
	  wh(1,ii,kk)=wh(ny-1,ii,kk)

	  um(ii,kk)=um(ii,kk)+(umean-um(ii,kk))/samples !update mean
	  vm(ii,kk)=vm(ii,kk)+(vmean-vm(ii,kk))/samples 
	  wm(ii,kk)=wm(ii,kk)+(wmean-wm(ii,kk))/samples 
	enddo
      enddo
       	  
c      2 -- transform	      
      call fft99(uh,work,trigsy,ifaxy,1,ny,ny-2,nzl*nxl,-1)
      call fft99(vh,work,trigsy,ifaxy,1,ny,ny-2,nzl*nxl,-1)
      call fft99(wh,work,trigsy,ifaxy,1,ny,ny-2,nzl*nxl,-1)
c      3 -- compute correlation
      do k=1,nzl
       do i=1,nxl
        do j=1,ny/2
	 jj=2*j-1
	 uh(jj,i,k)=uh(jj,i,k)**2+uh(jj+1,i,k)**2
	 uh(jj+1,i,k)=0.
	 vh(jj,i,k)=vh(jj,i,k)**2+vh(jj+1,i,k)**2
	 vh(jj+1,i,k)=0.
	 wh(jj,i,k)=wh(jj,i,k)**2+wh(jj+1,i,k)**2
	 wh(jj+1,i,k)=0.
	enddo
       enddo
      enddo        	 
c      3 -- inverse transform
      call fft99(uh,work,trigsy,ifaxy,1,ny,ny-2,nzl*nxl,1)
      call fft99(vh,work,trigsy,ifaxy,1,ny,ny-2,nzl*nxl,1)
      call fft99(wh,work,trigsy,ifaxy,1,ny,ny-2,nzl*nxl,1)
c      4 -- copy into permanent array
      do k=1,nzl
       do i=1,nxl
        do j=1,ny
         u2(i,j,k)=u2(i,j,k)+(uh(j,i,k)-u2(i,j,k))/samples
         v2(i,j,k)=v2(i,j,k)+(vh(j,i,k)-v2(i,j,k))/samples
         w2(i,j,k)=w2(i,j,k)+(wh(j,i,k)-w2(i,j,k))/samples
	enddo
       enddo
      enddo        	 
      return
      end       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_corr(u2pt,v2pt,w2pt,um2pt,vm2pt,wm2pt,nx_act,nz_act,xn,ycf,zn,npt,type)
      implicit none
      include 'headers/dimension.h'
      include 'headers/two_pt.h'
      include 'headers/common.h'
      include 'mpif.h'
      
      integer i,j,k,ierr,n,nn
      integer nx_act, nz_act, npt(0:ncpu-1)      
      real xn(nxl), zn((nzg-2+zfreq-1)/zfreq),ycf(ny) 
      integer type, status(MPI_STATUS_SIZE)
      real um2pt(nxl,nzl), vm2pt(nxl,nzl), wm2pt(nxl,nzl)
      real u2pt(nxl,ny,nzl),v2pt(nxl,ny,nzl),w2pt(nxl,ny,nzl)
      real u2tmp(nxl,ny,nzl),v2tmp(nxl,ny,nzl),w2tmp(nxl,ny,nzl)
      character*80 fname

      integer nlen
      nlen=index(corrname,' ')-1
      if(nlen.le.0) nlen=80
      fname=corrname(1:nlen)
    
c          subtract mean from correlations	    
      do k=1,nzl
       do i=1,nxl
        do j=1,ny
	 u2tmp(i,j,k)=u2pt(i,j,k)-um2pt(i,k)**2
	 v2tmp(i,j,k)=v2pt(i,j,k)-vm2pt(i,k)**2
	 w2tmp(i,j,k)=w2pt(i,j,k)-wm2pt(i,k)**2
	enddo
c        normalize to one
c	 u2tmp(i,1,k)=u2tmp(i,1,k)/u2tmp(i,2,k)
c	 v2tmp(i,1,k)=v2tmp(i,1,k)/v2tmp(i,2,k)
c	 w2tmp(i,1,k)=w2tmp(i,1,k)/w2tmp(i,2,k)
c        do j=ny,2,-1
c	 u2tmp(i,j,k)=u2tmp(i,j,k)/u2tmp(i,2,k)
c	 v2tmp(i,j,k)=v2tmp(i,j,k)/v2tmp(i,2,k)
c	 w2tmp(i,j,k)=w2tmp(i,j,k)/w2tmp(i,2,k)
c	enddo
       enddo
      enddo 	 
c
      if (myrank.eq.0) then 
           open(62,file=corrname(1:nlen),status='unknown',form='unformatted')
           write(62) nx_act,ny-2,nz_act 
           write(62) (xn(i),i=1,nx_act)
           write(62) (ycf(j),j=2,ny-1)		 
           write(62) (zn(k),k=1,nz_act)
           do k=1,npt(myrank)
              write(62) ((
     &                 u2tmp(i,j,k),
     &                 v2tmp(i,j,k),
     &                 w2tmp(i,j,k),
     &                 j=2,ny-1),i=1,nx_act)
           enddo
           do n=1,mysize-1
            nn=n*10
              call MPI_Recv(u2tmp,nxl*nzl*ny,type,n,nn+0,mpi_comm_hydro,status,ierr)
              call MPI_Recv(v2tmp,nxl*nzl*ny,type,n,nn+1,mpi_comm_hydro,status,ierr)
              call MPI_Recv(w2tmp,nxl*nzl*ny,type,n,nn+2,mpi_comm_hydro,status,ierr)
     
              do k=1,npt(n)
                write(62) ((
     &                 u2tmp(i,j,k),
     &                 v2tmp(i,j,k),
     &                 w2tmp(i,j,k),
     &                 j=2,ny-1),i=1,nx_act)
              enddo
           enddo		           	       
           close(62)
      else
c         now send	    
       nn=myrank*10
       call MPI_Send(u2tmp,nxl*nzl*ny,type,0,nn,  MPI_COMM_HYDRO,ierr)
       call MPI_Send(v2tmp,nxl*nzl*ny,type,0,nn+1,MPI_COMM_HYDRO,ierr)
       call MPI_Send(w2tmp,nxl*nzl*ny,type,0,nn+2,MPI_COMM_HYDRO,ierr)
      endif 
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      