      integer nspec,specfreq,skip
      character*80 specname
      parameter (nspec=3,specfreq=5,skip=10) !spectra to be computed at nspec locations; 
                                             !data sampled once in every specfreq timesteps at every "skip"'th z-point
      real specpos(nspec)
      data specname/'spec/t3bh.s'/
      data specpos/40,50,65/


!cccc VARIABLES FOR STREAMWISE SPECTRA
      include 'headers/spec.h'
      integer spec_cpu(nspec), spec_k(nspec)
!cccccccccccccccccccccccccccccccccccccccccccccccccc      


!
      if (spec) call spec_setup(spec_cpu,spec_k,xu,zwg,index)
!


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine spec_setup(spec_cpu,spec_k,xu,zwg,ind)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/spec.h'
      
      real xu(nx), zwg(nzg)
      integer i,j,k,spec_cpu(nspec), spec_k(nspec)
      character*3 ind(0:999)
      integer nlen ! length of specname before spaces
c
      nlen=index(specname,' ')-1
      if (nlen.le.0) nlen=80 

      if (myrank.eq.0) then 
      	   if (lcfl) write(0,*) 'WARNING: spectra sampled with cfl ON!' 
      endif
c         
c        specpos specified in spec.h
c	 
         do i=1,nspec !same code as for 'save'
      	 
          spec_k(i)=1
	  spec_cpu(i)=0 
          k=1
          do while (zwg(k).lt.(specpos(i)))
c           write(0,*) i,specpos(i),k,zwg(k)    
	   k=k+1
          enddo
	  if (k.gt.1) then 
	    spec_cpu(i)=(k-2)/(nz-2)
	    spec_k(i)=mod(k-2,nz-2)+1
          endif
	  if (spec_cpu(i).gt.mysize-1) then 
	    spec_cpu(i)=spec_cpu(i)-1
	    spec_k(i)=nz
	  endif  
          if (myrank.eq.spec_cpu(i)) then 
      	     open(600+i,file=specname(1:nlen)//'.'//ind(int(specpos(i))),form='unformatted',status='unknown')
             write(600+i) specpos(i),nx,(xu(k),k=1,nx,skip)
       	  endif      
         enddo
      return
      end
ccccccccccccccccccccccccccccccccccc
      subroutine spec_write(uo,vo,wo,p,spec_cpu,spec_k,tlevel,j)
      implicit none
      include 'headers/common.h'
      include 'headers/spec.h'
      include 'headers/dimension.h'
      
      integer i,j,n
      integer spec_k(nspec), spec_cpu(nspec)
      real uo(nx,ny,nz), vo(nx,ny,nz), wo(nx,ny,nz), p(nx,ny,nz), tlevel
      
        do n=1,nspec
	 if (myrank.eq.spec_cpu(n)) then
	   write(600+n) tlevel,(wo(i,j,spec_k(n)),i=2,nx-1,skip),
     &	                       (vo(i,j,spec_k(n)),i=2,nx-1,skip), 
     &                         (wo(i,j,spec_k(n)),i=2,nx-1,skip),
     &	                       (p (i,j,spec_k(n)),i=2,nx-1,skip)
	 endif
	enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccc
      subroutine spec_close(spec_cpu)
      implicit none
      include 'headers/common.h'
      include 'headers/spec.h'
      
      integer i,spec_cpu(nspec)
             
      do i=1,nspec
       if (myrank.eq.spec_cpu(i)) close(600+i)
      enddo  

      return
      end
      


      