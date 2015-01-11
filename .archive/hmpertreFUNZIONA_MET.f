      PROGRAM HMPERTRE
      implicit none

      integer L,nloop,iloop,i,j,k,imet,lenw,lens,ier,inc,lenc,n
      parameter(L=4*4096)
      parameter(n=2*L,lenw=2*n,lens=2*n*int(log(real(n)))+4)
      
      double precision fmesh,D,U
      double precision dt,pi,beta,ex
      double precision xmu,w,sig,zerop,t0,zeta,zetapp,w0
      double precision V,t,ed0,ome(2*L),ep0,xnd,xnp,xntot,tum,tuma
      double precision dome,mues0,reg,img,dsig
      double precision zzerop,zzerom,dsigma
      double precision tau,sq,sq0,ggr,ggi,gg
      double precision work(lenw),wsave(lens)

      double complex xi,one,zero,D2,iome
      double complex root,a,b,iw
      double complex fg0(2*L),fg(2*L),gc(2*L)
      double complex xg0t(-L:L),g0t(-L:L)    
      double complex sigma(2*L),xs(-L:L)
      double complex xsigmat(2*L),sf,gf

!     improvements WIP
      double complex sigmat(-L:L)
      complex*16 green0,green,sqroot,sqroot0


      open(99,file='inputIPT.in',status='old')
      read(99,*)fmesh,D,U
      read(99,*)nloop
      close(99)

      one=(1.d0,0.d0)
      xi=(0.d0,1.d0)
      zero=(0.d0,0.d0)
      D2=D*one
      pi=datan(1.d0)*4.d0
      dt=2.d0*pi/fmesh          
      dt=dt/dfloat(n)

!     initialize some functions
      do i=1,2*L
         fg0(i)  =zero
         fg(i)   =zero
         sigma(i)=zero
      enddo

      
      do i=1,n
         if(i.le.L)then
            ome(i)=dfloat(i-1)  !>0 [0,wmax]
         else
            ome(i)=dfloat(i-n-1) !<0 [-wmax,0-eta]
         endif
         ome(i)=ome(i)*fmesh
      enddo
      imet=1
      if(U.gt.3.3d0)imet=0
      

!     Starts DMFT-loop
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do iloop=1,nloop

         do i=1,L  
            w=ome(i)            !>0
            sig=1.d0
            sf=sigma(i)
            if(iloop.eq.1)then  !Guess G
               fg0(i)=one/(w+xi*(dfloat(imet)+1.d-3)/2.d0) !2.d0=pinning  
               if(abs(w).lt.1.d-9)fg0(i)=0.d0*one
            else
               iw=w-sf
               sqroot=cdsqrt(iw**2-D2)
               gf=2.d0/(iw+sqroot)
               gf=one/gf+sf
               reg=real(gf)
               img=dimag(gf)
               dsig=img/(reg**2+img**2)
               if(dsig.lt.0.d0)sig=-1.d0
               fg0(i)=2.d0/(w+sf+sig*sqroot)
c     fg0(i)=one/fg(i)+sigma(i)
c     fg0(i)=one/fg0(i)
               fg0(1)=zero
            endif
         enddo


!We need CAUSAL Green's function (G''>0 for w<0)
         do i=1,L-1
            fg0(2*L-i)=-fg0(i+1)
         enddo
         fg0(n)=fg0(n-1)


!     
!     FFT w-->t
!-------------------
         do i=1,n
            if(iloop.ge.nloop)then
               write(80,*)ome(i),dimag(fg0(i)),real(fg0(i))
            endif
            gc(i)=fg0(i)
         enddo

         call four1(gc,n,1)
!     starts manipulations of the arrays:
!     1) [1,2*L=n]---> [-L,L]
         do i=1,n
            xg0t(i-L-1)=fmesh/2.d0/pi*gc(i)
         enddo
!     2) g[0,L-1]<--- x[-L,-1]
         do i=-L,-1
            g0t(i+L)=xg0t(i)
         enddo
!     3) g[-L,-1]<--- x[0,L]
         do i=0,L-1
            g0t(i-L)=xg0t(i)   
         enddo         

!     get Sigma: Impurity Solver
         do i=-L+1,L-1
            sigmat(i)=-(U**2)*(g0t(i)**2)*g0t(-i)
            if(iloop.ge.nloop)then
               write(81,*)dfloat(i)*dt,dimag(g0t(i))
               write(82,*)dfloat(i)*dt,dimag(sigmat(i))
            endif
         enddo
         sigmat(-L)=-(U**2)*(g0t(L)**2)*g0t(L-1)

         
!     FFT t-->w + starts manipulations
         do i=1,n
            xsigmat(i)=sigmat(i-L-1)
         enddo
         CALL four1(xsigmat,n,1)
         ex=-1.d0
         do i=1,n
            ex=-ex
            sigma(i)=ex*dt*xsigmat(i)
         enddo

         if(iloop.ge.nloop)then
            do i=2,n
               fg(i)=one/fg0(i)-sigma(i)
               fg(i)=one/fg(i)
               write(83,*)ome(i),dimag(sigma(i)),real(sigma(i))
               write(84,*)ome(i),dimag(fg(i))
            enddo
         endif

      enddo                     !here the dmft loops end up
c=======================================================================
      end
      include 'routines.f'
c      include 'fftpack5.f'
