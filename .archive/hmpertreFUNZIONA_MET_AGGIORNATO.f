      PROGRAM HMPERTRE
      implicit none

      integer L,nloop,iloop,i,j,k,imet,lenw,lens,ier,inc,lenc,n
      parameter(L=4*4096,n=2*L)
      
      double precision fmesh,D,U
      double precision dt,pi,beta,ex
      double precision xmu,w,sig,zerop,t0,zeta,zetapp,w0
      double precision V,t,ed0,ome(2*L),ep0,xnd,xnp,xntot,tum,tuma
      double precision dome,mues0,reg,img,dsig
      double precision zzerop,zzerom,dsigma
      double precision tau,sq,sq0,ggr,ggi,gg,xp,xm

      double complex xi,one,zero,D2,iome
      double complex root,a,b,iw
      double complex fg0(2*L),fg(2*L),gc(2*L)
      double complex xg0t(-L:L),g0t(-L:L)    
      double complex sigma(2*L),xs(-L:L)
      double complex xsigmat(2*L),sf,gf,gii,gi,gp,gm

!     improvements WIP
      double complex sigmat(-L:L)
      complex*16 green0,green,sqroot,sqroot0


      open(99,file='inputIPT.in',status='old')
      read(99,*)fmesh,D,U
      read(99,*)nloop,imet
      close(99)

      one=(1.d0,0.d0)
      xi=(0.d0,1.d0)
      zero=(0.d0,0.d0)
      D2=D*one
      pi=datan(1.d0)*4.d0
      dt=2.d0*pi/fmesh          
      dt=dt/dfloat(n)

      print*, '       U = ',U
      print*, '    Mesh = ',fmesh
c     =======================================================================

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

!     Starts DMFT-loop
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do iloop=1,nloop
         print*,'iloop',iloop,'/',nloop

         do i=1,L  
            w=ome(i)            !>0
            sig=1.d0
            sf=sigma(i)
!     First loop: guess
            if(iloop.eq.1)then
               fg0(i)=one/(w+xi*dfloat(imet)/2.d0) !2.d0=pinning  
               if(abs(w).lt.1.d-9)fg0(i)=0.d0*one
!     other DMFT loops
            else
               iw=w-sf
               sqroot=cdsqrt(iw**2-D2)

!     Old marcelo's way
c     if(i.le.5)then   !this is always true
c     print*,'ciao'
c     if(real(sqroot).gt.0.d0.and.imet.eq.0)sig=1.d0
c     else

!     CONTINUITY criteria
               gp=2.d0/(w+sf+sqroot)
               gm=2.d0/(w+sf-sqroot)
               xp=abs(gp+gii-2.d0*gi)**2
               xm=abs(gm+gii-2.d0*gi)**2
               if(xp.gt.xm)sig=-1.d0
!     BRANCH CUT                  
               gf=(iw+sqroot)/2.d0
               if(dimag(gf).lt.0.d0)sig=-1.d0

!     Marcelo's way.
c     gf=2.d0/(iw+sqroot) !try to get G
c     gf=one/gf+sf     !G0
c     reg=real(gf)
c     img=dimag(gf)
c     dsig=img/(reg**2+img**2)
c     if(dsig.lt.0.d0)sig=-1.d0 !If Go''>0 change branch
c     endif
               if(i.le.5)then
                  if(real(sqroot).gt.0.d0.and.imet.eq.0)sig=1.d0
               endif
               fg0(i)=2.d0/(w+sf+sig*sqroot)
               if(i.eq.1) fg0(i)=zero
            endif
            gi=fg0(i)
            gii=fg0(i-1)
         enddo


!     We need CAUSAL Green's function (G''>0 for w<0)         
         fg0(1)=zero
         do i=1,L-1
            fg0(2*L-i)=-fg0(i+1)
         enddo
         fg0(n)=fg0(n-1)
         

!     
!     FFT w-->t
!-------------------
         do i=1,n
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



!     
!     Impurity Solver: Sigma=U^2G0^3
!-----------------------------
!     get Sigma
         do i=-L+1,L-1
            sigmat(i)=-(U**2)*(g0t(i)**2)*g0t(-i)
         enddo
         sigmat(-L)=-(U**2)*(g0t(L)**2)*g0t(L-1)

!     FFT t-->w 
!-------------------
!     starts manipulations
         do i=1,n
            xsigmat(i)=sigmat(i-L-1)
         enddo
         CALL four1(xsigmat,n,1)
         ex=-1.d0
         do i=1,n
            ex=-ex
            sigma(i)=ex*dt*xsigmat(i)
         enddo



         if(iloop.ge.nloop)then !write
            open(10,file='Sigma_t.ipt')
            open(11,file='G_t.ipt')
            do i=-L,L
               write(10,*)dfloat(i)*dt,dimag(sigmat(i)),real(sigmat(i))
               write(11,*)dfloat(i)*dt,dimag(g0t(i)),real(g0t(i))
            enddo
            close(10)
            close(11)

            open(10,file='G_realw.ipt')
            open(11,file='Sigma_realw.ipt')
            do i=L+1,2*L,50
               fg(i)=one/fg0(i)-sigma(i)
               fg(i)=one/fg(i)
               write(10,*)ome(i),dimag(fg(i)),real(fg(i))
               write(11,*)ome(i),dimag(sigma(i)),real(sigma(i))
            enddo
            do i=1,L,50
               fg(i)=one/fg0(i)-sigma(i)
               fg(i)=one/fg(i)
               fg(1)=-fg(L+1)
               write(10,*)ome(i),dimag(fg(i)),real(fg(i))
               write(11,*)ome(i),dimag(sigma(i)),real(sigma(i))
            enddo
            close(10)
            close(11)

            open(10,file='DOS.ipt')
            do i=L+1,2*L,50
               write(10,*)ome(i),dimag(fg(i)))/pi
            enddo
            do i=1,2*L,50
               write(10,*)ome(i),-dimag(fg(i))/pi
            enddo
            close(10)
         endif

      enddo                     !here the dmft loops end up
c=======================================================================
      end
      include 'routines.f'
c     include 'fftpack5.f'
