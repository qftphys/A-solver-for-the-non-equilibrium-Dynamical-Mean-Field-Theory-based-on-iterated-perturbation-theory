      PROGRAM HMPERTRE
      implicit none

      integer L,nloop,iloop,i,n,imet
      parameter(L=4*4096,n=2*L)
      
      double precision fmesh,D,U,zerop
      double precision dt,pi,beta,ex,fermi
      double precision w,sig
      double precision ome(2*L)
      double precision gg,xp,xm,sq

      double complex xi,one,zero,D2
      double complex iw
      double complex fg0(2*L),fg(2*L)
      double complex gc(2*L),g0t(-L:L),xAt(-L:L),xAbt(-L:L)
      double complex sigma(2*L),At(-L:L),Abt(-L:L)
      double complex xsigmat(2*L),sf,gf,gii,gi,gp,gm

      double precision A,Ab,pi2,si,sr,s0
      double complex Adummy(2*L),Abdummy(2*L)
      
      double complex g0tmm(-L:L),g0tpm(-L:L),g0tmp(-L:L)
      double complex g0mm(2*L)
      double complex sigmatmm(-L:L),sigmatpm(-L:L) !,xg0tadv(2*L)
c     double complex xsigmatpm(2*L),xsigmatadv(2*L)
c     double complex sigmamm(2*L),sigmapm(2*L),sigmaadv(2*L)
!     improvements WIP
      double complex sigmat(-L:L)
      complex*16 sqroot


      open(99,file='inputIPT.in',status='old')
      read(99,*)U
      read(99,*)nloop
      read(99,*)beta
      close(99)

      fmesh=0.0005d0
      D=1.d0
      one=(1.d0,0.d0)
      xi=(0.d0,1.d0)
      zero=(0.d0,0.d0)
      D2=D*one
      pi=datan(1.d0)*4.d0
      pi2=2.d0*pi
      dt=2.d0*pi/fmesh          
      dt=dt/dfloat(n)
c     U=2.d0*U      

      print*, '       U = ',U
      print*, '    beta = ',beta
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
            ome(i)=dfloat(i-n) !<0 [-wmax,0-eta]
         endif
         ome(i)=ome(i)*fmesh
      enddo

!     Starts DMFT-loop
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do iloop=1,nloop
         print*,'iloop',iloop,'/',nloop

!     First loop: guess
         if(iloop.eq.1)then
            do i=1,n
               w=ome(i)
               sqroot=cdsqrt(w**2-D2)
               fg0(i)=2.d0/(w-sqroot) !G_0=G for Sigma=0!converge piu' piano
!     marcelo's guess
c     fg0(i)=one/(w+xi*dfloat(imet)/2.d0) !2.d0=pinning  
c     if(abs(w).lt.1.d-9)fg0(i)=0.d0*one
            enddo
         else
            do i=1,L 
               w=ome(i)         !>0
c     SELF CONSISTENCY CONDITION
c---------------------------------------------------------------------
               fg0(i)=w-xi*0.01d0-0.25d0*fg(i) !small xi*eta needed4 INS sol.
               fg0(i)=one/fg0(i)
            enddo         
         endif

!     Get w<0 function
         do i=1,L-1
            fg0(2*L-i)=fg0(i+1)
         enddo
         fg0(n)=fg0(1)
         

!Build the spectral densities
         open(11,file='A_realw.ipt')
         do i=1,n
            w=ome(i)
            A=dimag(fg0(i))/pi
            Ab=A*fermi(w,0.d0,beta)
            Adummy(i)=dcmplx(A,0.d0)
            Abdummy(i)=dcmplx(Ab,0.d0)
            if(iloop.eq.nloop)then
               write(11,*)w,A,Ab
            endif
         enddo
         close(11)
!     
!     FFT w-->t
!-------------------
         call four1(Adummy,n,1)
         call four1(Abdummy,n,1)
!     starts manipulations of the arrays:
!     1) [1,2*L=n]---> [-L,L]
         do i=1,n
            xAt(i-L-1)=fmesh/pi2*Adummy(i)
            xAbt(i-L-1)=fmesh/pi2*Abdummy(i)
         enddo
!     2) g[0,L-1]<--- x[-L,-1]
         do i=-L,-1
            At(i+L)=xAt(i)
            Abt(i+L)=xAbt(i)
         enddo
!     3) g[-L,-1]<--- x[0,L]
         do i=0,L-1
            At(i-L)=xAt(i)   
            Abt(i-L)=xAbt(i)   
         enddo  
         if(iloop.eq.nloop)then
            open(81,file='A_t.ipt')
            open(91,file='Ab_t.ipt')
            do i=-L,L
               write(81,*)dfloat(i)*dt,real(At(i)),dimag(At(i))
               write(91,*)dfloat(i)*dt,real(Abt(i)),dimag(Abt(i))
            enddo
            close(81)
            close(91)
         endif
         
!     Build the G^{--} & G^{+-} \= G_{\tbar} & G_{>}
!--------------------------
         open(13,file='Ggreater_t.ipt')
         open(14,file='Glesser_t.ipt')
         open(15,file='G_tordered_t.ipt')
         do i=-L,L
            g0tpm(i)=-pi2*xi*(At(i)-Abt(i))
            g0tmp(i)=pi2*xi*Abt(i)
            if(i.lt.0)then
               g0tmm(i)=g0tmp(i)
            else
               g0tmm(i)=g0tpm(i)
            endif
            g0t(i)=g0tmm(i)-g0tpm(i)
            if(i.eq.0)g0t(i)=-g0tpm(i)

            if(iloop.eq.nloop)then
               write(13,*)dfloat(i)*dt,dimag(g0tpm(i)),real(g0tpm(i))
               write(14,*)dfloat(i)*dt,dimag(g0tmp(i)),real(g0tmp(i))
               write(15,*)dfloat(i)*dt,dimag(g0tmm(i)),real(g0tmm(i))
            endif
         enddo
         do i=13,15
            close(i)
         enddo
         
!     Impurity Solver: Sigma=U^2G0^3
!-----------------------------
!     get Sigma
         do i=-L+1,L-1
            sigmatmm(i)=(U**2)*(g0tmm(i)**2)*g0tmm(-i)
            sigmatpm(i)=-(U**2)*(g0tpm(i)**2)*g0tmp(-i)
         enddo
         sigmatmm(L)=(U**2)*(g0mm(L)**2)*g0tmm(L-1)
         sigmatpm(L)=-(U**2)*(g0tpm(L)**2)*g0tmp(L-1)

         do i=-L,L
            sigmat(i)=sigmatmm(i)+sigmatpm(i)
         enddo
         sigmat(0)=sigmatmm(0)

         if(iloop.eq.nloop)then
            open(73,file='Sigma_tordered_t.ipt')
            open(83,file='Sigma_greater_t.ipt')
            do i=-L,L
               write(73,*)dfloat(i)*dt,dimag(sigmatmm(i)),
     *              real(sigmatmm(i))
               write(83,*)dfloat(i)*dt,dimag(sigmatpm(i)),
     *              real(sigmatpm(i))
            enddo
            close(73)
            close(83)
         endif

!     FFT t-->w 
!-------------------
!     starts manipulations
         do i=1,n
            xsigmat(i)=sigmat(i-L-1)
            gc(i)=g0t(i-L-1) 
         enddo
         
         CALL four1(xsigmat,n,1)
         CALL four1(gc,n,1)
         ex=-1.d0
         do i=1,n
            ex=-ex
            sigma(i)=ex*dt*xsigmat(i)
            fg0(i)=ex*dt*gc(i)  !\cal G_0Adv.
            fg(i)=one/fg0(i)-sigma(i) !G^-1=G0^-1-Sigma
            fg(i)=one/fg(i)
         enddo

!     End DMFT-ITERATION
!===========================================================

!     
!     Write out
!-----------------------
         if(iloop.ge.nloop)then
            open(10,file='Sigmaadv_t.ipt')
            open(11,file='Gadv_t.ipt')
            do i=-L,L
               write(10,*)dfloat(i)*dt,dimag(sigmat(i)),real(sigmat(i))
               write(11,*)dfloat(i)*dt,dimag(g0t(i)),real(g0t(i))
            enddo
            close(10)
            close(11)

            open(10,file='Gadv_realw.ipt')
            open(11,file='Sigmaadv_realw.ipt')
            open(12,file='G0adv_realw.ipt')
            do i=L+1,2*L
               write(10,*)ome(i),dimag(fg(i)),real(fg(i))
               write(11,*)ome(i),dimag(sigma(i)),real(sigma(i))
               write(12,*)ome(i),dimag(fg0(i)),real(fg0(i))
            enddo
            do i=1,L
               write(10,*)ome(i),dimag(fg(i)),real(fg(i))
               write(11,*)ome(i),dimag(sigma(i)),real(sigma(i))
               write(12,*)ome(i),dimag(fg0(i)),real(fg0(i))
            enddo
            close(10)
            close(11)
            close(12)

            open(10,file='DOSadv.ipt')
            do i=L+1,2*L
               gg=dimag(fg(i))
               if(abs(ome(i)).lt.5.d0)then
                  write(10,*)ome(i),gg
               endif
            enddo
            do i=1,L
               gg=dimag(fg(i))
               if(abs(ome(i)).lt.5.d0)then
                  write(10,*)ome(i),gg
               endif
            enddo
            close(10)
         endif

      enddo                     !here the dmft loops end up
c=======================================================================
      end
      include 'routines.f'
c     include 'fftpack5.f'
