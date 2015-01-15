!########################################################################
!     PROGRAM  : EQ_SOLUTION_KELDYSH
!     TYPE     : module
!     PURPOSE  : Solve the Hubbard model at EQUILIBRIUM using DMFT and 
!                Keldysh technique
!     AUTHORS  : Adriano Amaricci
!     LAST UPDATE: 10/2009: based on the F77 version (working)
!########################################################################
module SOLUTION_EQ
  use FFTW_MOD

  implicit none
  private 
  public eq_hmkeldysh,fermi,dens_bethe,dens_hyperc
contains 
  subroutine eq_hmkeldysh
    integer :: nloop,iloop,i
    integer, parameter :: L=2*4096, n=2*L

    real(8) :: fmesh,D,U,dt,pi,beta,ex,w,sig,gg,xp,xm,sq
    real(8) :: A,Ab,pi2,si,sr,s0
    real(8),dimension(2*L) :: ome 

    complex(8) :: xi,one,zero,D2,iw, sqroot
    complex(8),dimension(2*L) :: fg0,fg,gc,g0mm,sigma
    complex(8),dimension(2*L) :: xsigmat,Adummy,Abdummy
    complex(8),dimension(-L:L):: g0tmm,g0tpm,g0tmp,g0t
    complex(8),dimension(-L:L):: xAt,xAbt,At,Abt
    complex(8),dimension(-L:L):: sigmatmm,sigmatpm,sigmat


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

    print*, '       U = ',U
    print*, '    beta = ',beta
    print*, '    Mesh = ',fmesh
    !===============================================================

    !initialize some functions
    fg0=zero
    fg=zero
    sigma=zero
    do i=1,n
       ome(i)=real(i-L-1,8)
    enddo
    ome=ome*fmesh

    !Starts DMFT-loop
    !===============================================================
    do iloop=1,nloop
       print*,'iloop',iloop,'/',nloop

       !     First loop: guess
       if(iloop.eq.1)then
          do i=1,n
             w=ome(i)
             sqroot=cdsqrt(w**2-D2)
             fg0(i)=2.d0/(w-sqroot) !G_0=G for Sigma=0
             !fg0(i)=one/(w+xi*dfloat(imet)/2.d0)!marcelo's guess
             !if(abs(w).lt.1.d-9)fg0(i)=0.d0*one
          enddo
       else
          !SELF CONSISTENCY CONDITION
          !-------------------------------------------------------------
          fg0=ome-xi*0.01d0-fg/4.d0
          fg0=one/fg0
       endif
       !Get w<0 function
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
          Adummy(i)=cmplx(A,0.d0)
          Abdummy(i)=cmplx(Ab,0.d0)
          if(iloop.eq.nloop)then
             write(11,*)w,A,Ab
          endif
       enddo
       close(11)

       !FFT w-->t
       !Manipulations required by the structure of FFTW in real space/real time
       call cfft_rw2rt(Adummy,n)
       call manip_fftrw2rt(Adummy,At,L)
       call cfft_rw2rt(Abdummy,n)
       call manip_fftrw2rt(Abdummy,Abt,L)

       At=fmesh/pi2*At
       Abt=fmesh/pi2*Abt
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

       !Build the G^{--} & G^{+-} \= G_{\tbar} & G_{>}
       open(13,file='Gbigger_t.ipt')
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

       !Impurity Solver: Sigma=U^2G0^3
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
             write(73,*)dfloat(i)*dt,dimag(sigmatmm(i)),real(sigmatmm(i))
             write(83,*)dfloat(i)*dt,dimag(sigmatpm(i)),real(sigmatpm(i))
          enddo
          close(73)
          close(83)
       endif

       !FFT t-->w 
       do i=1,n
          xsigmat(i)=sigmat(i-L-1)
          gc(i)=g0t(i-L-1) 
       enddo

       call cfft_rt2rw(xsigmat,n)
       call cfft_rt2rw(gc,n)
       do i=1,n
          sigma(i)=dt*xsigmat(i)
          fg0(i)=dt*gc(i)  !\cal G_0Adv.
          fg(i)=one/fg0(i)-sigma(i) !G^-1=G0^-1-Sigma
          fg(i)=one/fg(i)
       enddo
       !End DMFT-ITERATION
       !===========================================================
       !Write out
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
          do i=1,2*L
             write(10,*)ome(i),dimag(fg(i)),real(fg(i))
             write(11,*)ome(i),dimag(sigma(i)),real(sigma(i))
             write(12,*)ome(i),dimag(fg0(i)),real(fg0(i))
          enddo
          close(10)
          close(11)
          close(12)

          open(10,file='DOSadv.ipt')
          do i=1,2*L
             gg=dimag(fg(i))
             if(abs(ome(i)).lt.5.d0)then
                write(10,*)ome(i),gg
             endif
          enddo
          close(10)
       endif
    enddo                     !here the dmft loops end up
  end subroutine eq_hmkeldysh

  !##################################################################
  !     FUNCTIONS
  !##################################################################
  !+-------------------------------------------------------------------+
  !PROGRAM  : FERMI
  !TYPE     : function
  !PURPOSE  : calculate the Fermi-Dirac distribution
  !VERSION  : 31-01-2006
  !+-------------------------------------------------------------------+
  FUNCTION fermi(x, mu,beta)
    real(8) :: fermi, x, mu, xeff, beta
    xeff=x-mu
    if(abs(xeff).le.0.5d0)then
       fermi = 1.d0/(1.d0+exp(beta*(xeff)))      
    else
       if(xeff.gt.0.d0)fermi=0.d0
       if(xeff.lt.0.d0)fermi=1.d0
    endif
    return
  END FUNCTION fermi

  !+-------------------------------------------------------------------+
  !PROGRAM  : DENS_BETHE
  !TYPE     : function
  !PURPOSE  : calculate the non-interacting dos for BETHE lattice 
  !VERSION  : 31-01-2006
  !+-------------------------------------------------------------------+
  FUNCTION dens_bethe(x)
    REAL(8) :: dens_bethe,x,t1,D
    complex(8):: root
    D=(1.d0,0.d0)      
    root=dcmplx((1.d0-1.d0*((x/D))**2),0.d0)
    root=sqrt(root)
    dens_bethe=(2.d0/(3.141592653589793238d0*D))*root
    return
  END FUNCTION dens_bethe

  !+-------------------------------------------------------------------+
  !PROGRAM  : DENS_HYPERC
  !TYPE     : function
  !PURPOSE  : calculate the non-interacting dos for HYPERCUBIC lattice 
  !VERSION  : 31-01-2006
  !+-------------------------------------------------------------------+
  FUNCTION dens_hyperc(x)
    REAL(8):: dens_hyperc,x,t1,D
    complex(8):: root
    t1=1.d0/sqrt(2.d0)
    dens_hyperc = (1/(t1*sqrt(2*3.141592653589793238d0)))*exp(-(x**2)/(2*t1**2))
    return
  END FUNCTION dens_hyperc

end module SOLUTION_EQ

