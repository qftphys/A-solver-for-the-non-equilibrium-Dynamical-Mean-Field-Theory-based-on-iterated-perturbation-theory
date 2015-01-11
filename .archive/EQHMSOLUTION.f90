MODULE EQ_SOLUTION
  !####################################################################
  !PROGRAM  : EQ_SOLUTION_KELDYSH
  !TYPE     : module
  !PURPOSE  : Solve the Hubbard model at EQUILIBRIUM using DMFT and 
  !Keldysh technique
  !AUTHORS  : Adriano Amaricci
  !LAST UPDATE: 10/2009: based on the F77 version (working)
  !COMMENTS : On output the 
  ! - real frequency Weiss field \cal{G}_0(\w) [BATH GF]
  ! - Sigma_Adv on real axis
  !####################################################################
  USE VARS_GLOBAL !definisce le variabili globali sharate dove serve
  USE FUNX_GLOBAL !definice funzioni utili
  USE FFTW        !MKL Fast Fourier Transforms 
  implicit none
  private 
  public eq_hmkeldysh
contains 
  subroutine eq_hmkeldysh(neqloop,sigma,sigmatau,fg0,fg00tau)
    integer :: ieqloop,neqloop,i,j
    real(8) :: w,time,shift
    real(8) :: A,Ab,e
    complex(8) :: zetan
    complex(8),dimension(2*L) :: fg0,fg,sigma,fg0iw,fgiw,sigmaiw
    real(8),dimension(0:L)    :: fg0tau
    real(8),dimension(0:Ltau) :: fg00tau,fgtau,sigmatau
    complex(8),dimension(2*L) :: Adummy,Abdummy
    complex(8),dimension(-L:L):: g0tmm,g0tpm,g0tmp,g0t,g0tret
    complex(8),dimension(-L:L):: At,Abt
    complex(8),dimension(-L:L):: sigmatmm,sigmatpm,sigmat
    character(len=160) :: cmd
    cmd = "mkdir EQsolution";call system(trim(cmd))   

    print*, '       U = ',U
    print*, '    beta = ',beta
    print*, '    Mesh = ',fmesh
    !===============================================================

    !initialize some functions
    fg0=zero
    fg=zero
    sigma=zero

    !Starts DMFT-loop
    !===============================================================
    do ieqloop=1,neqloop
       fg=zero
       !Get G_loc=\int d\e \rho_0(\e) G(\e,w)
       do i=1,n
          w=wmin+dble(i-1)*fmesh!ome(i)
          zetan=w-xi*0.1d0-sigma(i)
          do j=1,Lepsi
             e=emin+dble(j-1)*de
             fg(i)=fg(i)+dens_hyperc_1D(e)/(zetan-e)*de
          enddo
       enddo

       select case (iloop)
       case (1) 
          fg0=fg           
       case default
          fg0=one/fg + sigma
          fg0=one/fg0
       end select

       !Get w<0 function
       do i=1,L-1
          fg0(2*L-i)=fg0(i+1)
       enddo
       !Fix the endpoint
       fg0(n)=fg0(1)

       !Build the spectral densities
       do i=1,n
          w=wmin+dble(i-1)*fmesh!ome(i)
          A=dimag(fg0(i))/pi
          Ab=A*fermi(w)
          Adummy(i)=cmplx(A,0.d0)
          Abdummy(i)=cmplx(Ab,0.d0)
       enddo

       !FFT w-->t
       call cfft_rw2rt(Adummy,n)
       call cfft_rw2rt(Abdummy,n)
       !Manipulations required by the structure of FFTW in real space/real time
       call manip_fftrw2rt(Adummy,At,L)
       call manip_fftrw2rt(Abdummy,Abt,L)
       At=fmesh/pi2*At
       Abt=fmesh/pi2*Abt

       !Build the G^{--} & G^{+-} \= G_{\tbar} & G_{>}
       do i=-L,L
          time=dble(i)*dt
          g0tpm(i)=-pi2*xi*(At(i)-Abt(i)) !G^>
          g0tmp(i)=pi2*xi*Abt(i)          !G^<
          g0tmm(i)=heaviside(-time)*g0tmp(i) + heaviside(time)*g0tpm(i)
          g0tret(i)=heaviside(time)*(g0tpm(i) - g0tmp(i))
       enddo
       g0tmm(0)=g0tpm(0)
       g0t=g0tmm-g0tpm !; g0t=g0t/2.d0
       g0t(0)=g0tmp(0)!-g0tpm(0)

       !Impurity Solver: Sigma=U^2G0^3
       do i=-L,L
          sigmatmm(i)=(U**2)*(g0tmm(i)**2)*g0tmm(-i)
          sigmatpm(i)=-(U**2)*(g0tpm(i)**2)*g0tmp(-i)
       enddo
       sigmatmm(L)=(U**2)*(g0tmm(L)**2)*g0tmm(L-1)
       sigmatpm(L)=-(U**2)*(g0tpm(L)**2)*g0tmp(L-1)
       sigmat=sigmatmm+sigmatpm
       sigmat(0)=sigmatmm(0)

       !FFT t-->w 
       do i=1,n
          Adummy(i)=sigmat(i-L-1)
          Abdummy(i)=g0t(i-L-1) 
       enddo

       call cfft_rt2rw(Adummy,n)
       call cfft_rt2rw(Abdummy,n)

       !Get G_imp & Sigma_imp
       sigma=dt*Adummy
       if(aimag(sigma(L)) < 0.d0)then
          shift=-aimag(sigma(L))
          sigma=sigma+xi*shift
       endif
       fg0=dt*Abdummy  !\cal G_0Adv.
       fg=one/fg0 - sigma
       fg=one/fg

       !End DMFT-ITERATION
       if(ieqloop == neqloop)then
          open(10,file='EQsolution/Sigmaadv_t.ipt')
          open(11,file='EQsolution/G0adv_t.ipt')
          open(12,file='EQsolution/G0less_t.ipt')
          open(13,file='EQsolution/G0gtr_t.ipt')
          open(14,file='EQsolution/G0ret_t.ipt')
          open(15,file='EQsolution/Sigmagtr_t.ipt')
          do i=-L,L
             A=1.d0;if(i >0)A=-1.d0             
             write(10,*)dfloat(i)*dt,abs(aimag(sigmat(i))),A*(real(sigmat(i)))
             write(11,*)dfloat(i)*dt,abs(aimag(g0t(i))),A*(real(g0t(i)))
             write(12,*)dfloat(i)*dt,abs(aimag(g0tmp(i))),A*(real(g0tmp(i)))
             write(13,*)dfloat(i)*dt,-abs(aimag(g0tpm(i))),A*(real(g0tpm(i)))
             write(14,*)dfloat(i)*dt,abs(aimag(g0tret(i))),A*(real(g0tret(i)))
             write(15,*)dfloat(i)*dt,abs(aimag(sigmatpm(i))),&
                  A*(real(sigmatpm(i)))
          enddo
          do i=10,15
             close(i)
          enddo

          call getGmats(fg0,fg0iw,beta)
          call getGmats(fg,fgiw,beta)
          call getGmats(sigma,sigmaiw,beta)
          open(10,file='EQsolution/G_iw.ipt')
          open(11,file='EQsolution/G0_iw.ipt')
          open(12,file='EQsolution/Sigma_iw.ipt')
          do i=1,n
             w=pi/beta*(2.d0*dble(i)-1.d0)
             write(10,*)w,aimag(fgiw(i)),real(fgiw(i))
             write(11,*)w,aimag(fg0iw(i)),real(fg0iw(i))
             write(12,*)w,aimag(sigmaiw(i)),real(sigmaiw(i))
          enddo
          close(10)
          close(11)

          call cfft_iw2it(fg0iw,fg0tau,beta)
          call extract(fg0tau,fg00tau);fg0tau=0.d0
          call cfft_iw2it(fgiw,fg0tau,beta)
          call extract(fg0tau,fgtau);fg0tau=0.d0
          call cfft_iw2it(sigmaiw,fg0tau,beta)
          call extract(fg0tau,sigmatau)

          open(10,file='EQsolution/G_tau.ipt')
          open(11,file='EQsolution/G0_tau.ipt')
          open(12,file='EQsolution/Sigma_tau.ipt')
          do i=0,Ltau
             write(10,*)dble(i)*beta/dble(Ltau),fgtau(i)
             write(11,*)dble(i)*beta/dble(Ltau),fg00tau(i)
             write(12,*)dble(i)*beta/dble(Ltau),sigmatau(i)
          enddo
          close(10);close(11);close(12)

          open(10,file='EQsolution/Gadv_realw.ipt')
          open(11,file='EQsolution/Sigmaadv_realw.ipt')
          open(12,file='EQsolution/G0adv_realw.ipt')
          open(13,file='EQsolution/DOSadv.ipt')         
          do i=1,n
             w=wmin+dble(i-1)*fmesh
             write(10,*)w,aimag(fg(i)),real(fg(i))
             write(11,*)w,aimag(sigma(i)),real(sigma(i))
             write(12,*)w,aimag(fg0(i)),real(fg0(i))
             write(13,*)w,aimag(fg(i))/pi
          enddo
          close(10)
          close(11)
          close(12)
          close(13)
       endif
    enddo                     !here the dmft loops end up

    write(*,*)"Got solution at EQUILIBRIUM "
    write(*,*)"----------------------------"
    return
  end subroutine eq_hmkeldysh
end module EQ_SOLUTION


