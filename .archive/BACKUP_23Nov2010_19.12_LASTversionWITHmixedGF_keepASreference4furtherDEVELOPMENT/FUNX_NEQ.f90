!###############################################################
!     PROGRAM  : FUNCS_NEQ
!     TYPE     : Module
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!     LAST UPDATE: 07/2009
!###############################################################
module FUNX_NEQ
  !LIBRARY:  
  USE MATRIX
  !LOCAL:
  USE VARS_GLOBAL
  implicit none
  private
  public :: guess_G0  
  public :: get_Bath
  public :: update_G0_Dyson
  public :: get_Sigma
  public :: obtain_Gimp
  public :: get_Gloc_equilibrium

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : GUESS_G0
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial Weiss Fields G0^{<,>} as a non-interacting
  !guess. In absence of interaction the guess is time translation invariant. 
  !This is used to obtain the first Self-energy:
  !+-------------------------------------------------------------------+
  subroutine guess_G0()
    integer    :: i,j,ik
    real(8)    :: en,intE,A
    complex(8) :: peso
    real(8)    :: nless,ngtr
    complex(8),dimension(0:nstep,0:nstep) :: G0ret
    call dump("Get G0guess(t,t'):")    
    call system("mkdir GUESS")

    g0tgtr=zero; g0tless=zero
    G0gtr=zero ; G0less=zero
    if(Efield==0.d0)then        !If time-translation invariance is present
       if(irdeq)then            !Read from equilibrium solution
          do ik=1,2*L
             en   = eqwr(ik)
             nless= fermi0(en,beta)
             ngtr = fermi0(en,beta)-1.d0
             A    = -aimag(icG0w(ik))/pi
             do i=-nstep,nstep
                peso=exp(-xi*en*t(i))
                g0tless(i)=g0tless(i) + xi*nless*A*peso*rfmesh
                g0tgtr(i) =g0tgtr(i)  + xi*ngtr*A*peso*rfmesh
             enddo
          enddo

       else                     !Create your first guess in t-t'

          do ik=1,Lk
             en   = epsik(ik)
             nless= fermi0(en,beta)
             ngtr = fermi0(en,beta)-1.d0
             do i=-nstep,nstep
                peso=exp(-xi*en*t(i))
                g0tless(i)=g0tless(i) + xi*nless*peso*wt(ik)
                g0tgtr(i) =g0tgtr(i)  + xi*ngtr*peso*wt(ik)
             enddo
          enddo
       endif

       forall(i=0:nstep,j=0:nstep)
          G0less(i,j)=g0tless(i-j)
          G0gtr(i,j) =g0tgtr(i-j)
       end forall

    elseif(Efield/=0.d0)then    !If time-translation invariance is NOT present

       do ik=1,Lk
          call print_bar(ik,Lk)
          en   = epsik(ik)
          nless= fermi0(en,beta)
          ngtr = fermi0(en,beta)-1.d0
          do i=0,nstep
             do j=0,nstep
                intE=int_Ht(ik,i,j)
                peso=exp(-xi*intE)
                G0less(i,j)= G0less(i,j) + xi*nless*peso*wt(ik)
                G0gtr(i,j) = G0gtr(i,j)  + xi*ngtr*peso*wt(ik)
             enddo
          enddo
       enddo
       call close_bar
    endif

    forall(i=0:nstep,j=0:nstep)g0tret(i-j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
    call cfft_rt2rw(g0tret,g0fret,nstep) ; g0fret=g0fret*dt ; call swap_fftrt2rw(g0fret)
    forall(i=0:nstep,j=0:nstep)
       g0tless(i-j)=G0less(i,j)
       g0tgtr(i-j)  =G0gtr(i,j) 
    end forall
    call splot("GUESS/guessG0less_t.ipt",t,g0tless)
    call splot("GUESS/guessG0gtr_t.ipt",t,g0tgtr)
    call splot("GUESS/guessG0ret_t.ipt",t,g0tret)
    call splot("GUESS/guessG0ret_realw.ipt",wr,g0fret)
    if(plot3D)then
       call plot_3D("guessG0less3D","$\Delta t$","$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
       call plot_3D("guessG0gtr3D","$\Delta t","$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0gtr(0:nstep,0:nstep))
       call system("mv guess*3D GUESS/")
    endif

#ifdef _mix
    forall(i=0:nstep,j=0:nstep)G0ret(i,j)=heaviside(t(i-j))*(G0gtr(i,j)-G0less(i,j))
    forall(i=0:nstep)G0ret(i,i)=-xi
    forall(i=0:nstep,j=0:Ltau)G0lceil(i,j) = -G0ret(i,0)*(-icGtau(Ltau-j))
    forall(j=0:Ltau)G0rceil(j,:)=conjg(G0lceil(:,Ltau-j))
    call plot_3D("guessG0lceil3D","$\Delta t$","$\tau$","Z",t(0:nstep)/dt,tau,G0lceil(0:nstep,0:Ltau))
    call plot_3D("guessG0rceil3D","$\tau$","$\Delta t$","Z",tau,t(0:nstep)/dt,G0rceil(0:Ltau,0:nstep))
    call system("mv guess*3D GUESS/")
#endif
    call dump("")

    return

  contains
    function int_Ht(ik,it,jt)
      real(8) :: int_Ht
      integer :: i,j,ii,ik,it,jt,sgn
      type(vect2D) :: kt,Ak
      int_Ht=0.d0
      if(it==jt)return
      sgn=1 
      if(jt > it)sgn=-1
      i=ik2ix(ik); j=ik2iy(ik)
      do ii=jt,it,sgn
         Ak=Afield(t(ii),Ek)
         kt=kgrid(i,j) - Ak     !t(ii)*Ek
         int_Ht=int_Ht + sgn*epsk(kt)*dt
      enddo
      return
    end function int_Ht
  end subroutine guess_G0
  !********************************************************************
  !********************************************************************
  !********************************************************************









  !+-------------------------------------------------------------------+
  !PROGRAM  : BUILD_BATH
  !TYPE     : Subroutine
  !PURPOSE  : Build the Bath part of the system using exclusively time
  !dependent formalism. The bath is not interacting so it is 
  ! time translation invariant.
  !+-------------------------------------------------------------------+
  subroutine get_Bath()
    integer :: i,im,j
    real(8) :: en,w
    complex(8) :: peso,iw
    real(8) :: ngtr,nless,arg
    complex(8),dimension(-nstep:nstep)   :: Gbret
#ifdef _mix
    complex(8),dimension(2*L)            :: Gbiw
    real(8),dimension(0:L)               :: Gbt(0:L),Gbtau(0:Ltau)
    complex(8),dimension(0:nstep,0:Ltau) :: Gblceil
    complex(8),dimension(0:Ltau,0:nstep) :: Gbrceil
#endif

    call dump("Get Bath:")    
    S0less=zero ; S0gtr=zero
    do im=1,Lmu 
       en   = epsimu(im)
       nless= fermi0(en,beta)
       ngtr = fermi0(en,beta)-1.d0
       do i=-nstep,nstep
          peso=exp(-xi*t(i)*en)
          S0less(i)=S0less(i)+ xi*Vpd**2*nless*densmu(im)*peso*de
          S0gtr(i) =S0gtr(i) + xi*Vpd**2*ngtr*densmu(im)*peso*de
       enddo
    enddo
    call system("mkdir BATH")
    call splot("BATH/S0less_t.ipt",t,S0less)
    call splot("BATH/S0gtr_t.ipt",t,S0gtr)

#ifdef _mix
    Gbiw=zero
    do i=1,2*L
       w=pi/beta*dble(2*i-1);iw=cmplx(0.d0,w)
       do im=0,Lmu
          Gbiw(i)=Gbiw(i)+de*densmu(im)/(iw-epsimu(im))!/dble(Lmu)
       enddo
    enddo
    call cfft_iw2tau(Gbiw,Gbt,beta) ; call extract(Gbt,Gbtau)
    forall(i=-nstep:nstep)Gbret(i)=heaviside(t(i))*(S0gtr(i)-S0less(i))/Vpd**2;Gbret(0)=-xi
    forall(i=0:nstep,j=0:Ltau)Gblceil(i,j)=Gbret(i)*Gbtau(Ltau-j)
    forall(j=0:Ltau)Gbrceil(j,:)=conjg(Gblceil(:,Ltau-j))
    S0lceil=Vpd**2*Gblceil
    S0rceil=Vpd**2*Gbrceil
    if(Vpd/=0.d0)then
       call plot_3D("S0lceil3D","$\Delta t","$\tau","Z",t(0:nstep)/dt,tau,S0lceil(0:nstep,0:Ltau))
       call plot_3D("S0rceil3D","$\tau","$\Delta t","Z",tau,t(0:nstep)/dt,S0rceil(0:Ltau,0:nstep))
       call system("mv S0*3D BATH/")
    endif
#endif

    call dump("")
    return
  end subroutine Get_Bath
  !********************************************************************
  !********************************************************************
  !********************************************************************









  !+-------------------------------------------------------------------+
  !PROGRAM  : GET_SIGMA
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine get_Sigma(method_)
    character(len=*)                      :: method_
    integer                               :: i,j,itau
    complex(8),dimension(0:nstep,0:nstep) :: impGless,impGgtr
    call dump("Get Sigma(t,t')")
    forall(i=0:nstep,j=0:nstep)
       Sless(i,j) = (U**2)*(G0less(i,j)**2)*G0gtr(j,i) + S0less(i-j) !add the HF-like  bath contribution
       Sgtr (i,j) = (U**2)*(G0gtr(i,j)**2)*G0less(j,i) + S0gtr(i-j)  !add the HF-like  bath contribution
    end forall

#ifdef _mix
    Slceil=zero ; Srceil=zero
    forall(i=0:nstep,itau=0:Ltau)Slceil(i,itau) = U**2*(G0lceil(i,itau)**2)*G0rceil(itau,i)
    forall(itau=0:Ltau)Srceil(itau,:)=conjg(Slceil(:,Ltau-itau))
#endif

    !SPT impurity:
    if(trim(method_) == "spt")then
       call dump("Doing SPT")
       call obtain_Gimp(impGless,impGgtr)
       forall(i=0:nstep,j=0:nstep)
          Sless(i,j) = U**2*(impGless(i,j)**2)*impGgtr(j,i)
          Sgtr (i,j) = U**2*(impGgtr(i,j)**2)*impGless(j,i)                
       end forall
    endif
    call dump("")
    return
  end subroutine Get_Sigma
  !********************************************************************
  !********************************************************************
  !********************************************************************










  !+-------------------------------------------------------------------+
  !PROGRAM  : UPDATE_G0_DYSON 
  !TYPE     : Subroutine
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine update_G0_Dyson()
    integer :: M,i,j,k,itau,jtau,NN
    real(8) :: R,deg
    real(8) :: w,A,An
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet
    complex(8),dimension(0:nstep,0:nstep) :: G0ret,locGret,Sret
    complex(8),dimension(0:nstep,0:nstep) :: G0adv,locGadv,Sadv
    complex(8),dimension(0:nstep,0:nstep) :: dG0less,dG0gtr
    complex(8),dimension(0:nstep,0:nstep) :: G0kel,locGkel,Skel
    complex(8),dimension(0:nstep,0:nstep) :: locGtt,locGat,Stt,Sat
    complex(8),dimension(:,:),allocatable :: locGmat,Smat,G0mat,GammaMat,UnoMat
#ifdef _mix
    complex(8),dimension(0:Ltau,0:Ltau)   :: GammaM,locGtau,Stau
    complex(8),dimension(2*L)             :: GammafM
    real(8),dimension(0:L)                :: GammatM_
    real(8),dimension(-Ltau:Ltau)         :: GammatM
    complex(8) :: Op1lceil(0:nstep,0:Ltau)
    complex(8) :: Op2lceil(0:nstep,0:Ltau)
    complex(8) :: Op3less(0:nstep,0:Ltau)
    complex(8) :: Op1less(0:Ltau,0:nstep)
    complex(8) :: Op2less(0:Ltau,0:nstep)
#endif


    call dump("Update WF: Dyson")
    !Frequency domain procedure: works for E=0, time-translation invariant
    !========================================================
    if(Efield==0.d0 .AND. wfftw)then
       include "update_G0_equilibrium.f90"
    else
       include "update_G0_nonequilibrium.f90"
#ifdef _mix
       include "update_G0hookd_dyson.f90"
#endif
    endif
    call splot("dG0less_t.ipt",t,g0tless)
    call splot("dG0gtr_t.ipt",t,g0tgtr)
    call splot("dG0ret_t.ipt",t,g0tret)
    call splot("dG0ret_realw.ipt",wr,g0fret)
    if(plot3D)then
       call plot_3D("dG0less3D","$\Delta t$","$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
       call plot_3D("dG0gtr3D","$\Delta t","$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0gtr(0:nstep,0:nstep))
       call plot_3D("dG0ret3D","$\Delta t","$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0ret(0:nstep,0:nstep))
    endif
    call dump("")
    return
  end subroutine update_G0_Dyson
  !********************************************************************
  !********************************************************************
  !********************************************************************
























  !+-------------------------------------------------------------------+
  !PROGRAM  : OBTAIN_GIMP
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine obtain_Gimp(dGless,dGgtr)
    integer                               :: i,j
    real(8)                               :: A,w
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet,Gamma0Ret
    complex(8),dimension(0:nstep,0:nstep) :: dG0ret
    complex(8),dimension(0:nstep,0:nstep) :: dGless,dGgtr,dGret
    complex(8),dimension(0:nstep,0:nstep) :: dSret

    if(Efield==0.d0 .AND. wfftw)then
       include "obtain_Gimp_equilibrium.f90"
       forall(i=0:nstep,j=0:nstep)
          dGless(i,j)= gtless(i-j)
          dGgtr(i,j) = gtgtr(i-j)
       end forall
    else
       dSret=zero ; dG0ret=zero ; dGret=zero ; dGless=zero ; dGgtr=zero
       GammaRet=zero ; Gamma0Ret=zero
       !1 - get the Ret components of G_0 && \Sigma:
       forall(i=0:nstep,j=0:nstep)
          dG0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j) - G0less(i,j))
          dSret(i,j) =heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
       end forall
       !2 - get the  operator: \Gamma_0^R = \Id - \Sigma^R\circ G_0^R && invert it
       Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt
       Gamma0Ret(0:nstep,0:nstep) = Uno-matmul(dSret(0:nstep,0:nstep),dG0ret(0:nstep,0:nstep))*dt
       Gamma0Ret(0:nstep,0:nstep)=Gamma0Ret(0:nstep,0:nstep)*dt**2
       !call mat_inversion(Gamma0Ret(0:nstep,0:nstep))
       call mat_inversion_GJ(Gamma0Ret(0:nstep,0:nstep))
       !3 - get G_imp^R, G_imp^{>,<} using Dyson equations:
       dGret(0:nstep,0:nstep)    = matmul(dG0ret(0:nstep,0:nstep),Gamma0Ret(0:nstep,0:nstep))*dt 
       GammaRet(0:nstep,0:nstep) = Uno+matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt
       dGless(0:nstep,0:nstep)   = matmul(GammaRet(0:nstep,0:nstep),matmul(G0less(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2 +&
            matmul(dGret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
       dGgtr(0:nstep,0:nstep)    = matmul(GammaRet(0:nstep,0:nstep),matmul(G0gtr(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2  +&
            matmul(dGret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
    endif
    return
  end subroutine obtain_Gimp
  !********************************************************************
  !********************************************************************
  !********************************************************************







  !+-------------------------------------------------------------------+
  !PROGRAM  : GET_GLOC_EQUILIBRIUM
  !TYPE     : function
  !PURPOSE  : Evaluate the k-sum to get the local Green's functions (R,>,<)
  !+-------------------------------------------------------------------+
  subroutine get_Gloc_equilibrium()
    integer    :: i,j,ik
    complex(8) :: A,zetan
    real(8)    :: w,n
    complex(8) :: funcM(2*L),sigma(2*L)
    real(8)    :: funcT(0:L) 

    !Get Sret(w) = FFT(Sret(t-t'))
    forall(i=0:nstep,j=0:nstep) stret(i-j)=heaviside(t(i-j))*(Sgtr(i,j)-Sless(i,j))
    stret=exa*stret ; call cfft_rt2rw(stret,sfret,nstep) ; sfret=dt*sfret

    !Get locGret(w),locG<(w),locG>(w)
    gfret=zero
    do i=1,2*nstep
       w=wr(i)
       zetan=cmplx(w,eps)-sfret(i) !-eqsbfret(i)
       do ik=1,Lk
          gfret(i)=gfret(i)+wt(ik)/(zetan-epsik(ik))
       enddo
       A=-aimag(gfret(i))/pi
       gfless(i)= pi2*xi*fermi0(w,beta)*A
       gfgtr(i) = pi2*xi*(fermi0(w,beta)-1.d0)*A
    enddo
    call splot("locGret_realw.ipt",wr,gfret)

    !Get locG<(t),locG>(t)
    call cfft_rw2rt(gfless,gtless,nstep)  ; gtless=fmesh/pi2*gtless
    call cfft_rw2rt(gfgtr,gtgtr,nstep)    ; gtgtr=fmesh/pi2*gtgtr
    gtless=exa*gtless
    gtgtr =exa*gtgtr
    call splot("locGless_t.ipt",t,gtless)
    call splot("locGgtr_t.ipt",t,gtgtr)

    forall(i=0:nstep,j=0:nstep)
       locGless(i,j) = gtless(i-j)
       locGgtr(i,j)  = gtgtr(i-j)
       gtret(i-j)  = heaviside(t(i-j))*(locGgtr(i,j)-locGless(i,j))
    end forall
    call splot("locGret_t.ipt",t,gtret)

    call getGmats(wr,sfret,sigma,beta)
    call splot('locSM_iw.ipt',wm,sigma)
    do ik=1,Lk
       funcM=zero
       do i=1,2*L
          w=pi/beta*dble(2*i-1) ; zetan=cmplx(0.d0,w) - sigma(i)
          funcM(i)=one/(zetan - epsik(ik))
       enddo
       call cfft_iw2tau(funcM,funcT,beta)
       n=-funcT(L) !from G(-tau) = -G(beta-tau)==> G(tau=0-)=-G(beta)
       nk(:,ik)=n
    enddo

    return
  end subroutine get_Gloc_equilibrium
  !******************************************************************
  !******************************************************************
  !******************************************************************
































  !+-------------------------------------------------------------------+
  !PROGRAM  : BUILD_BATH_W
  !TYPE     : Subroutine
  !PURPOSE  : Build the Bath part of the system using frequency domain
  !formalism. 
  !COMMENT  : Construct the Bath-Sigma functions from freq. domain
  !(that is 'cause you know the analytic expression, same
  !procedure as in IPTkeldysh @equilibrium)
  !This turns out to be the same provided one takes 
  !the limit \eta\goto 0.d0
  !+-------------------------------------------------------------------+
  subroutine get_Bath_w
    integer :: i,im
    real(8) :: A,An,w
    complex(8) :: iw,fg
    complex(8),dimension(2*nstep)     :: gb0fless,gb0fgtr
    complex(8),dimension(-nstep:nstep):: gb0tgtr,gb0tless

    call dump("Get Bath:")    
    do i=1,2*nstep
       w = wr(i) ;iw= cmplx(w,eps)
       fg=zero
       do im=0,Lmu
          fg=fg + de*densmu(im)/(iw-epsimu(im))!/dble(Lmu)
       enddo
       A = -aimag(fg)/pi
       An= A*fermi0(w,beta)
       gb0fless(i)= pi2*xi*An
       gb0fgtr(i) = pi2*xi*(An-A)
    enddo
    call cfft_rw2rt(gb0fless,gb0tless,nstep)  ; gb0tless=fmesh/pi2*gb0tless
    call cfft_rw2rt(gb0fgtr, gb0tgtr,nstep)   ; gb0tgtr =fmesh/pi2*gb0tgtr
    gb0tgtr =exa*gb0tgtr
    gb0tless=exa*gb0tless
    !Get the Self-energy contribution of the bath:
    !=================================================
    S0gtr=zero
    S0less=zero
    S0less=Vpd**2*gb0tless
    S0gtr =Vpd**2*gb0tgtr
    call system("mkdir BATH")   
    call splot("BATH/bathG0less_t.ipt",t(-nstep:nstep),gb0tless)
    call splot("BATH/bathG0gtr_t.ipt",t(-nstep:nstep),gb0tgtr)
    call splot("BATH/S0less_t.ipt",t(-nstep:nstep),S0less)
    call splot("BATH/S0gtr_t.ipt",t(-nstep:nstep),S0gtr)
    call dump("")
    return
  end subroutine Get_Bath_w
  !********************************************************************
  !********************************************************************
  !********************************************************************




end module FUNX_NEQ
