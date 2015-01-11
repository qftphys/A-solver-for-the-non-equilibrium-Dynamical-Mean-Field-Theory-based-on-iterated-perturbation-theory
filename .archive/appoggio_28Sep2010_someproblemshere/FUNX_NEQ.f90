!###############################################################
!     PROGRAM  : FUNCS_NEQ
!     TYPE     : Module
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!     LAST UPDATE: 07/2009
!###############################################################
module FUNX_NEQ
  !LIBRARY:  
  USE TOOLS
  USE GRIDS
  USE LATTICE
  USE PROGRESS
  USE MATRIX
  USE DLPLOT
  USE SPLINE
  USE FFTW
  !LOCAL:
  USE VARS_GLOBAL
  implicit none
  private
  public :: guess_G0  
  public :: get_Bath
  public :: update_G0_Dyson
  public :: get_Sigma
  public :: save_solution     
  public :: obtain_Gimp
  save
contains

  !+------------------------------------------------------------------+
  !PROGRAM  : Xc
  !TYPE     : function
  !PURPOSE  : build  the function X^c(e,t,t') equal to the term in front
  !of the exponential in the expressions for G_0(t,t') 
  !COMMENT  : The function is defined WITHOUT the imaginary unit xi
  !+------------------------------------------------------------------+
  function Xc(cc,x)
    real(8)          :: Xc,x
    character(len=1) :: cc
    Xc=0.d0    
    if(cc == '<')then
       Xc=fermi(x,beta)
    elseif(cc == '>')then
       Xc=fermi(x,beta)-1.d0
    endif
  end function Xc
  !********************************************************************
  !********************************************************************
  !********************************************************************













  !+-------------------------------------------------------------------+
  !PROGRAM  : GETG0FREE
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial Weiss Fields G0^{<,} as a non-interacting
  !guess. In absence of interaction the guess is time translation invariant. 
  !This is used to obtain the first Self-energy:
  !+-------------------------------------------------------------------+
  subroutine guess_G0()
    integer    :: i,j,ik
    real(8)    :: en
    complex(8) :: peso
    real(8)    :: nless,ngtr
    complex(8),dimension(0:nstep,0:nstep) :: G0ret
    call dump("Get G0guess(t,t'):")    
    call system("mkdir GUESS")
    G0gtr= zero ; G0less= zero
    if(eqflag)then !if eqflag you already have the equilibrium G0t_less/gtr(t-t')
       forall(i=0:nstep,j=0:nstep)
          G0less(i,j)=icG0tless(i-j)
          G0gtr(i,j) =icG0tgtr(i-j)
       end forall
       call splot("GUESS/guessG0less_t.ipt",t,icG0tless(-nstep:nstep))
       call splot("GUESS/guessG0gtr_t.ipt",t,icG0tgtr(-nstep:nstep))
    else
       g0tgtr=zero ; g0tless=zero
       do ik=1,Lk
          en   = epsik(ik)
          ngtr = Xc('>',en)
          nless= Xc('<',en)
          do i=-nstep,nstep
             peso=exp(-xi*en*t(i))
             g0tless(i)=g0tless(i) + xi*nless*peso*wt(ik)
             g0tgtr(i) =g0tgtr(i)  + xi*ngtr*peso*wt(ik)
          enddo
       enddo
       forall(i=0:nstep,j=0:nstep)
          G0less(i,j)=g0tless(i-j)
          G0gtr(i,j) =g0tgtr(i-j)
       end forall
       call splot("GUESS/guessG0less_t.ipt",t,g0tless)
       call splot("GUESS/guessG0gtr_t.ipt",t,g0tgtr)
    endif


#ifdef _mix
    ! G0lceil= zero
    ! do ik=1,Lk
    !    en   = epsik(ik)
    !    nless= Xc('<',en)
	!	 do i=0,nstep
    !      do j=0,Ltau
    !         peso=exp(-xi*en*t(i))
    !         peso=peso*(exp(en*tau(j))+exp(en*(beta-tau(j))))*nless
    !         G0lceil(i,j)= G0lceil(i,j) + xi*nless*peso*wt(ik)
    !      enddo
    ! 	 enddo
	! enddo


    do i=0,nstep
       do j=0,nstep
          G0ret(i,j)=heaviside(t(i-j))*(G0gtr(i,j)-G0less(i,j))
       enddo
       G0ret(i,i)=-xi
    enddo

    do i=0,nstep
       do j=0,Ltau
          G0lceil(i,j) = -xi*G0ret(i,0)*icGtau(Ltau-j)
       enddo
    enddo
    forall(j=0:Ltau)G0rceil(j,:)=conjg(G0lceil(:,Ltau-j))
    call plot_3D("guessG0lceil3D","$\Delta t$","$\tau$","Z",t(0:nstep)/dt,tau,G0lceil(0:nstep,0:Ltau))
    call plot_3D("guessG0rceil3D","$\tau$","$\Delta t$","Z",tau,t(0:nstep)/dt,G0rceil(0:Ltau,0:nstep))
    call system("rm GUESS/guess*3D 2>/dev/null; mv guess*3D GUESS/")
#endif

    call dump("")
    return
  end subroutine guess_G0
  !********************************************************************
  !********************************************************************
  !********************************************************************







  !+-------------------------------------------------------------------+
  !PROGRAM  : BUILD_BATH_T
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
    complex(8),dimension(-nstep:nstep)   :: Gbless,Gbgtr,Gbret
#ifdef _mix
    complex(8),dimension(2*L)            :: Gbiw
    real(8),dimension(0:L)               :: Gbt(0:L),Gbtau(0:Ltau)
    complex(8),dimension(0:nstep,0:Ltau) :: Gblceil
    complex(8),dimension(0:Ltau,0:nstep) :: Gbrceil
#endif

    call dump("Get Bath:")    
    Gbless=zero ; Gbgtr=zero
    do im=1,Lmu 
       en   = epsimu(im)
       nless= Xc("<",en)
       ngtr = Xc(">",en)
       do i=-nstep,nstep
          arg=t(i)*en
          peso=exp(-xi*arg)
          Gbless(i)=Gbless(i)+ xi*nless*peso/dble(Lmu)
          Gbgtr(i) =Gbgtr(i) + xi*ngtr*peso/dble(Lmu)
       enddo
    enddo
    S0less=Vpd**2*Gbless
    S0gtr =Vpd**2*Gbgtr
    call system("mkdir BATH")
    call splot("BATH/bathG0less_t.ipt",t,Gbless)
    call splot("BATH/bathG0gtr_t.ipt",t,Gbgtr)
    call splot("BATH/S0less_t.ipt",t,S0less)
    call splot("BATH/S0gtr_t.ipt",t,S0gtr)


#ifdef _mix
    ! Gblceil=zero
    ! do im=1,Lmu
    !    en   = epsimu(im)
    !    nless= Xc('<',en)
    !    do i=0,nstep
    !       do j=0,Ltau
    !          peso=exp(-xi*en*t(i))
    !          peso=peso*(exp(-abs(en)*tau(j))+exp(-abs(en)*(beta-tau(j))))*nless
    !          Gblceil(i,j) = Gblceil(i,j) + xi*nless*peso/dble(Lmu)
    !       enddo
    !    enddo
    ! enddo

    Gbiw=zero
    do i=1,2*L
       w=pi/beta*dble(2*i-1);iw=cmplx(0.d0,w)
       do im=1,Lmu
          Gbiw(i)=Gbiw(i)+one/(iw-epsimu(im))/dble(Lmu)
       enddo
    enddo
    call cfft_iw2tau(Gbiw,Gbt,beta) ; call extract(Gbt,Gbtau)
    do i=-nstep,nstep
       Gbret(i)=heaviside(t(i))*(Gbgtr(i)-Gbless(i))
    enddo
    Gbret(0)=-xi
    do i=0,nstep
       do j=0,Ltau
          Gblceil(i,j)=-xi*Gbret(i)*Gbtau(Ltau-j)
       enddo
    enddo
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
  !PROGRAM  : UPDATE_G0_DYSON 
  !TYPE     : Subroutine
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine update_G0_Dyson()
    integer :: i,j
    complex(8),dimension(0:nstep,0:nstep) :: Uno,GammaRet
    complex(8),dimension(0:nstep,0:nstep) :: dG0ret,dGret,dSret    
    real(8)                               :: w,A,An
#ifdef _mix
    complex(8),dimension(0:Ltau,0:Ltau)   :: GammaM,locGtau,Stau
    complex(8),dimension(2*L)             :: GammafM
    real(8),dimension(0:L)                :: GammatM_
    real(8),dimension(-Ltau:Ltau)         :: GammatM
    complex(8) :: Op0(0:nstep,0:Ltau)
    complex(8) :: Op1(0:Ltau,0:nstep)
    complex(8) :: Op2(0:Ltau,0:nstep)
#endif

    call dump("Update WF: Dyson")
    !Frequency domain procedure: works for E=0, time-translation invariant
    !========================================================
    if(Efield==0.d0 .AND. wfftw)then
       include "update_G0_dyson_equilibrium.f90"
       forall(i=0:nstep,j=0:nstep)
          G0less(i,j)= g0tless(i-j)
          G0gtr(i,j) = g0tgtr(i-j)
       end forall

    else
       !Real time domain procedure: always works, in principle:
       !========================================================
       dGret=zero 
       dSret=zero 
       dG0ret=zero
       GammaRet=zero
       do i=0,nstep
          do j=0,nstep
             dGret(i,j)=heaviside(t(i)-t(j))*(locGgtr(i,j) - locGless(i,j))
             dSret(i,j)=heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
          enddo
       enddo

       Uno=zero   ; forall(i=0:nstep)Uno(i,i)=One/dt
       GammaRet(0:nstep,0:nstep) = Uno+matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt
       GammaRet(0:nstep,0:nstep) = GammaRet(0:nstep,0:nstep)*dt**2
       call mat_inversion_GJ(GammaRet(0:nstep,0:nstep))
       dG0ret(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),dGret(0:nstep,0:nstep))*dt 

       !G0less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sless * gA
       G0less(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),matmul(locGless(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt -&
            matmul(dG0ret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),conjg(transpose(dG0ret(0:nstep,0:nstep))))*dt)*dt

       !G0gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sgtr * gA
       G0gtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),matmul(locGgtr(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt  -&
            matmul(dG0ret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),conjg(transpose(dG0ret(0:nstep,0:nstep))))*dt)*dt

#ifdef _mix
       GammafM = One + icSiw*icGiw ; GammafM=one/GammafM
       call cfft_iw2tau(GammafM,GammatM_,beta) ; call extract(GammatM_,GammatM(0:Ltau))
       forall(i=0:Ltau) GammatM(-i) = -GammatM(Ltau-i)
       forall(i=0:Ltau,j=0:Ltau) 
          GammaM(i,j)  = GammatM(i-j)
          locGtau(i,j) = icGtau(i-j)
          Stau(i,j)    = icStau(i-j)
       end forall

       !G0lceil  = GammaR^-1 * Glceil * GammaM^-1   -   gR * Slceil * gM
       G0lceil(0:nstep,0:Ltau) = matmul(GammaRet(0:nstep,0:nstep), matmul(locGlceil(0:nstep,0:Ltau),GammaM(0:Ltau,0:Ltau)) )*dt*dtau -&
            matmul(dG0ret(0:nstep,0:nstep), matmul(Slceil(0:nstep,0:Ltau),locGtau(0:Ltau,0:Ltau)) )*dt*dtau
       forall(i=0:Ltau)G0rceil(i,:)=conjg(G0lceil(:,Ltau-i))

       !G0less   = G0less -   gR * Slceil * grceil   -   GammaR^-1 * Glceil * (Srceil * gA + Sm * grceil)
       G0less(0:nstep,0:nstep) = G0less(0:nstep,0:nstep) + matmul(dG0ret(0:nstep,0:nstep),matmul(Slceil(0:nstep,0:Ltau),G0rceil(0:Ltau,0:nstep)))*dt*dtau-&
            matmul(matmul(GammaRet(0:nstep,0:nstep),locGlceil(0:nstep,0:Ltau))*dt,&
            (matmul(Srceil(0:Ltau,0:nstep),conjg(transpose(dG0ret(0:nstep,0:nstep))))*dt + matmul(Stau(0:Ltau,0:Ltau),G0rceil(0:Ltau,0:nstep))*dtau) )*dtau
       !G0gtr    = G0gtr  -   gR * Slceil * grceil   -  GammaR^-1 * Glceil * (Srceil * gA + Sm * grceil)
       G0gtr(0:nstep,0:nstep) =   G0gtr(0:nstep,0:nstep) + matmul(dG0ret(0:nstep,0:nstep),matmul(Slceil(0:nstep,0:Ltau),G0rceil(0:Ltau,0:nstep)))*dt*dtau-&
            matmul(matmul(GammaRet(0:nstep,0:nstep),locGlceil(0:nstep,0:Ltau))*dt,&
            (matmul(Srceil(0:Ltau,0:nstep),conjg(transpose(dG0ret(0:nstep,0:nstep))))*dt + matmul(Stau(0:Ltau,0:Ltau),G0rceil(0:Ltau,0:nstep))*dtau) )*dtau

       call plot_3D("dysonGtau3D","$\Delta t","$\tau","Z",tau,tau,locGtau(0:Ltau,0:Ltau))
       call plot_3D("dysonStau3D","$\Delta t","$\tau","Z",tau,tau,Stau(0:Ltau,0:Ltau))
       call plot_3D("dysonG0lceil3D","$\Delta t","$\tau","Z",t(0:nstep)/dt,tau,G0lceil(0:nstep,0:Ltau))
       call plot_3D("dysonG0less3D","$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
#endif

    endif
    call dump("")
    return
  end subroutine update_G0_Dyson
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
    !IPT impurity:
    Sless =zero
    Sgtr  =zero
    forall(i=0:nstep,j=0:nstep)
       Sless(i,j) = U**2*(G0less(i,j)**2)*G0gtr(j,i)
       Sgtr (i,j) = U**2*(G0gtr(i,j)**2)*G0less(j,i)
    end forall
#ifdef _mix
    Slceil=zero
    Srceil=zero
    forall(i=0:nstep,itau=0:Ltau)
       Slceil(i,itau) = U**2*(G0lceil(i,itau)**2)*G0rceil(itau,i)
    end forall
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
       dSret=zero 
       dG0ret=zero
       dGret=zero 
       dGless=zero 
       dGgtr=zero
       GammaRet=zero
       Gamma0Ret=zero
       !1 - get the Ret components of G_0 && \Sigma:
       do i=0,nstep
          do j=0,nstep
             dG0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j) - G0less(i,j))
             dSret(i,j) =heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
          enddo
       enddo
       !2 - get the  operator: \Gamma_0^R = \Id - \Sigma^R\circ G_0^R && invert it
       Uno=zero   
       forall(i=0:nstep)Uno(i,i)=One/dt
       Gamma0Ret(0:nstep,0:nstep) = Uno-matmul(dSret(0:nstep,0:nstep),dG0ret(0:nstep,0:nstep))*dt
       Gamma0Ret(0:nstep,0:nstep)=Gamma0Ret(0:nstep,0:nstep)*dt**2
       !call mat_inversion(Gamma0Ret(0:nstep,0:nstep))
       call mat_inversion_GJ(Gamma0Ret(0:nstep,0:nstep))
       !3 - get G_imp^R, G_imp^{>,<} using Dyson equations:
       dGret(0:nstep,0:nstep)    = matmul(dG0ret(0:nstep,0:nstep),Gamma0Ret(0:nstep,0:nstep))*dt 
       GammaRet(0:nstep,0:nstep) = Uno+matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt
       dGless(0:nstep,0:nstep)   = matmul(GammaRet(0:nstep,0:nstep),matmul(G0less(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2 +&
            matmul(dGret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
       dGgtr(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),matmul(G0gtr(0:nstep,0:nstep),conjg(transpose(GammaRet(0:nstep,0:nstep)))))*dt**2  +&
            matmul(dGret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),conjg(transpose(dGret(0:nstep,0:nstep)))))*dt**2
    endif
    return
  end subroutine obtain_Gimp
  !********************************************************************
  !********************************************************************
  !********************************************************************




















  !+-------------------------------------------------------------------+
  !PROGRAM  : SAVE_WF
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine save_solution()
    integer :: i,j,ik
    logical :: IOfile
    call system("if [ ! -d SEED ]; then mkdir SEED; fi")
    inquire(file="SEED/G0_t1t2.restart",exist=IOfile)
    if(IOfile)call system("rm -fv SEED/G0_t1t2.restart")
    inquire(file="SEED/G_t1t2.restart",exist=IOfile)
    if(IOfile)call system("rm -fv SEED/G_t1t2.restart")
    inquire(file="SEED/Sigma_t1t2.restart",exist=IOfile)
    if(IOfile)call system("rm -fv SEED/Sigma_t1t2.restart")
    inquire(file="SEED/nk.restart",exist=IOfile)
    if(IOfile)call system("rm -fv SEED/nk.restart")
    open(30,file="SEED/G0_t1t2.restart")
    open(31,file="SEED/G_t1t2.restart")
    open(32,file="SEED/Sigma_t1t2.restart")
    write(30,*)nstep
    do i=0,nstep
       do j=0,nstep
          write(30,*)G0less(i,j)
          write(31,*)locGless(i,j)
          write(32,*)Sless(i,j)
       enddo
    enddo
    write(30,*)"";write(31,*)"";write(32,*)""
    do i=0,nstep
       do j=0,nstep
          write(30,*)G0gtr(i,j)
          write(31,*)locGgtr(i,j)
          write(32,*)Sgtr(i,j)
       enddo
    enddo
    close(30);close(31);close(32)
    open(30,file="SEED/nk.restart",form="formatted")
    do i=0,nstep
       do ik=1,Lk
          write(30,*)nk(i,ik)
       enddo
    enddo
    close(30)
    call system("mv -vf SEED/G0_t1t2.restart neqWF.restart.new")
    return
  end subroutine save_solution
  !********************************************************************
  !********************************************************************
  !********************************************************************


















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
       w = wr(i)
       iw= cmplx(w,eps)
       fg=zero
       do im=1,Lmu
          fg=fg+1.d0/(iw-epsimu(im))/dble(Lmu)
       enddo
       A = -aimag(fg)/pi
       An= A*fermi(w,beta)
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
