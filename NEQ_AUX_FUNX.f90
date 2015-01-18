module NEQ_AUX_FUNX
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_VARS_GLOBAL  
  USE CONSTANTS
  USE IOTOOLS
  USE FUNCTIONS
  USE DMFT_TOOLS
  private

  interface neq_continue_equilibirum
     module procedure neq_continue_equilibirum_bethe!,neq_continue_equilibirum_
  end interface neq_continue_equilibirum

  public :: neq_continue_equilibirum     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.
  public :: neq_setup_weiss_field        ! continue the weiss field to the next time step.
  public :: neq_setup_initial_conditions ! setup the initial conditions for the e/k dependent GF.

contains



  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_bethe(g0,dg0,g,self,params,wband)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_dgf)                :: dg0
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8),optional                    :: wband
    real(8)                             :: wband_
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: GxG0
    real(8)                             :: Scoeff(2),Gcoeff(4)
    wband_=1d0;if(present(wband))wband_=wband
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.dg0%status)stop "init_functions: dg0 is not allocated"
    if(.not.g%status)stop "init_functions: g is not allocated"
    if(.not.self%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    N = params%Nt
    L = params%Ntau
    Lf= params%Lf
    !
    !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
    inquire(file=trim(g0file),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial G0(iw) from file "//reg(g0file)
       unit = free_unit()
       open(unit,file=reg(g0file),status='old')
       i = file_length(reg(g0file))
       if(i/=Lf)then
          print*,"init_equilibrium_weiss_field: Lfreq in +g0file does not correspond",i
          stop
       endif
       do i=1,Lf
          read(unit,*)wm,ims,res
          g0%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       write(*,"(A)")"Start from Non-interacting G0(iw)"
       do i=1,Lf
          wm    = pi/beta*dble(2*i-1)
          zeta  = xi*wm
          g0%iw(i) = gfbethe(wm,zeta,wband)
       enddo
    endif
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    Gcoeff = tail_coeff_glat(U,0.5d0,0d0,0d0)
    call fft_gf_iw2tau(g0%iw,g0%mats(0:),params%beta,Gcoeff)
    g0%less(1,1) = -xi*g0%mats(L)
    g0%ret(1,1)  = -xi
    forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    !(this step depends on the imp. solv.)
    ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
    ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
    ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
    ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
    ! self^R(0,0) = self^> - self^<
    do j=0,L
       Self%mats(j) = Ui*Ui*g0%mats(j)*g0%mats(L-j)*g0%mats(j)
    end do
    Scoeff  = tail_coeff_sigma(Ui,0.5d0)
    call fft_sigma_tau2iw(Self%iw,Self%mats(0:),beta,Scoeff)
    Self%iw = xi*dimag(self%iw) !!ACTHUNG: imposing half-filling symmetry
    Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
    Self_gtr      =(xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
    do j=0,L
       Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
    end do
    Self%ret(1,1) = Self_gtr - Self%less(1,1)
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    do i=1,Lf
       wm    = pi/beta*dble(2*i-1)
       zeta  = xi*wm - Self%iw(i)
       G%iw(i) = gfbethe(wm,zeta,wband)
    enddo
    Gcoeff      = tail_coeff_glat(U,0.5d0,0d0,0d0)
    call fft_gf_iw2tau(G%iw,G%mats,beta,Gcoeff)     !get G(tau)
    G%ret(1,1)  = -xi                               !get G^R(0,0)=-xi
    G%less(1,1) = -xi*G%mats(L)                  !get G^<(0,0)= xi*G^M(0-)
    forall(i=0:L)G%lmix(1,i)=-xi*G%mats(L-i) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
    !Derivatives
    allocate(GxG0(0:L))
    !
    !get d/dt G0^R = 0.d0
    dG0%ret(1)  = zero
    !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
    do k=0,L
       GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k))
    end do
    dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L) 
    !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
    dG0%lmix(0:)= zero
    do j=0,L
       do k=0,j
          GxG0(k)=G%lmix(1,k)*G0%mats(k+L-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
       do k=j,L
          GxG0(k)=G%lmix(1,k)*G0%mats(k-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
    enddo
    deallocate(GxG0)
    return
  end subroutine neq_continue_equilibirum_bethe






  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine neq_setup_initial_conditions(Gk,dGk,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk,Self
    type(kb_contour_dgf)                :: dGk
    real(8)                             :: Hk
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    complex(8),dimension(Lfreq)         :: Gkiw
    real(8),dimension(Lfreq)            :: wm
    real(8)                             :: nk
    complex(8)                          :: epsk
    complex(8),allocatable,dimension(:) :: SxG
    L = params%Ntau
    Lf= params%Lf
    do i=1,Lf
       Gk%iw(i) = one/(xi*params%wm(i) - Hk - Self%iw(i))
    enddo
    call fft_gf_iw2tau(Gk%iw,Gk%mats(0:),beta)        
    Gk%less(1,1) = -xi*Gk%mats(L)
    Gk%ret(1,1)  = -xi
    forall(i=0:L)Gk%lmix(1,i)=-xi*Gk%mats(L-i)
    !
    !Derivatives
    allocate(SxG(0:L))
    !get d/dt G_k^R = -i e(k,0)G_k^R
    dGk%ret(1)  = -xi*Hk*Gk%ret(1,1)
    !get d/dt G_k^< = -i e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    do k=0,L
       SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k))
    end do
    dGk%less(1) = -xi*Hk*Gk%less(1,1)-&
         xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
    !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(0:)= -xi*Hk*Gk%lmix(1,0:)
    do j=0,L
       do k=0,j
          SxG(k)=self%lmix(1,k)*Gk%mats(k+L-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       do k=j,L
          SxG(k)=self%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
    enddo
  end subroutine neq_setup_initial_conditions







  !+-------------------------------------------------------------------+
  !PURPOSE: setup the Weiss Field G0 for the next time-step
  !+-------------------------------------------------------------------+
  subroutine neq_setup_weiss_field(g0,params)
    type(kb_contour_gf)                   :: g0
    type(kb_contour_params)               :: params
    integer                               :: i,j,k,N,L
    complex(8),allocatable,dimension(:)   :: SxG
    if(.not.g0%status)stop "init_g0: g0 is not allocated"
    if(.not.params%status)stop "init_g0: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    !
    select case(N)
    case(1)
       return
    case(2)
       !GUESS G0 AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
       do j=1,N
          g0%ret(N,j) =g0%ret(1,1)
          g0%less(N,j)=g0%less(1,1)
       end do
       do i=1,N-1
          g0%less(i,N)=g0%less(1,1)
       end do
       do j=0,L
          g0%lmix(N,j)=g0%lmix(1,j)
       end do
    case default
       !EXTEND G0 FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
       !USING QUADRATIC EXTRAPOLATION
       do k=1,N-1
          g0%less(N,k)=2.d0*g0%less(N-1,k)-g0%less(N-2,k)
          g0%less(k,N)=2.d0*g0%less(k,N-1)-g0%less(k,N-2)
       end do
       g0%less(N,N)=2.d0*g0%less(N-1,N-1)-g0%less(N-2,N-2)
       !
       do k=0,L
          g0%lmix(N,k)=2.d0*g0%lmix(N-1,k)-g0%lmix(N-2,k)
       end do
       !
       g0%ret(N,N)=-xi
       do k=1,N-2
          g0%ret(N,k)=2.d0*g0%ret(N-1,k)-g0%ret(N-2,k)
       end do
       g0%ret(N,N-1)=0.5d0*(g0%ret(N,N)+g0%ret(N,N-2))
    end select
  end subroutine neq_setup_weiss_field

end module NEQ_AUX_FUNX
