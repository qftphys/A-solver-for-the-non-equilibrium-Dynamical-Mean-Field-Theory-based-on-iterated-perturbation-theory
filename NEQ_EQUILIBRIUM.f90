module NEQ_EQUILIBRIUM
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_INPUT_VARS  
  USE SF_CONSTANTS
  USE SF_IOTOOLS
  USE SF_SPECIAL
  USE DMFT_FFTGF
  USE DMFT_FFTAUX
  USE DMFT_MISC
  private

  interface neq_continue_equilibirum
     module procedure neq_continue_equilibirum_normal_lattice
     module procedure neq_continue_equilibirum_normal_bethe
     module procedure neq_continue_equilibirum_superc_lattice
     module procedure neq_continue_equilibirum_superc_bethe
  end interface neq_continue_equilibirum


  public  :: neq_continue_equilibirum     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.
  public  :: neq_setup_initial_conditions ! setup the initial conditions for the e/k dependent GF.


  ! public :: read_equilibrium_weiss
  ! public :: read_equilibrium_sigma
  ! public :: fft_extract_gtau

  real(8),public :: h1,h2,hDC,dens

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the NORMAL equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_normal_lattice(g0,gk,dgk,g,self,Hk,Wtk,params)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_gf)                 :: gk(:)
    type(kb_contour_dgf)                :: dgk(size(gk))
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: Hk(size(gk)),Wtk(size(gk))
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: SxG
    Lk=size(gk)
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.g%status)stop "init_functions: g is not allocated"
    do ik=1,Lk
       if(.not.gk(ik)%status)stop "init_functions: gk(ik) is not allocated"
       if(.not.dgk(ik)%status)stop "init_functions: dgk(ik) is not allocated"
    enddo
    if(.not.self%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_normal(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    ! ! !(this step depends on the imp. solv.)
    ! ! ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
    ! ! ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
    ! ! ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
    ! ! ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
    ! ! ! self^R(0,0) = self^> - self^<
    ! ! do j=0,Lf
    ! !    self%tau(j) = Ui*Ui*g0%tau(j)*g0%tau(Lf-j)*g0%tau(j)
    ! ! end do
    ! ! call extract_gtau_(self%tau,self%mats)
    ! ! ! Scoeff  = tail_coeff_sigma(Ui,0.5d0)
    ! ! call fft_sigma_tau2iw(self%iw,self%tau(0:),beta)!Scoeff)
    ! ! if(Ui==0d0)self%iw=zero
    ! ! Self%iw = xi*dimag(self%iw)                                   !imposing half-filling symmetry
    ! ! ! Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
    ! ! ! Self_gtr      =(xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
    ! ! ! Self%ret(1,1) = Self_gtr - Self%less(1,1)
    ! ! ! do j=0,L
    ! ! !    Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
    ! ! ! end do
    ! Scoeff  = tail_coeff_sigma(Ui,dens)
    ! call fft_sigma_iw2tau(self%iw,self%tau,beta,Scoeff)
    ! call extract_gtau_(self%tau,self%mats)
    ! if(Ui==0d0)self%tau=0d0
    ! if(Ui==0d0)self%mats=0d0
    ! self%less(1,1) = -xi*self%mats(L)                  !OK
    ! self%ret(1,1) =   xi*(self%mats(L)-self%mats(0))   !OK
    ! forall(i=0:L)self%lmix(1,i)=-xi*self%mats(L-i)     !small errors near 0,beta
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    call read_equilibrium_weiss_normal(g0,params,Hk=Hk,Wtk=Wtk)!<== get G0^{x=iw,tau,M,<,R,\lmix}
    ! Gcoeff = tail_coeff_glat(Ui,h1,h2,hDC)
    ! call fft_gf_iw2tau(g0%iw,g0%tau(0:),params%beta,Gcoeff)
    ! call extract_gtau_(g0%tau,g0%mats)
    ! g0%less(1,1) = -xi*g0%mats(L)
    ! g0%ret(1,1)  = -xi
    ! forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    !
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    G=zero
    do ik=1,Lk
       call neq_setup_initial_conditions_normal(gk(ik),dgk(ik),self,hk(ik),params)
       G%mats(0:)  = G%mats(0:)  + wtk(ik)*gk(ik)%mats(0:)
       G%tau(0:)   = G%tau(0:)   + wtk(ik)*gk(ik)%tau(0:)
       G%iw(:)     = G%iw(:)     + wtk(ik)*gk(ik)%iw(:)
       G%ret(1,1)  = G%ret(1,1)  + wtk(ik)*gk(ik)%ret(1,1)
       G%less(1,1) = G%less(1,1) + wtk(ik)*gk(ik)%less(1,1)
       G%lmix(1,0:)= G%lmix(1,0:)+ wtk(ik)*gk(ik)%lmix(1,0:)
    enddo
    return
  end subroutine neq_continue_equilibirum_normal_lattice

  subroutine neq_continue_equilibirum_normal_bethe(g0,dg0,g,self,params,wband)
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
    Lf= params%Niw
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_normal(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    ! !(this step depends on the imp. solv.)
    ! ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
    ! ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
    ! ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
    ! ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
    ! ! self^R(0,0) = self^> - self^<
    ! do j=0,L
    !    Self%mats(j) = Ui*Ui*g0%mats(j)*g0%mats(L-j)*g0%mats(j)
    ! end do
    ! Scoeff  = tail_coeff_sigma(Ui,0.5d0)
    ! call fft_sigma_tau2iw(Self%iw,Self%mats(0:),beta,Scoeff)
    ! Self%iw = xi*dimag(self%iw) !!ACTHUNG: imposing half-filling symmetry
    ! Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
    ! Self_gtr      =(xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
    ! do j=0,L
    !    Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
    ! end do
    ! Self%ret(1,1) = Self_gtr - Self%less(1,1)
    !
    !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
    call read_equilibrium_weiss_normal(g0,params,wband=wband)!<== get G0^{x=iw,tau,M,<,R,\lmix}
    ! !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    ! Gcoeff = tail_coeff_glat(U,0.5d0,0d0,0d0)
    ! call fft_gf_iw2tau(g0%iw,g0%mats(0:),params%beta,Gcoeff)
    ! g0%less(1,1) = -xi*g0%mats(L)
    ! g0%ret(1,1)  = -xi
    ! forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    do i=1,Lf
       wm    = pi/beta*dble(2*i-1)
       zeta  = xi*wm - Self%iw(i)
       G%iw(i) = gfbethe(wm,zeta,wband)
    enddo
    Gcoeff      = tail_coeff_glat(Ui,dens,xmu,hloc=0d0)
    call fft_gf_iw2tau(G%iw,G%tau,beta,Gcoeff)     !get G(tau)
    call fft_extract_gtau(G%tau,G%mats)
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
  end subroutine neq_continue_equilibirum_normal_bethe








  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the SUPERC equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_superc_lattice(g0,gk,gkaux,dgkaux,g,self,Hk,Wtk,params)
    type(kb_contour_gf)                 :: g0(2,2)
    type(kb_contour_gf)                 :: gk(:,:)               ![2][Lk]
    type(kb_contour_gf)                 :: gkaux(2,size(gk,2))   ![2][Lk]
    type(kb_contour_dgf)                :: dgkaux(2,size(gk,2))   ![2][Lk]
    type(kb_contour_gf)                 :: g(2,2)
    type(kb_contour_gf)                 :: self(2,2)
    real(8)                             :: Hk(size(gk,2)),Wtk(size(gk,2))
    type(kb_contour_params)             :: params
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: SxG
    !
    Lk=size(gk,2)
    do i=1,2
       do j=1,2
          if(.not.g0(i,j)%status)stop "neq_continue_equilibirum_superc_lattice error: g0 component is not allocated"
          if(.not.g(i,j)%status)stop "neq_continue_equilibirum_superc_lattice error: g component is not allocated"
          if(.not.self(i,j)%status)stop "neq_continue_equilibirum_superc_lattice error: self component is not allocated"
       enddo
    enddo
    do i=1,2
       do ik=1,Lk
          if(.not.gk(i,ik)%status)stop "neq_continue_equilibirum_superc_lattice error: gk component is not allocated"
          if(.not.gkaux(i,ik)%status)stop "neq_continue_equilibirum_superc_lattice error: gk_aux component is not allocated"
          if(.not.dgkaux(i,ik)%status)stop "neq_continue_equilibirum_superc_lattice error: Dgk_aux component is not allocated"
       enddo
    enddo
    if(.not.params%status)stop "neq_continue_equilibirum_superc_lattice error: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_superc(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE WEISS FIELD G0//F0^{x=M,<,R,\lmix}
    call read_equilibrium_weiss_superc(g0,params,Hk=Hk,Wtk=Wtk)!<== get G0^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    !<TODO
    ! this can be simplified by overloading the equality for kb_contour_gf.
    !>TODO
    do i=1,2
       do j=1,2
          g(i,j)=zero
       enddo
    enddo
    do ik=1,Lk
       call neq_setup_initial_conditions_superc(gk(:,ik),gkaux(:,ik),dgkaux(:,ik),self,hk(ik),params)
       G(1,1)%mats(0:)  = G(1,1)%mats(0:)  + wtk(ik)*gk(1,ik)%mats(0:)
       G(1,1)%tau(0:)   = G(1,1)%tau(0:)   + wtk(ik)*gk(1,ik)%tau(0:)
       G(1,1)%iw(:)     = G(1,1)%iw(:)     + wtk(ik)*gk(1,ik)%iw(:)
       G(1,1)%ret(1,1)  = G(1,1)%ret(1,1)  + wtk(ik)*gk(1,ik)%ret(1,1)
       G(1,1)%less(1,1) = G(1,1)%less(1,1) + wtk(ik)*gk(1,ik)%less(1,1)
       G(1,1)%lmix(1,0:)= G(1,1)%lmix(1,0:)+ wtk(ik)*gk(1,ik)%lmix(1,0:)
       !
       G(1,2)%mats(0:)  = G(1,2)%mats(0:)  + wtk(ik)*gk(2,ik)%mats(0:)
       G(1,2)%tau(0:)   = G(1,2)%tau(0:)   + wtk(ik)*gk(2,ik)%tau(0:)
       G(1,2)%iw(:)     = G(1,2)%iw(:)     + wtk(ik)*gk(2,ik)%iw(:)
       G(1,2)%ret(1,1)  = G(1,2)%ret(1,1)  + wtk(ik)*gk(2,ik)%ret(1,1)
       G(1,2)%less(1,1) = G(1,2)%less(1,1) + wtk(ik)*gk(2,ik)%less(1,1)
       G(1,2)%lmix(1,0:)= G(1,2)%lmix(1,0:)+ wtk(ik)*gk(2,ik)%lmix(1,0:)
       !
    enddo
    !
    !Get the other two components by symmetry
    call get_bar(G(2,2),G(1,1),params)
    call get_bar(G(2,1),G(1,2),params)
    !
    return
  end subroutine neq_continue_equilibirum_superc_lattice

  !<TODO
  ! Extend the Bethe lattice analytic case to  the superconducting channel.
  ! Sketch the algorithm for the self-consistency in this case and change 
  ! the corresponding NORMAL routine as for the lattice case.
  !>TODO








  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine neq_setup_initial_conditions_normal(Gk,dGk,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk,Self
    type(kb_contour_dgf)                :: dGk
    real(8)                             :: Hk
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    real(8)                             :: nk
    complex(8)                          :: epsk
    complex(8),allocatable,dimension(:) :: SxG
    real(8),dimension(:),allocatable    :: ftau
    L = params%Ntau
    Lf= params%Niw
    do i=1,Lf
       Gk%iw(i) = one/(xi*params%wm(i) + xmu - Hk - Self%iw(i))
    enddo
    call fft_gf_iw2tau(Gk%iw,Gk%tau(0:),beta)     !get G_k(tau)
    call fft_extract_gtau(Gk%tau,Gk%mats)            !
    Gk%less(1,1) = -xi*Gk%mats(L)                 !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk%ret(1,1)  = -xi                            !get G^R_k(0,0)=-xi
    forall(i=0:L)Gk%lmix(1,i)=-xi*Gk%mats(L-i)    !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    !Derivatives
    allocate(SxG(0:L))
    !get d/dt G_k^R = -i e(k,0)G_k^R
    dGk%ret(1)  = -xi*Hk*Gk%ret(1,1)
    !get d/dt G_k^< = -i e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    do k=0,L
       !      SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k))
       SxG(k)=self%lmix(1,k)*get_rmix(Gk,k,1,L)
    end do
    dGk%less(1) = -xi*Hk*Gk%less(1,1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
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
  end subroutine neq_setup_initial_conditions_normal


  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine neq_setup_initial_conditions_superc(Gk,Gkaux,dGkaux,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk(2),Gkaux(2),Self(2,2)
    type(kb_contour_dgf)                :: dGkaux(2)
    real(8)                             :: Hk
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    real(8)                             :: nk
    complex(8)                          :: epsk,zita
    complex(8),allocatable,dimension(:) :: SxG
    real(8),dimension(:),allocatable    :: ftau
    !
    L = params%Ntau
    Lf= params%Niw
    !
    !Construct the Gk functions
    do i=1,Lf
       zita = xi*params%wm(i) + xmu - Self(1,1)%iw(i)
       det  = abs(zita-Hk)**2 + Self(1,2)%iw(i)*Self(2,1)%iw(i)
       Gk(1)%iw(i) = (conjg(zita)-Hk)/det
       Gk(2)%iw(i) =-Self(2,1)%iw(i)/det
    enddo
    call fft_gf_iw2tau(Gk(1)%iw,Gk(1)%tau(0:),beta)     !get G_k(tau)
    call fft_gf_iw2tau(Gk(2)%iw,Gk(2)%tau(0:),beta)     !get F_k(tau)
    !ACTHUNG! it may be reasonable to perform the FT with C==0
    ! call fft_gf_iw2tau(Gk(2)%iw,Gk(2)%tau(0:),beta,C=[0d0,0d0,0d0,0d0])     !get F_k(tau)
    !
    call fft_extract_gtau(Gk(1)%tau,Gk(1)%mats)            !
    call fft_extract_gtau(Gk(2)%tau,Gk(2)%mats)            !
    !
    Gk(1)%less(1,1) = -xi*Gk(1)%mats(L)                 !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk(1)%ret(1,1)  = -xi                               !get G^R_k(0,0)=-xi
    forall(i=0:L)Gk(1)%lmix(1,i)=-xi*Gk(1)%mats(L-i)    !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    Gk(2)%less(1,1) = -xi*Gk(2)%mats(L)                 !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk(2)%ret(1,1)  =  zero                             !get G^R_k(0,0)=-xi
    forall(i=0:L)Gk(2)%lmix(1,i)=-xi*Gk(2)%mats(L-i)    !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    !Construct the auxiliary Gk functions== Gk_aux
    !Construct the derivatives of the auxiliary Gk functions == dGk_aux
    call neq_setup_initial_conditions_normal(Gkaux(1),dGkaux(1),self(1,1), hk,params)
    call neq_setup_initial_conditions_normal(Gkaux(2),dGkaux(2),self(2,2),-hk,params)
    !
  end subroutine neq_setup_initial_conditions_superc

















  !+-------------------------------------------------------------------+
  !PURPOSE: Read equilibrium solution and initialize the corresponding
  ! function.
  ! - read_equilibrium_weiss_normal: read and init Normal Weiss Field G0
  ! - read_equilibrium_sigma_normal: read and init Normal Self-energy function Self
  ! - read_equilibrium_weiss_superc: read and init Normal Weiss Field G0
  ! - read_equilibrium_sigma_superc: read and init Normal Self-energy function Self
  !+-------------------------------------------------------------------+
  subroutine read_equilibrium_sigma_normal(self,params)
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    real(8)                 :: wm,res,ims
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8)                 :: u_
    real(8),dimension(0:1)  :: Scoeff
    if(.not.self%status)  stop "read_equilibrium_sigma: sigma is not allocated"
    if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    inquire(file=trim(sigfile),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial Sigma(iw) from file "//reg(sigfile)
       unit = free_unit()
       i = file_length(trim(sigfile)) - 1
       open(unit,file=trim(sigfile),status='old')
       if(i/=Lf)then
          print*,"read_equilibrium_sigma: Liw in "//reg(sigfile)//" different from the input:",Lf
          print*,"read_equilibrium_sigma: check the header of the file u_i, <n>"
          stop
       endif
       read(unit,*)u_,dens
       if(u_/=Ui)then
          print*,"read_equilibrium_sigma: U_eq in "//reg(sigfile)//" different from the input:",u0
          stop
       endif
       write(*,"(3A)")"Header of the file:",reg(txtfy(u_)),reg(txtfy(dens))
       do i=1,Lf
          read(unit,*)wm,ims,res
          self%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Sigma(iw)=Ui*(n-1/2)"
       dens=0.5d0
       self%iw=Ui*(dens-0.5d0)    !set self-energy to zero or whatever is the initial HF solution
    endif
    Scoeff  = tail_coeff_sigma(Ui,dens)
    call fft_sigma_iw2tau(self%iw,self%tau(0:),beta,Scoeff)
    call fft_extract_gtau(self%tau(0:),self%mats(0:))
    if(Ui==0d0)self%tau=0d0
    if(Ui==0d0)self%mats=0d0
    self%less(1,1) = -xi*self%mats(L)                  !OK
    self%ret(1,1) =   xi*(self%mats(L)-self%mats(0))   !OK
    forall(i=0:L)self%lmix(1,i)=-xi*self%mats(L-i)     !small errors near 0,beta
  end subroutine read_equilibrium_sigma_normal

  subroutine read_equilibrium_weiss_normal(g0,params,wband,Hk,Wtk)
    type(kb_contour_gf)     :: g0
    type(kb_contour_params) :: params
    real(8),optional        :: Hk(:)
    real(8),optional        :: Wtk(:)
    real(8),optional        :: wband
    real(8)                 :: wband_
    real(8)                 :: wm,res,ims,mu
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8)                 :: Gcoeff(4)
    wband_=1d0;if(present(wband))wband_=wband
    !
    if(.not.g0%status)stop "read_equilibrium_g0: g0 is not allocated"
    if(.not.params%status)stop "read_equilibrium_g0: params is not allocated"
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
    inquire(file=trim(g0file),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial G0(iw) from file"//trim(g0file)
       unit = free_unit()
       i = file_length(trim(g0file))-1
       open(unit,file=trim(g0file),status='old')
       if(i/=Lf)then
          print*,"read_equilibrium_weiss: Liw in "//reg(g0file)//" different from the input:",Lf
          print*,"read_equilibrium_weiss: check the header of the file xmu, <h>, <h**2>, <h_dc>"
          stop
       endif
       read(unit,*)mu,h1,h2,hDC
       if(mu/=xmu)then
          print*,"read_equilibrium_weiss: mu in "//reg(g0file)//" different from the input:",xmu
          stop
       endif
       write(*,"(A)")"Header of the file:",reg(txtfy(mu)),reg(txtfy(h1)),reg(txtfy(h2)),reg(txtfy(hDC))
       do i=1,Lf
          read(unit,*)wm,ims,res
          g0%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       mu=xmu ; h1=0d0 ; h2=0d0 ; hDC=0d0
       write(*,"(5A)")"Header of the file:",reg(txtfy(mu)),reg(txtfy(h1)),reg(txtfy(h2)),reg(txtfy(hDC))
       if(present(Hk))then
          if(.not.present(Wtk))stop "read_equilibrium_weiss: Wtk not present"
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu
             g0%iw(i) = sum_overk_zeta(zeta,Hk,Wtk)
          enddo
       else
          write(*,"(A)")"Start from Non-interacting Bethe-model G0(iw)"
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu
             g0%iw(i) = gfbethe(wm,zeta,wband_)
          enddo
       endif
    endif
    !
    Gcoeff = tail_coeff_g0(Ui,h1,h2,hDC)
    call fft_gf_iw2tau(g0%iw,g0%tau(0:),params%beta,Gcoeff)
    call fft_extract_gtau(g0%tau,g0%mats)
    g0%less(1,1) = -xi*g0%mats(L)
    g0%ret(1,1)  = -xi
    forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
  end subroutine read_equilibrium_weiss_normal














  subroutine read_equilibrium_sigma_superc(self,params)
    type(kb_contour_gf)     :: self(2,2)
    type(kb_contour_params) :: params
    real(8)                 :: wm,res,ims
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8)                 :: u_
    real(8),dimension(0:1)  :: Scoeff
    do i=1,2
       do j=1,2
          if(.not.self(i,j)%status)  stop "read_equilibrium_sigma: sigma is not allocated"
       enddo
    enddo
    if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !NORMAL COMPONENT
    inquire(file=trim(sigfile),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial Sigma(iw) from file "//reg(sigfile)
       unit = free_unit()
       i = file_length(trim(sigfile))
       open(unit,file=trim(sigfile),status='old')
       if(i/=Lf)stop "read_equilibrium_sigma: Liw in "//reg(sigfile)//" different from the input:",Lf
       do i=1,Lf
          read(unit,*)wm,ims,res
          self(1,1)%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Sigma(iw)=Ui*(n-1/2)"
       dens=0.5d0
       self(1,1)%iw = Ui*(dens-0.5d0)    !set self-energy to zero or whatever is the initial HF solution
    endif
    !ANOMAL COMPONENT
    inquire(file=trim(selfile),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial Self(iw) from file "//reg(selfile)
       unit = free_unit()
       i = file_length(trim(selfile))
       open(unit,file=trim(selfile),status='old')
       if(i/=Lf)stop "read_equilibrium_sigma: Liw in "//reg(selfile)//" different from the input:",Lf
       do i=1,Lf
          read(unit,*)wm,ims,res
          self(1,2)%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Self(iw)=-Delta_SC"
       delta = deltasc
       self(1,2)%iw = -delta    !set self-energy to zero or whatever is the initial HF solution
    endif
    !
    call fft_sigma_iw2tau(self(1,1)%iw,self(1,1)%tau(0:),beta) !,C=[0d0,0d0])
    call fft_sigma_iw2tau(self(1,2)%iw,self(1,2)%tau(0:),beta) !,C=[0d0,0d0])
    call fft_extract_gtau(self(1,1)%tau(0:),self(1,1)%mats(0:))
    call fft_extract_gtau(self(1,2)%tau(0:),self(1,2)%mats(0:))
    if(Ui==0d0)self(1,1)%tau=0d0
    if(Ui==0d0)self(1,2)%tau=0d0
    if(Ui==0d0)self(1,1)%mats=0d0
    if(Ui==0d0)self(1,2)%mats=0d0
    !ACTHUNG BABY:
    self(1,1)%less(1,1) = -xi*self(1,1)%mats(L) !OK
    self(1,1)%ret(1,1) =   xi*(self(1,1)%mats(L)-self(1,1)%mats(0))   !OK
    forall(i=0:L)self(1,1)%lmix(1,i)=-xi*self(1,1)%mats(L-i)     !small errors near 0,beta
    !
    self(1,2)%less(1,1) = -xi*self(1,2)%mats(L) !OK
    self(1,2)%ret(1,1) =   xi*(self(1,2)%mats(L)-self(1,2)%mats(0))   !OK
    forall(i=0:L)self(1,2)%lmix(1,i)=-xi*self(1,2)%mats(L-i)     !small errors near 0,beta
    !
    call get_bar(self(2,2),self(1,1),params)
    call get_bar(self(2,1),self(1,2),params)
    return
  end subroutine read_equilibrium_sigma_superc

  subroutine read_equilibrium_weiss_superc(g0,params,wband,Hk,Wtk)
    type(kb_contour_gf)     :: g0(2,2)
    type(kb_contour_params) :: params
    real(8),optional        :: Hk(:)
    real(8),optional        :: Wtk(size(Hk))
    real(8),optional        :: wband
    real(8)                 :: wband_
    real(8)                 :: wm,res,ims,mu
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf,Lk
    real(8)                 :: Gcoeff(4)
    real(8),allocatable  :: Wt,Epsik
    !
    wband_=1d0;if(present(wband))wband_=wband
    !
    do i=1,2
       do j=1,2
          if(.not.g0(i,j)%status)stop "read_equilibrium_weiss_superc: g0 is not allocated"
       enddo
    enddo
    if(.not.params%status)stop "read_equilibrium_weiss_superc: params is not allocated"
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    Lk= 1000
    if(present(Hk))Lk=size(Hk)
    !
    !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
    !NORMAL COMPONENT
    inquire(file=trim(g0file),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial G0(iw) from file"//trim(g0file)
       unit = free_unit()
       i = file_length(trim(g0file))
       open(unit,file=trim(g0file),status='old')
       if(i/=Lf)stop "read_equilibrium_weiss: Liw in "//reg(g0file)//" different from the input:",Lf
       do i=1,Lf
          read(unit,*)wm,ims,res
          g0%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       if(present(Hk))then
          if(.not.present(Wtk))stop "read_equilibrium_weiss_superc: Wtk not present"
          g0(1,1)%iw = zero
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu
             do ik=1,Lk
                det   = abs(zeta-Hk(ik))**2 + deltasc**2
                g0(1,1)%iw(i) = g0(1,1)%iw(i) + Wtk(ik)*(conjg(zeta)-Hk(ik))/det
             enddo
          enddo
       else
          write(*,"(A)")"Start from Non-interacting Bethe-model G0(iw)"
          allocate(Wt(Lk),Epsik(Lk))
          call bethe_lattice(Wt,Epsik,Lk,wband_)
          g0(1,1)%iw = zero
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu
             do ik=1,Lk
                det   = abs(zeta-Epsik(ik))**2 + deltasc**2
                g0(1,1)%iw(i) = g0(1,1)%iw(i) + Wt(ik)*(conjg(zeta)-Epsik(ik))/det
             enddo
          enddo
          deallocate(Wt,Epsik)
       endif
    endif
    !ANOMAL COMPONENT
    inquire(file=trim(f0file),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial F0(iw) from file"//trim(f0file)
       unit = free_unit()
       i = file_length(trim(f0file))
       open(unit,file=trim(f0file),status='old')
       if(i/=Lf)stop "read_equilibrium_weiss: Liw in "//reg(f0file)//" different from the input:",Lf
       do i=1,Lf
          read(unit,*)wm,ims,res
          g0(1,2)%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       if(present(Hk))then
          if(.not.present(Wtk))stop "read_equilibrium_weiss_superc: Wtk not present"
          g0(1,2)%iw = zero
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu
             do ik=1,Lk
                det   = abs(zeta-Hk(ik))**2 + deltasc**2
                g0(1,2)%iw(i) = g0(1,2)%iw(i) + Wtk(ik)*deltasc/det
             enddo
          enddo
       else
          write(*,"(A)")"Start from Non-interacting Bethe-model G0(iw)"
          allocate(Wt(Lk),Epsik(Lk))
          call bethe_lattice(Wt,Epsik,Lk,wband_)
          g0(1,2)%iw = zero
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu
             do ik=1,Lk
                det   = abs(zeta-Epsik(ik))**2 + deltasc**2
                g0(1,2)%iw(i) = g0(1,2)%iw(i) + Wt(ik)*deltasc/det
             enddo
          enddo
          deallocate(Wt,Epsik)
       endif
    endif
    !
    !Continuation to the t=t`=0 (N=1)
    call fft_gf_iw2tau(g0(1,1)%iw,g0(1,1)%tau(0:),params%beta)
    call fft_gf_iw2tau(g0(1,2)%iw,g0(1,2)%tau(0:),params%beta)
    call fft_extract_gtau(g0(1,1)%tau,g0(1,1)%mats)
    call fft_extract_gtau(g0(1,2)%tau,g0(1,2)%mats)
    !
    g0(1,1)%less(1,1) = -xi*g0(1,1)%mats(L)
    g0(1,1)%ret(1,1)  = -xi
    forall(i=0:L)g0(1,1)%lmix(1,i)=-xi*g0(1,1)%mats(L-i)
    !
    g0(1,2)%less(1,1) = -xi*g0(1,2)%mats(L)
    g0(1,2)%ret(1,1)  = zero
    forall(i=0:L)g0(1,2)%lmix(1,i)=-xi*g0(1,2)%mats(L-i)
    !
    call get_bar(g0(2,2),g0(1,1),params)
    call get_bar(g0(2,1),g0(1,2),params)
    return
  end subroutine read_equilibrium_weiss_superc


























  subroutine fft_extract_gtau(g,gred)
    real(8),dimension(0:) :: g
    real(8),dimension(0:) :: gred
    integer               :: N,Nred
    integer               :: i,ip
    real(8)               :: p,mismatch
    N   =size(g)-1
    Nred=size(gred)-1
    gred(0)=g(0)
    gred(Nred)=g(N)
    mismatch=dble(N)/dble(Nred)
    do i=1,Nred-1
       p=dble(i)*mismatch
       ip=int(p)
       gred(i)=g(ip)
    enddo
  end subroutine fft_extract_gtau



  !   subroutine neq_continue_equilibirum__(g0,gk,dgk,g,self,Hk,Wtk,params,Kerk)
  !     type(kb_contour_gf)                 :: g0(2,2)
  !     type(kb_contour_gf)                 :: gk(:,:),gk_aux(size(gk,1),2),Kerk(size(gk,1),4)
  !     type(kb_contour_dgf)                :: dgk(size(gk,1),2)
  !     type(kb_contour_gf)                 :: g(2,2)
  !     type(kb_contour_gf)                 :: self(2,2)
  !     type(kb_contour_params)             :: params
  !     real(8)                             :: Hk(size(gk)),Wtk(size(gk))
  !     real(8)                             :: wm,res,ims
  !     logical                             :: bool
  !     integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
  !     complex(8)                          :: zeta
  !     complex(8)                          :: Self_gtr
  !     complex(8),allocatable,dimension(:) :: SxG
  !     Lk=size(gk,1)
  !     !   if(.not.g0%status)stop "init_functions: g0 is not allocated"
  !     !   if(.not.g%status)stop "init_functions: g is not allocated"
  !     !   do ik=1,Lk
  !     !      if(.not.gk(ik)%status)stop "init_functions: gk(ik) is not allocated"
  !     !      if(.not.dgk(ik)%status)stop "init_functions: dgk(ik) is not allocated"
  !     !   enddo
  !     !   if(.not.self%status)stop "init_functions: self is not allocated"
  !     !   if(.not.params%status)stop "init_functions: params is not allocated"
  !     !
  !     N = params%Nt
  !     L = params%Ntau
  !     Lf= params%Niw
  !     !
  !     !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
  !     !   call read_equilibrium_weiss(g0,params,Hk=Hk,Wtk=Wtk)!<== get G0^{x=iw,tau,M,<,R,\lmix}
  !     !
  !     !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
  !     !   call read_equilibrium_sigma(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
  !     !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
  !     !   G=zero
  !     do ik=1,Lk
  !        call neq_setup_initial_conditions(gk_aux(ik,1),dgk(ik,1),self(1,1),hk(ik),params)
  !        call neq_setup_initial_conditions(gk_aux(ik,2),dgk(ik,2),self(1,1),-hk(ik),params)

  !        call convolute_kb_contour_gf(Kerk(ik,1),gk_aux(ik,1),self(1,2),params)
  !        call convolute_kb_contour_gf(Kerk(ik,2),Kerk(ik,1),gk_aux(ik,2),params)
  !        call convolute_kb_contour_gf(Kerk(ik,3),Kerk(ik,2),self(2,1),params)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !        ! Now the subroutine convolute, if N==0, does the convolution also in
  !        ! Matsubara frequencies (a mere multiplication) and in Matsubara
  !        ! imaginary times. In this way we can use the equation satisfied by G and
  !        ! F_bar also to initialise them.

  !        do i=1,Lf
  !           Gk(ik,1)%iw(i) = gk_aux(ik,1)%iw(i)/(1d0-Kerk(ik,3)%iw(i))
  !        enddo
  !        call fft_gf_iw2tau(Gk(ik,1)%iw,Gk(ik,1)%tau(0:),beta)
  !        call fft_extract_gtau(Gk(ik,1)%tau,Gk(ik,1)%mats)

  !        call convolute_kb_contour_gf(Kerk(ik,4),gk_aux(ik,2),self(2,1),params)
  !        call convolute_kb_contour_gf(Gk(ik,2),Kerk(ik,4),gk(ik,1),params)

  !        ! Alternatively, one can do the initialisation 'by hand'

  !        do i=i,Lf
  !           Gk(ik,1)%iw(i)=-(xi*params%wm(i)-hk(ik))/&
  !                (params%wm(i)**2+(hk(ik)+self(1,1))**2+self(1,2)*self(2,1))
  !           Gk(ik,2)%iw(i)=-self(2,1)/(params%wm(i)**2+(hk(ik)+self(1,1))**2+self(1,2)*self(2,1))
  !        enddo
  !        call fft_gf_iw2tau(Gk(ik,1)%iw,Gk(ik,1)%tau(0:),beta)
  !        call fft_extract_gtau(Gk(ik,1)%tau,Gk(ik,1)%mats)
  !        call fft_gf_iw2tau(Gk(ik,2)%iw,Gk(ik,2)%tau(0:),beta)
  !        call fft_extract_gtau(Gk(ik,2)%tau,Gk(ik,2)%mats)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !        G(1,1)%mats(0:)  = G(1,1)%mats(0:)  + wtk(ik)*gk(ik,1)%mats(0:)
  !        G(1,1)%tau(0:)   = G(1,1)%tau(0:)   + wtk(ik)*gk(ik,1)%tau(0:)
  !        G(1,1)%iw(:)     = G(1,1)%iw(:)     + wtk(ik)*gk(ik,1)%iw(:)
  !        G(1,1)%ret(1,1)  = G(1,1)%ret(1,1)  + wtk(ik)*gk(ik,1)%ret(1,1)
  !        G(1,1)%less(1,1) = G(1,1)%less(1,1) + wtk(ik)*gk(ik,1)%less(1,1)
  !        G(1,1)%lmix(1,0:)= G(1,1)%lmix(1,0:)+ wtk(ik)*gk(ik,1)%lmix(1,0:)

  !        G(2,1)%mats(0:)  = G(2,1)%mats(0:)  + wtk(ik)*gk(ik,2)%mats(0:)
  !        G(2,1)%tau(0:)   = G(2,1)%tau(0:)   + wtk(ik)*gk(ik,2)%tau(0:)
  !        G(2,1)%iw(:)     = G(2,1)%iw(:)     + wtk(ik)*gk(ik,2)%iw(:)
  !        G(2,1)%ret(1,1)  = G(2,1)%ret(1,1)  + wtk(ik)*gk(ik,2)%ret(1,1)
  !        G(2,1)%less(1,1) = G(2,1)%less(1,1) + wtk(ik)*gk(ik,2)%less(1,1)
  !        G(2,1)%lmix(1,0:)= G(2,1)%lmix(1,0:)+ wtk(ik)*gk(ik,2)%lmix(1,0:)

  !        call get_bar(G(2,2),G(1,1),params)
  !        call get_bar(G(1,2),G(2,1),params)
  !     enddo
  !     return
  !   end subroutine neq_continue_equilibirum__

end module NEQ_EQUILIBRIUM
