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
     module procedure neq_continue_equilibirum_superc_lattice
     module procedure neq_continue_equilibirum_normal_bethe
  end interface neq_continue_equilibirum
  public  :: neq_continue_equilibirum     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.


  interface neq_setup_initial_conditions
     module procedure neq_setup_initial_conditions_normal
     module procedure neq_setup_initial_conditions_superc
  end interface neq_setup_initial_conditions
  public  :: neq_setup_initial_conditions ! setup the initial conditions for the e/k dependent GF.

contains




  !+-------------------------------------------------------------------+
  !PURPOSE: Read equilibrium solution and initialize the corresponding
  ! function.
  ! The code expects to read Sigma_reg(tau), that is the regular part 
  ! of the self-energy in imaginary time. Regular here means deprived 
  ! of the Hartree-Fock or Hartree-Fock-Bogoliubov term (1st order).
  ! The latter term is reconstructed analytically from the knowledge 
  ! of the static observables in the header of the self-energy files:
  ! NORMAL SIGMA:
  ! # U_i n [with U_i = interactino strenght at equilibrium, n=density
  ! ....
  ! ANOMALOUS SELF:
  ! # U_i Delta [with U_i as before, Delta = U*phi, phi=<cc> order parameter.
  ! ....
  !+-------------------------------------------------------------------+
  subroutine read_equilibrium_sigma_normal(self,params)
    type(kb_contour_sigma)           :: self
    type(kb_contour_params)          :: params
    real(8)                          :: tau
    logical                          :: bool
    complex(8)                       :: zeta
    integer                          :: i,j,k,ik,unit,Len,N,L,Lf,Ltau
    real(8)                          :: u_,dens
    real(8),dimension(0:1)           :: Scoeff
    real(8),allocatable,dimension(:) :: Stau !dummy function for FFT to tau
    !
    if(.not.self%status)  stop "read_equilibrium_sigma: sigma is not allocated"
    if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    inquire(file=reg(sigma_file),exist=bool)
    if(bool)then
       write(*,"(A)")"read_equilibrium_sigma_normal: reading Sigma(tau) from file "//reg(sigma_file)
       dens = read_header(reg(sigma_file))
       Len  = file_length(reg(sigma_file))
       Ltau = Len-1-1           !-1 for the header, -1 for the 0
       if(Ltau < L)stop "read_equilibrium_sigma_normal error: Ltau < params%Ntau"
       allocate(Stau(0:Ltau))
       unit = free_unit()
       open(unit,file=reg(sigma_file),status='old')
       do i=0,Ltau
          read(unit,*)tau,Stau(i)
       enddo
       close(unit)
       call fft_extract_gtau(Stau(0:),Self%reg%mats(0:))       !<=== Get Sigma(tau)= xtract(tmp_selt(tau_))
       self%hf(1) = Ui*(dens-0.5d0)                            !<=== Get Sigma_HF  = Re[Sigma(iw-->infty)]
       Scoeff  = tail_coeff_sigma(Ui,dens)
       call fft_sigma_tau2iw(Self%reg%iw,Stau(0:),beta,Scoeff) !<=== Get Sigma(iw) = fft(tmp_self(tau_))
       deallocate(Stau)
    else
       write(*,"(A)")"read_equilibrium_sigma_normal: start from Hartree-Fock Sigma(iw)=Ui*(n-1/2)"
       dens=0.5d0
       self%hf(1)   = Ui*(dens-0.5d0)
       self%reg%iw  = zero
       self%reg%mats= zero
    endif
    self%reg%less(1,1)  = -xi*self%reg%mats(L)                     !OK
    self%reg%ret(1,1)   =  xi*(self%reg%mats(0)+self%reg%mats(L))  !OK
    self%reg%lmix(1,0:L)= -xi*self%reg%mats(L:0:-1)                !small errors near 0,beta
  end subroutine read_equilibrium_sigma_normal

  subroutine read_equilibrium_sigma_superc(self,params)
    type(kb_contour_sigma)           :: Self(2,2)
    type(kb_contour_params)          :: params
    real(8)                          :: tau
    logical                          :: bool
    complex(8)                       :: zeta
    integer                          :: i,j,k,ik,unit,Len,N,L,Lf,Ltau
    real(8)                          :: u_,dens,delta
    real(8),dimension(0:1)           :: Scoeff
    real(8),allocatable,dimension(:) :: Stau !dummy function for FFT to tau
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    do i=1,2
       do j=1,2
          if(.not.Self(i,j)%status)  stop "read_equilibrium_sigma: Sigma is not allocated"
       enddo
    enddo
    if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
    !
    !
    !NORMAL COMPONENT
    inquire(file=reg(sigma_file),exist=bool)
    if(bool)then
       write(*,"(A)")"read_equilibrium_sigma_superc: reading Sigma(tau) from file "//reg(sigma_file)
       dens = read_header(reg(sigma_file))
       Len  = file_length(reg(sigma_file))
       Ltau = Len-1-1           !-1 for the header, -1 for the 0
       if(Ltau < L)stop "read_equilibrium_sigma_superc error: Ltau < params%Ntau"
       allocate(Stau(0:Ltau))
       unit = free_unit()
       open(unit,file=reg(sigma_file),status='old')
       do i=0,Ltau
          read(unit,*)tau,Stau(i)
       enddo
       close(unit)
       call fft_extract_gtau(Stau(0:),Self(1,1)%reg%mats(0:))       !<=== Get S_reg_11(tau)= xtract(tmp_self(tau_))
       self(1,1)%hf(1) = Ui*(dens-0.5d0)                            !<=== Get S_HF_11  = U*(n-1/2)
       self(2,2)%hf(1) =-Ui*(dens-0.5d0)                            !<=== Get S_HF_22  =-U*(n-1/2) = -conjg(S_HF_11)
       !<DEBUG 
       !Scoeff  = tail_coeff_sigma(Ui,dens)
       call fft_sigma_tau2iw(self(1,1)%reg%iw,Stau(0:),beta)!,Scoeff) !<=== Get S_reg_11(iw) = fft(tmp_self(tau_))
       !>DEBUG
       deallocate(Stau)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Sigma(iw)=Ui*(n-1/2)"
       dens=0.5d0
       self(1,1)%hf(1) = Ui*(dens-0.5d0)    !set self-energy to zero or whatever is the initial HF solution
       self(2,2)%hf(1) =-Ui*(dens-0.5d0)    !-conjg(selfHFB(1,1,1))
       self(1,1)%reg%iw   = zero
       self(1,1)%reg%mats = zero
    endif
    !
    !
    !ANOMAL COMPONENT
    inquire(file=trim(self_file),exist=bool)
    if(bool)then
       write(*,"(A)")"read_equilibrium_sigma_superc: reading Self(tau) from file "//reg(self_file)
       delta = read_header(reg(self_file))
       Len   = file_length(reg(self_file))
       Ltau  = Len-1-1           !-1 for the header, -1 for the 0
       if(Ltau < L)stop "read_equilibrium_sigma_superc error: Ltau < params%Ntau"
       allocate(Stau(0:Ltau))
       unit = free_unit()
       open(unit,file=reg(self_file),status='old')
       do i=0,Ltau
          read(unit,*)tau,Stau(i)
       enddo
       close(unit)
       call fft_extract_gtau(Stau(0:),Self(1,2)%reg%mats(0:)) !<=== Get S_reg_12(tau)= xtract(tmp_self(tau_))
       self(1,2)%hf(1) = -delta                               !<=== Get S_HF_12  = -Delta = -U*Phi
       self(2,1)%hf(1) = -delta                               !<=== Get S_HF_22  = -Delta = S_HF_12
       call fft_sigma_tau2iw(Self(1,2)%reg%iw,Stau(0:),beta)      !<=== Get S_reg_12(iw) = fft(tmp_self(tau_))
       deallocate(Stau)
    else
       !
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Self(iw)=-Delta_SC"
       delta = deltasc
       self(1,2)%hf(1) = -delta    !set self-energy to zero or whatever is the initial HF solution
       self(2,1)%hf(1) = -delta
       self(1,2)%reg%iw   = zero
       self(1,2)%reg%mats = zero
    endif
    !
    !
    self(1,1)%reg%less(1,1)  = -xi*self(1,1)%reg%mats(L)
    self(1,1)%reg%ret(1,1)   =  xi*(self(1,1)%reg%mats(L)+self(1,1)%reg%mats(0))
    self(1,1)%reg%lmix(1,0:L)= -xi*self(1,1)%reg%mats(L:0:-1)
    call get_bar(self(2,2),self(1,1),params)
    !
    self(1,2)%reg%less(1,1)  = -xi*self(1,2)%reg%mats(L)
    self(1,2)%reg%ret(1,1)   =  xi*(self(1,2)%reg%mats(L)+self(1,2)%reg%mats(0))
    self(1,2)%reg%lmix(1,0:L)=  xi*self(1,2)%reg%mats(L:0:-1)
    call get_bar(self(2,1),self(1,2),params)
    !
    return
  end subroutine read_equilibrium_sigma_superc








  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine neq_setup_initial_conditions_normal(Gk,dGk,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk
    type(kb_contour_dgf)                :: dGk
    type(kb_contour_sigma)              :: Self
    real(8)                             :: Hk
    real(8)                             :: Hk_hf
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    complex(8),allocatable,dimension(:) :: SxG
    real(8),allocatable,dimension(:)    :: Stau !dummy function for FFT to tau
    complex(8),allocatable,dimension(:) :: Siw  !dummy function in iw_n
    !
    L = params%Ntau
    Lf= params%Niw
    allocate(Siw(Lf),Stau(0:Lf))
    allocate(SxG(0:L))
    !
    Siw   = self%reg%iw + self%hf(1)
    Hk_hf = Hk - self%hf(1)
    do i=1,Lf
       Gk%iw(i) = one/(xi*params%wm(i) - Hk - Siw(i))
    enddo
    call fft_gf_iw2tau(Gk%iw,Stau(0:),beta) !get G_k(tau)
    call fft_extract_gtau(Stau,Gk%mats)     !
    Gk%less(1,1)  = -xi*Gk%mats(L)          !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk%ret(1,1)   = -xi                     !get G^R_k(0,0)=-xi
    Gk%lmix(1,0:L)= -xi*Gk%mats(L:0:-1)     !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    !Derivatives
    !get d/dt G_k^R = -i [e(k,0)-Sigma_HF(0)]G_k^R
    dGk%ret(1)  = -xi*Hk_hf*Gk%ret(1,1)
    !
    !get d/dt G_k^< = -i [e(k,0)-Sigma_HF(0)]G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    SxG(0:)=zero
    do k=0,L
       SxG(k)=self%reg%lmix(1,k)*conjg(Gk%lmix(1,L-k))
    end do
    dGk%less(1) = -xi*Hk_hf*Gk%less(1,1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
    !
    !get d/dt G_k^\lmix = -xi*[e(k,0)-Sigma_HF(0)]*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(0:)= -xi*Hk_hf*Gk%lmix(1,0:)
    do j=0,L
       SxG(0:)=zero
       do k=0,j
          SxG(k)=Self%reg%lmix(1,k)*Gk%mats(k+L-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       SxG(0:)=zero
       do k=j,L
          SxG(k)=Self%reg%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
    enddo
    !
    deallocate(Siw,Stau,SxG)
  end subroutine neq_setup_initial_conditions_normal


  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine neq_setup_initial_conditions_superc(Gk,Gkaux,dGkaux,Self,Hk,params)
    type(kb_contour_gf)                   :: Gk(2)
    type(kb_contour_gf)                   :: Gkaux(2)
    type(kb_contour_dgf)                  :: dGkaux(2)
    type(kb_contour_sigma)                :: Self(2,2)
    real(8)                               :: Hk
    type(kb_contour_params)               :: params
    integer                               :: i,j,k,L,Lf
    real(8)                               :: nk,det
    complex(8)                            :: epsk,zita
    complex(8),allocatable,dimension(:)   :: SxG
    real(8),allocatable,dimension(:,:)    :: Stau !dummy function for FFT to tau
    complex(8),allocatable,dimension(:,:) :: Siw  !dummy function in iw_n
    !
    L = params%Ntau
    Lf= params%Niw
    allocate(Siw(2,Lf),Stau(2,0:Lf))
    !
    !Construct the Gk functions: Sigma(iw) here includes HF correction: Sigma=Sigma_reg + Sigma_HF
    Siw(1,:) = Self(1,1)%reg%iw + Self(1,1)%hf
    Siw(2,:) = Self(1,2)%reg%iw + Self(1,2)%hf
    do i=1,Lf
       zita = xi*params%wm(i) + xmu - Siw(1,i)
       det  = abs(zita-Hk)**2 + Siw(2,i)**2
       Gk(1)%iw(i) = (conjg(zita)-Hk)/det
       Gk(2)%iw(i) =-Siw(2,i)/det
    enddo
    !get G_k(tau)
    call fft_gf_iw2tau(Gk(1)%iw,Stau(1,0:),beta)
    call fft_extract_gtau(Stau(1,0:),Gk(1)%mats(0:))
    !get F_k(tau)
    call fft_gf_iw2tau(Gk(2)%iw,Stau(2,0:),beta,C=[0d0,0d0,0d0,0d0])
    call fft_extract_gtau(Stau(2,0:),Gk(2)%mats(0:))
    !
    Gk(1)%less(1,1)  = -xi*Gk(1)%mats(L)      !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk(1)%ret(1,1)   = -xi                    !get G^R_k(0,0)=-xi
    Gk(1)%lmix(1,0:L)= -xi*Gk(1)%mats(L:0:-1) !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    Gk(2)%less(1,1)  = -xi*Gk(2)%mats(L)      !get F^<_k(0,0)= xi*F^M_k(0-)
    Gk(2)%ret(1,1)   =  zero                  !get F^R_k(0,0)=0
    Gk(2)%lmix(1,0:L)=  xi*Gk(2)%mats(L:0:-1) !get F^\lmix_k(0,tau)=xi*F_k(tau<0)= xi*F_k(beta-tau>0)
    !
    !Construct the auxiliary Gk functions== Gk_aux
    !Construct the derivatives of the auxiliary Gk functions == dGk_aux
    call neq_setup_initial_conditions_normal(Gkaux(1),dGkaux(1),Self(1,1), hk,params)
    call neq_setup_initial_conditions_normal(Gkaux(2),dGkaux(2),Self(2,2),-hk,params)
    !
  end subroutine neq_setup_initial_conditions_superc





  !+-------------------------------------------------------------------+
  !PURPOSE: continue the NORMAL/SUPERC equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_normal_lattice(g0,gk,dgk,g,self,Hk,Wtk,params)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_gf)                 :: gk(:)
    type(kb_contour_dgf)                :: dgk(size(gk))
    type(kb_contour_gf)                 :: g
    type(kb_contour_sigma)              :: self
    real(8)                             :: Hk(size(gk)),Wtk(size(gk))
    type(kb_contour_params)             :: params
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
    real(8),allocatable,dimension(:)    :: G0tau !dummy function for FFT to tau
    complex(8),allocatable,dimension(:) :: Siw  !dummy function in iw_n
    !
    Lk=size(gk)
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    allocate( Siw(Lf),G0tau(0:Lf) )
    !
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.g%status)stop "init_functions: g is not allocated"
    if(.not.self%status)stop "init_functions: self is not allocated"
    do ik=1,Lk
       if(.not.gk(ik)%status)stop "init_functions: gk(ik) is not allocated"
       if(.not.dgk(ik)%status)stop "init_functions: dgk(ik) is not allocated"
    enddo
    if(.not.params%status)stop "init_functions: params is not allocated"
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_normal(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    g=zero
    do ik=1,Lk
       call neq_setup_initial_conditions_normal(gk(ik),dgk(ik),self,hk(ik),params)
    enddo
    call sum_kb_contour_gf(gk(:),wtk(:),g,params)
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    Siw   = self%reg%iw + self%hf(1)       !Sigma = Sigma_reg + SigmaHF
    ! g0%iw = one/( one/g%iw + Siw )       !Dyson: G0^-1 = Gloc^-1 + Sigma
    g0%iw = one/( one/g%iw + self%reg%iw ) !Dyson: tildeG0^-1 = Gloc^-1 + Sigma_reg = Gloc^-1 + Sigma - Sigma_HF
    call fft_gf_iw2tau(g0%iw,G0tau(0:),params%beta)
    call fft_extract_gtau(G0tau,g0%mats)
    g0%less(1,1)  = -xi*g0%mats(L)
    g0%ret(1,1)   = -xi
    g0%lmix(1,0:L)= -xi*g0%mats(L:0:-1)
    !
    deallocate(G0tau,Siw)
  end subroutine neq_continue_equilibirum_normal_lattice

  subroutine neq_continue_equilibirum_superc_lattice(g0,gk,gkaux,dgkaux,g,self,Hk,Wtk,params)
    type(kb_contour_gf)                   :: g0(2,2)
    type(kb_contour_gf)                   :: gk(:,:)               ![2][Lk]
    type(kb_contour_gf)                   :: gkaux(2,size(gk,2))   ![2][Lk]
    type(kb_contour_dgf)                  :: dgkaux(2,size(gk,2))   ![2][Lk]
    type(kb_contour_gf)                   :: g(2,2)
    type(kb_contour_sigma)                :: self(2,2)
    real(8)                               :: Hk(size(gk,2)),Wtk(size(gk,2))
    type(kb_contour_params)               :: params
    logical                               :: bool
    integer                               :: i,j,k,ik,unit,len,N,L,Lf,Lk
    real(8),allocatable,dimension(:,:)    :: G0tau !dummy function for FFT to tau
    complex(8),allocatable,dimension(:,:) :: Siw  !dummy function in iw_n
    real(8),allocatable,dimension(:)      :: Gdet
    !
    Lk=size(gk,2)
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    allocate( Siw(2,Lf), G0tau(2,0:Lf), Gdet(Lf) )
    !
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
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_superc(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    do i=1,2
       do j=1,2
          G(i,j)=zero
       enddo
    enddo
    do ik=1,Lk
       call neq_setup_initial_conditions_superc(gk(:,ik),gkaux(:,ik),dgkaux(:,ik),self,hk(ik),params)
    enddo
    call sum_kb_contour_gf(gk(1,:),wtk(:),g(1,1),params); call get_bar(g(2,2),g(1,1),params)
    call sum_kb_contour_gf(gk(2,:),wtk(:),g(2,1),params); call get_bar(g(1,2),g(2,1),params)
    !
    !<DEBUG 
    print*,Ui*fft_get_density(g(1,2)%iw,beta,[0d0,0d0,0d0,0d0])
    !>DEBUG
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    ! get the elements of the inverse matrix
    ! tildeG0^-1 = Gloc^-1 + Sigma_Reg==Sigma-SigmaHFB
    Gdet       = abs(g(1,1)%iw)**2 + g(1,2)%iw**2                            !== g(1,1)*g(2,2) + g(1,2)*g(2,1)
    g0(1,1)%iw = conjg(g(1,1)%iw)/Gdet + self(1,1)%reg%iw
    g0(1,2)%iw = g(1,2)%iw/Gdet        + self(1,2)%reg%iw
    !get the elements of G0 matrix by inverting again:
    Gdet       = abs(g0(1,1)%iw)**2 + (g0(1,2)%iw)**2
    g0(1,1)%iw = conjg(g0(1,1)%iw)/Gdet
    g0(1,2)%iw = g0(1,2)%iw/Gdet
    !FFt
    call fft_gf_iw2tau(g0(1,1)%iw,G0tau(1,0:),params%beta)
    call fft_extract_gtau(G0tau(1,0:),g0(1,1)%mats)
    call fft_gf_iw2tau(g0(1,2)%iw,G0tau(2,0:),params%beta,[0d0,0d0,0d0,0d0])
    call fft_extract_gtau(G0tau(2,0:),g0(1,2)%mats)
    g0(1,1)%less(1,1)   = -xi*g0(1,1)%mats(L)
    g0(1,1)%ret(1,1)    = -xi
    g0(1,1)%lmix(1,0:L) = -xi*g0(1,1)%mats(L:0:-1)
    call get_bar(g0(2,2),g0(1,1),params)
    !
    g0(1,2)%less(1,1)   = -xi*g0(1,2)%mats(L)
    g0(1,2)%ret(1,1)    =  zero
    g0(1,2)%lmix(1,0:L) =  xi*g0(1,2)%mats(L:0:-1)
    call get_bar(g0(2,1),g0(1,2),params)
    !
    deallocate(Siw,G0tau,Gdet)
    return
  end subroutine neq_continue_equilibirum_superc_lattice



  subroutine neq_continue_equilibirum_normal_bethe(g0,dg0,g,self,params,wband)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_dgf)                :: dg0
    type(kb_contour_gf)                 :: g
    type(kb_contour_sigma)              :: self
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
    real(8),allocatable,dimension(:)    :: Gtau
    complex(8),allocatable,dimension(:) :: Siw

    !
    wband_=1d0;if(present(wband))wband_=wband
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.dg0%status)stop "init_functions: dg0 is not allocated"
    if(.not.g%status)stop "init_functions: g is not allocated"
    if(.not.self%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    allocate( Siw(Lf), Gtau(0:L), GxG0(0:L) )
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_normal(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    Siw = self%reg%iw + self%hf(1)
    do i=1,Lf
       wm    = pi/beta*dble(2*i-1)
       zeta  = xi*wm - Siw(i)
       G%iw(i) = gfbethe(wm,zeta,wband)
    enddo
    !Gcoeff      = tail_coeff_glat(Ui,dens,xmu,hloc=0d0)
    call fft_gf_iw2tau(G%iw,Gtau(0:),beta,Gcoeff)     !get G(tau)
    call fft_extract_gtau(Gtau,G%mats)
    G%ret(1,1)   = -xi                !get G^R(0,0)=-xi
    G%less(1,1)  = -xi*G%mats(L)      !get G^<(0,0)= xi*G^M(0-)
    G%lmix(1,0:L)= -xi*G%mats(L:0:-1) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    g0%iw = one/(one/g%iw + self%reg%iw) !Dyson: tideG0^-1 = Gloc^-1 + Sigma_reg
    call fft_gf_iw2tau(g0%iw,Gtau(0:),params%beta)
    call fft_extract_gtau(Gtau,g0%mats)
    g0%less(1,1)  = -xi*g0%mats(L)
    g0%ret(1,1)   = -xi
    g0%lmix(1,0:L)= -xi*g0%mats(L:0:-1)
    !
    !Derivatives d/dt G0:
    !get d/dt G0^R = 0.d0
    dG0%ret(1)  = zero
    !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
    GxG0(0:)=zero
    do k=0,L
       GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k))
    end do
    dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L) 
    !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
    dG0%lmix(0:)= zero

    do j=0,L
       GxG0(0:)=zero
       do k=0,j
          GxG0(k)=G%lmix(1,k)*G0%mats(k+L-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
       GxG0(0:)=zero
       do k=j,L
          GxG0(k)=G%lmix(1,k)*G0%mats(k-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
    enddo
    !
    deallocate(Siw,Gtau,GxG0)
    return
  end subroutine neq_continue_equilibirum_normal_bethe
  !<TODO
  ! Extend the Bethe lattice analytic case to  the superconducting channel.
  ! Sketch the algorithm for the self-consistency in this case and change 
  ! the corresponding NORMAL routine as for the lattice case.
  !>TODO


  




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: AUXILIARY ROUTINES:
  !  - read_header     : read the header of the seed file for the Self-energy S(tau)
  !  - fft_extract_gtau: extract a G(tau) from a dense to a sparse set.
  !+-----------------------------------------------------------------------------+!
  function read_header(file)  result(obs)
    character(len=*) :: file
    integer          :: unit
    real(8)          :: u_,obs
    unit = free_unit()
    open(unit,file=file,status='old')
    read(unit,*)u_,obs
    close(unit)
    if(u_/=Ui)then
       print*,"read_header error: U_eq in "//file//" different from the input:",Ui
       stop
    endif
    write(*,"(4A)")"Header of the file:",reg(txtfy(u_))," ",reg(txtfy(obs))
  end function read_header


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



end module NEQ_EQUILIBRIUM
