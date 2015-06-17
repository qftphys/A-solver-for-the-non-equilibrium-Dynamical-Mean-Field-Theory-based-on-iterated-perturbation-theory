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
  end interface neq_continue_equilibirum

  interface neq_setup_initial_conditions
     module procedure neq_setup_initial_conditions_normal
     module procedure neq_setup_initial_conditions_superc
  end interface neq_setup_initial_conditions

  public  :: neq_continue_equilibirum     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.
  public  :: neq_setup_initial_conditions ! setup the initial conditions for the e/k dependent GF.


  real(8)                          :: h1,h2,hDC,dens
  real(8),allocatable,dimension(:) :: Rtau

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the NORMAL equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_normal_lattice(g0,gk,dgk,g,selfHF,selfReg,self,Hk,Wtk,params)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_gf)                 :: gk(:)
    type(kb_contour_dgf)                :: dgk(size(gk))
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_gf)                 :: selfReg
    complex(8)                          :: selfHF(:)
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
    if(size(selfHF)/=params%Ntime)stop "init_functions: size(selfHF)!=Ntime"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_normal(selfHF,selfReg,self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    G=zero
    do ik=1,Lk
       call neq_setup_initial_conditions_normal(gk(ik),dgk(ik),selfHF(1),selfReg,self,hk(ik),params)
    enddo
    call sum_kb_contour_gf(gk(:),wtk(:),g,params)
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    g0%iw = one/(one/g%iw + self%iw) !Dyson: G0^-1 = Gloc^-1 + Sigma
    allocate(Rtau(0:Lf))
    call fft_gf_iw2tau(g0%iw,Rtau(0:),params%beta)
    call fft_extract_gtau(Rtau,g0%mats)
    deallocate(Rtau)
    g0%less(1,1) = -xi*g0%mats(L)
    g0%ret(1,1)  = -xi
    forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    return
  end subroutine neq_continue_equilibirum_normal_lattice

  subroutine neq_continue_equilibirum_normal_bethe(g0,dg0,g,selfHF,selfReg,self,params,wband)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_dgf)                :: dg0
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_gf)                 :: selfReg
    complex(8)                          :: selfHF(:)
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
    if(size(selfHF)/=params%Ntime)stop "init_functions: size(selfHF)!=Ntime"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_normal(selfHF,selfReg,self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    do i=1,Lf
       wm    = pi/beta*dble(2*i-1)
       zeta  = xi*wm - Self%iw(i)
       G%iw(i) = gfbethe(wm,zeta,wband)
    enddo
    !Gcoeff      = tail_coeff_glat(Ui,dens,xmu,hloc=0d0)
    allocate(Rtau(0:Lf))
    call fft_gf_iw2tau(G%iw,Rtau(0:),beta,Gcoeff)     !get G(tau)
    call fft_extract_gtau(Rtau(0:),G%mats)
    deallocate(Rtau)
    G%ret(1,1)  = -xi                               !get G^R(0,0)=-xi
    G%less(1,1) = -xi*G%mats(L)                  !get G^<(0,0)= xi*G^M(0-)
    forall(i=0:L)G%lmix(1,i)=-xi*G%mats(L-i) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    g0%iw = one/(one/g%iw + self%iw) !Dyson: G0^-1 = Gloc^-1 + Sigma
    allocate(Rtau(0:Lf))
    call fft_gf_iw2tau(g0%iw,Rtau(0:),params%beta)
    call fft_extract_gtau(Rtau,g0%mats)
    deallocate(Rtau)
    g0%less(1,1) = -xi*g0%mats(L)
    g0%ret(1,1)  = -xi
    forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    !
    !Derivatives d/dt G0:
    allocate(GxG0(0:L))
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
    !
    return
  end subroutine neq_continue_equilibirum_normal_bethe


  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the SUPERC equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_superc_lattice(g0,gk,gkaux,dgkaux,g,selfHF,selfReg,self,Hk,Wtk,params)
    type(kb_contour_gf)                 :: g0(2,2)
    type(kb_contour_gf)                 :: gk(:,:)               ![2][Lk]
    type(kb_contour_gf)                 :: gkaux(2,size(gk,2))   ![2][Lk]
    type(kb_contour_dgf)                :: dgkaux(2,size(gk,2))   ![2][Lk]
    type(kb_contour_gf)                 :: g(2,2)
    type(kb_contour_gf)                 :: self(2,2)
    type(kb_contour_gf)                 :: selfReg(2,2)
    complex(8)                          :: selfHF(:,:,:)
    real(8)                             :: Hk(size(gk,2)),Wtk(size(gk,2))
    type(kb_contour_params)             :: params
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
    real(8),allocatable,dimension(:)    :: det
    !
    Lk=size(gk,2)
    do i=1,2
       do j=1,2
          if(.not.g0(i,j)%status)stop "neq_continue_equilibirum_superc_lattice error: g0 component is not allocated"
          if(.not.g(i,j)%status)stop "neq_continue_equilibirum_superc_lattice error: g component is not allocated"
          if(.not.self(i,j)%status)stop "neq_continue_equilibirum_superc_lattice error: self component is not allocated"
          if(.not.selfReg(i,j)%status)stop "neq_continue_equilibirum_superc_lattice error: selfReg component is not allocated"
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
    if(size(selfHF,1)/=2)stop "init_functions: size(selfHF,1)!=2"
    if(size(selfHF,2)/=2)stop "init_functions: size(selfHF,2)!=2"
    if(size(selfHF,3)/=params%Ntime)stop "init_functions: size(selfHF,3)!=Ntime"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma_superc(selfHF,selfReg,self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    do i=1,2
       do j=1,2
          G(i,j)=zero
       enddo
    enddo
    do ik=1,Lk
       call neq_setup_initial_conditions_superc(gk(:,ik),gkaux(:,ik),dgkaux(:,ik),selfHF,selfReg,self,hk(ik),params)
    enddo
    call sum_kb_contour_gf(gk(1,:),wtk(:),G(1,1),params)
    call sum_kb_contour_gf(gk(2,:),wtk(:),G(2,1),params)
    !
    !Get the other two components by symmetry
    call get_bar(G(2,2),G(1,1),params)
    call get_bar(G(1,2),G(2,1),params)
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    allocate(det(Lf))
    det = abs(g(1,1)%iw)**2 + g(1,2)%iw**2 !== g(1,1)*g(2,2) + g(1,2)*g(2,1)
    g0(1,1)%iw = conjg(g(1,1)%iw)/det + selfReg(1,1)%iw
    g0(1,2)%iw =-g(1,2)%iw/det        + selfReg(1,2)%iw
    allocate(Rtau(0:Lf))
    call fft_gf_iw2tau(g0(1,1)%iw,Rtau(0:),params%beta)
    call fft_extract_gtau(Rtau,g0(1,1)%mats)
    call fft_gf_iw2tau(g0(1,2)%iw,Rtau(0:),params%beta,[0d0,0d0,0d0,0d0])
    call fft_extract_gtau(Rtau,g0(1,2)%mats)
    deallocate(Rtau)
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
  subroutine neq_setup_initial_conditions_normal(Gk,dGk,SelfHF,SelfReg,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk,SelfReg,Self
    type(kb_contour_dgf)                :: dGk
    complex(8)                          :: SelfHF
    real(8)                             :: Hk
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    real(8)                             :: nk
    complex(8)                          :: epsk
    complex(8),allocatable,dimension(:) :: SxG
    real(8),dimension(:),allocatable    :: ftau
    !
    L = params%Ntau
    Lf= params%Niw
    !
    do i=1,Lf
       Gk%iw(i) = one/(xi*params%wm(i) + xmu - Hk - Self%iw(i))
    enddo
    allocate(Rtau(0:Lf))
    call fft_gf_iw2tau(Gk%iw,Rtau(0:),beta)     !get G_k(tau)
    call fft_extract_gtau(Rtau,Gk%mats)            !
    deallocate(Rtau)
    Gk%less(1,1) = -xi*Gk%mats(L)                 !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk%ret(1,1)  = -xi                            !get G^R_k(0,0)=-xi
    forall(i=0:L)Gk%lmix(1,i)=-xi*Gk%mats(L-i)    !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    !Derivatives
    allocate(SxG(0:L))
    !get d/dt G_k^R = -i [e(k,0)+Sigma_HF(0)]G_k^R
    dGk%ret(1)  = -xi*(Hk+SelfHF)*Gk%ret(1,1)
    !get d/dt G_k^< = -i [e(k,0)+Sigma_HF(0)]G_k^< -xi(-xi)int_0^beta Sreg^\lmix*G_k^\rmix
    do k=0,L
       SxG(k)=SelfReg%lmix(1,k)*conjg(Gk%lmix(1,L-k))
    end do
    dGk%less(1) = -xi*(Hk+SelfHF)*Gk%less(1,1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
    !get d/dt G_k^\lmix = -xi*[e(k,0)+Sigma_HF(0)]*G_k^\lmix - xi*int_0^beta Sreg^\lmix*G_k^M
    dGk%lmix(0:)= -xi*(Hk+SelfHF)*Gk%lmix(1,0:)
    do j=0,L
       do k=0,j
          SxG(k)=SelfReg%lmix(1,k)*Gk%mats(k+L-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       do k=j,L
          SxG(k)=SelfReg%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
    enddo
  end subroutine neq_setup_initial_conditions_normal


  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine neq_setup_initial_conditions_superc(Gk,Gkaux,dGkaux,SelfHFB,SelfReg,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk(2),Gkaux(2),SelfReg(2,2),Self(2,2)
    type(kb_contour_dgf)                :: dGkaux(2)
    real(8)                             :: Hk
    complex(8)                          :: SelfHFB(2,2)
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    real(8)                             :: nk,det
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
    allocate(Rtau(0:Lf))
    !get G_k(tau)
    call fft_gf_iw2tau(Gk(1)%iw,Rtau(0:),beta)
    call fft_extract_gtau(Rtau,Gk(1)%mats)
    !get F_k(tau)
    call fft_gf_iw2tau(Gk(2)%iw,Rtau(0:),beta)!,C=[0d0,0d0,0d0,0d0])
    call fft_extract_gtau(Rtau,Gk(2)%mats)
    deallocate(Rtau)
    !
    Gk(1)%less(1,1) = -xi*Gk(1)%mats(L)                 !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk(1)%ret(1,1)  = -xi                               !get G^R_k(0,0)=-xi
    forall(i=0:L)Gk(1)%lmix(1,i)=-xi*Gk(1)%mats(L-i)    !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    Gk(2)%less(1,1) = -xi*Gk(2)%mats(L)                 !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk(2)%ret(1,1)  =  zero                             !get G^R_k(0,0)=0
    forall(i=0:L)Gk(2)%lmix(1,i)=-xi*Gk(2)%mats(L-i)    !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    !Construct the auxiliary Gk functions== Gk_aux
    !Construct the derivatives of the auxiliary Gk functions == dGk_aux
    call neq_setup_initial_conditions_normal(Gkaux(1),dGkaux(1),selfHFB(1,1),SelfReg(1,1),Self(1,1), hk,params)
    call neq_setup_initial_conditions_normal(Gkaux(2),dGkaux(2),selfHFB(2,2),SelfReg(2,2),self(2,2),-hk,params)
    !
  end subroutine neq_setup_initial_conditions_superc

















  !+-------------------------------------------------------------------+
  !PURPOSE: Read equilibrium solution and initialize the corresponding
  ! function.
  ! We assume that the read self-energy function INCLUDES the Hartree-Fock
  ! term, so we are able to extract it and separate regular from local part.
  ! - read_equilibrium_sigma: read and init Self-energy function Self
  ! - read_equilibrium_sigma_superc: read and init Normal Self-energy function Self
  !+-------------------------------------------------------------------+
  subroutine read_equilibrium_sigma_normal(selfHF,selfReg,self,params)
    complex(8)              :: selfHF(:)
    type(kb_contour_gf)     :: selfReg
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    real(8)                 :: wm,res,ims
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8)                 :: u_
    real(8),dimension(0:1)  :: Scoeff
    !
    if(.not.self%status)  stop "read_equilibrium_sigma: sigma is not allocated"
    if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
    if(size(selfHF)/=params%Ntime)stop "read_equilibrium_sigma: size(selfHF)!=Ntime"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    selfHF=zero
    !
    inquire(file=trim(sigfile),exist=bool)
    if(bool)then
       !
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
          print*,"read_equilibrium_sigma: U_eq in "//reg(sigfile)//" different from the input:",Ui
          stop
       endif
       write(*,"(3A)")"Header of the file:",reg(txtfy(u_))," ",reg(txtfy(dens))
       do i=1,Lf
          read(unit,*)wm,ims,res
          self%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
       selfHF(1) = dreal(self%iw(Lf))
       selfReg%iw = self%iw - selfHF(1)
       Scoeff  = tail_coeff_sigma(Ui,dens)
       allocate(Rtau(0:Lf))
       call fft_sigma_iw2tau(selfReg%iw,Rtau(0:),beta)!,Scoeff)
       call fft_extract_gtau(Rtau(0:),selfReg%mats(0:))
       deallocate(Rtau)
       !
    else
       !
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Sigma(iw)=Ui*(n-1/2)"
       dens=0.5d0
       selfHF(1) = Ui*(dens-0.5d0)    !set self-energy to zero or whatever is the initial HF solution
       selfReg   = zero
       !
    endif
    !
    selfReg%less(1,1) = -xi*selfReg%mats(L)                     !OK
    selfReg%ret(1,1) =   xi*(selfReg%mats(0)+selfReg%mats(L))   !OK
    forall(i=0:L)selfReg%lmix(1,i)=-xi*selfReg%mats(L-i)        !small errors near 0,beta
    !
    call add_kb_contour_gf(selfReg,selfHF,self,params)
    !
  end subroutine read_equilibrium_sigma_normal

  subroutine read_equilibrium_sigma_superc(selfHFB,selfReg,self,params)
    complex(8)              :: SelfHFB(:,:,:) ![2][2][Ntime]
    type(kb_contour_gf)     :: SelfReg(2,2)
    type(kb_contour_gf)     :: Self(2,2)
    type(kb_contour_params) :: params
    real(8)                 :: wm,res,ims,delta
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8)                 :: u_
    real(8),dimension(0:1)  :: Scoeff
    do i=1,2
       do j=1,2
          if(.not.Self(i,j)%status)  stop "read_equilibrium_sigma: Sigma is not allocated"
          if(.not.SelfReg(i,j)%status)  stop "read_equilibrium_sigma: Sigma_reg is not allocated"
       enddo
    enddo
    if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
    if(size(selfHFB,1)/=2)stop "read_equilibrium_sigma: size(selfHF,1)!=2"
    if(size(selfHFB,2)/=2)stop "read_equilibrium_sigma: size(selfHF,2)!=2"
    if(size(selfHFB,3)/=params%Ntime)stop "read_equilibrium_sigma: size(selfHF,3)!=Ntime"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    selfHFB=zero
    !
    !NORMAL COMPONENT
    inquire(file=trim(sigfile),exist=bool)
    if(bool)then
       !
       write(*,"(A)")"Reading initial Sigma(iw) from file "//reg(sigfile)
       unit = free_unit()
       i = file_length(trim(sigfile)) - 1
       open(unit,file=trim(sigfile),status='old')
       if(i/=Lf)then
          print*,"read_equilibrium_sigma: Liw in "//reg(sigfile)//" different from the input:",Lf
          stop
       endif
       do i=1,Lf
          read(unit,*)wm,ims,res
          self(1,1)%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
       selfHFB(1,1,1) = dreal(self(1,1)%iw(Lf))
       selfHFB(2,2,1) =-dreal(self(1,1)%iw(Lf)) !-conjg(selfHFB(1,1,1))
       selfReg(1,1)%iw = self(1,1)%iw - selfHFB(1,1,1)
       allocate(Rtau(0:Lf))
       call fft_sigma_iw2tau(selfReg(1,1)%iw,Rtau(0:),beta) !,C=[0d0,0d0])
       call fft_extract_gtau(Rtau(0:),selfReg(1,1)%mats(0:))
       deallocate(Rtau)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Sigma(iw)=Ui*(n-1/2)"
       dens=0.5d0
       selfHFB(1,1,1) = Ui*(dens-0.5d0)    !set self-energy to zero or whatever is the initial HF solution
       selfHFB(2,2,1) =-Ui*(dens-0.5d0)    !-conjg(selfHFB(1,1,1))
       selfReg(1,1)   = zero
    endif
    !ANOMAL COMPONENT
    inquire(file=trim(selfile),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial Self(iw) from file "//reg(selfile)
       unit = free_unit()
       i = file_length(trim(selfile))
       open(unit,file=trim(selfile),status='old')
       if(i/=Lf)then
          print*,"read_equilibrium_sigma: Liw in "//reg(selfile)//" different from the input:",Lf
          stop
       endif
       do i=1,Lf
          read(unit,*)wm,ims,res
          self(1,2)%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
       selfHFB(1,2,1) = dreal(self(1,2)%iw(Lf))
       selfHFB(2,1,1) = dreal(self(1,2)%iw(Lf))
       selfReg(1,2)%iw = self(1,2)%iw - selfHFB(1,2,1)
       allocate(Rtau(0:Lf))
       call fft_sigma_iw2tau(selfReg(1,2)%iw,Rtau(0:),beta,C=[0d0,0d0])
       call fft_extract_gtau(Rtau(0:),selfReg(1,2)%mats(0:))
       deallocate(Rtau)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Self(iw)=-Delta_SC"
       delta = deltasc
       selfHFB(1,2,1) = -delta    !set self-energy to zero or whatever is the initial HF solution
       selfHFB(2,1,1) = -delta
       selfReg(1,2)   = zero
    endif
    !
    !ACTHUNG BABY:
    selfReg(1,1)%less(1,1) = -xi*selfReg(1,1)%mats(L) !OK
    selfReg(1,1)%ret(1,1) =   xi*(selfReg(1,1)%mats(L)+selfReg(1,1)%mats(0))   !OK
    forall(i=0:L)selfReg(1,1)%lmix(1,i)=-xi*selfReg(1,1)%mats(L-i)     !small errors near 0,beta
    !
    selfReg(1,2)%less(1,1) = -xi*selfReg(1,2)%mats(L) !OK
    selfReg(1,2)%ret(1,1) =   xi*(selfReg(1,2)%mats(L)+self(1,2)%mats(0))   !OK
    forall(i=0:L)selfReg(1,2)%lmix(1,i)=-xi*selfReg(1,2)%mats(L-i)     !
    !
    call add_kb_contour_gf(selfReg(1,1),selfHFB(1,1,:),self(1,1),params)
    call add_kb_contour_gf(selfReg(1,2),selfHFB(1,2,:),self(1,2),params)
    !
    call get_bar(selfReg(2,2),selfReg(1,1),params)
    call get_bar(selfReg(2,1),selfReg(1,2),params)
    !
    call get_bar(self(2,2),self(1,1),params)
    call get_bar(self(2,1),self(1,2),params)
    !
    return
  end subroutine read_equilibrium_sigma_superc




  subroutine neq_fft_fg()
  end subroutine neq_fft_fg


  subroutine neq_fft_sigma()
  end subroutine neq_fft_sigma


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
