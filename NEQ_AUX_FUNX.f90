module NEQ_AUX_FUNX
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_INPUT_VARS
  USE NEQ_VARS_GLOBAL
  USE NEQ_IPT
  USE CONSTANTS
  USE FUNCTIONS
  USE IOTOOLS
  USE DMFT_TOOLS
  implicit none

  private

  interface neq_continue_equilibirum
     module procedure neq_continue_equilibirum_,neq_continue_equilibirum_bethe
  end interface neq_continue_equilibirum

  public :: neq_continue_equilibirum     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.
  public :: neq_setup_weiss_field        ! continue the weiss field to the next time step.
  public :: neq_setup_initial_conditions ! setup the initial conditions for the e/k dependent GF.
contains



  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  ! - neq_continue_equilibirum_: generic one (1band)
  ! - neq_continue_equilibirum_bethe: this is for the Bethe lattice
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_(g0,gk,dgk,g,self,Hk,Wtk,params)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_gf)                 :: gk(:)
    type(kb_contour_dgf)                :: dgk(size(gk))
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: Hk(size(gk)),Wtk(size(gk))
    real(8)                             :: hamk
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
    N = params%Itime
    L = params%Ntau

    Lf= params%Niw
    !
    call read_equilibrium_g0(g0,params,Hk=Hk,wtk=Wtk)
    !
    call neq_solve_ipt_first_step(g0,self,params)
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    do ik=1,Lk
       call neq_setup_initial_conditions(gk(ik),dgk(ik),self,Hk(ik),params)
       G%mats(:)   = G%mats(:)   + wtk(ik)*gk(ik)%mats(:)
       G%iw(:)     = G%iw(:)     + wtk(ik)*gk(ik)%iw(:)
       G%ret(1,1)  = G%ret(1,1)  + wtk(ik)*gk(ik)%ret(1,1)
       G%less(1,1) = G%less(1,1) + wtk(ik)*gk(ik)%less(1,1)
       G%lmix(1,:) = G%lmix(1,:) + wtk(ik)*gk(ik)%lmix(1,:)
    enddo
    return
  end subroutine neq_continue_equilibirum_

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
    real(8),dimension(0:3)              :: Gcoeff
    wband_=1d0;if(present(wband))wband_=wband
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.dg0%status)stop "init_functions: dg0 is not allocated"
    if(.not.g%status)stop "init_functions: g is not allocated"
    if(.not.self%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    N = params%Itime
    L = params%Ntau
    Lf= params%Niw
    !
    call read_equilibrium_g0(g0,params,wband=wband_)
    !
    call neq_solve_ipt_first_step(g0,self,params)
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    do i=1,Lf
       wm      = pi/beta*dble(2*i-1)
       zeta    = xi*wm - Self%iw(i)
       G%iw(i) = gfbethe(wm,zeta,wband_)
    enddo
    Gcoeff      = tail_coeff_glat(uloc,0.5d0,0d0,0d0)
    call fft_gf_iw2tau(G%iw,G%mats,beta,Gcoeff)     !get G(tau)
    G%ret(1,1)  = -xi                               !get G^R(0,0)=-xi
    G%less(1,1) = -xi*G%mats(L)                     !get G^<(0,0)= xi*G^M(0-)
    forall(i=1:L)G%lmix(1,i)=-xi*G%mats(L-i+1)      !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
    !Derivatives
    allocate(GxG0(L))
    !get d/dt G0^R = 0.d0
    dG0%ret(1)  = zero
    !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
    do k=1,L
       GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k+1))
    end do
    dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0,1,L) 
    !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
    dG0%lmix= zero
    do j=1,L
       do k=1,j
          GxG0(k)=G%lmix(1,k)*G0%mats(k + L-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0,1,j)
       do k=j,L
          GxG0(k)=G%lmix(1,k)*G0%mats(k-j+1)
       end do
       dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0,j,L) 
    enddo
    deallocate(GxG0)
    return
  end subroutine neq_continue_equilibirum_bethe







  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for e/k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine neq_setup_initial_conditions(Gk,dGk,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk,Self
    type(kb_contour_dgf)                :: dGk
    real(8)                             :: Hk
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    real(8)                             :: nk
    complex(8)                          :: epsk
    complex(8),allocatable,dimension(:) :: SxG
    real(8),dimension(0:3)              :: Gcoeff
    L = params%Ntau
    Gk%iw   = one/(xi*params%wm - Hk - Self%iw)
    Gcoeff  = tail_coeff_gk(uloc,0.5d0,0d0,Hk)
    call fft_gf_iw2tau(Gk%iw,Gk%mats,beta,Gcoeff)        !get G_k(tau)
    !                                                    !n(k,t=0)=-G^M_k(beta)=G^M_k(0-)
    Gk%less(1,1) = -xi*Gk%mats(L)                        !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk%ret(1,1)  = -xi                                   !get G^R_k(0,0)=-xi
    forall(i=1:L)Gk%lmix(1,i)=-xi*Gk%mats(L-i+1)         !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !Derivatives
    allocate(SxG(L))
    !get d/dt G_k^R = -i e(k,0)G_k^R
    dGk%ret(1)  = -xi*Hk*Gk%ret(1,1)
    !get d/dt G_k^< = -i e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    do k=1,L
       SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k+1))
    end do
    dGk%less(1) = -xi*Hk*Gk%less(1,1)-&
         xi*(-xi)*params%dtau*kb_trapz(SxG,1,L) 
    !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(:) = -xi*Hk*Gk%lmix(1,:)
    do j=1,L
       do k=1,j
          SxG(k)=self%lmix(1,k)*Gk%mats(k + L-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG,1,j)
       do k=j,L
          SxG(k)=self%lmix(1,k)*Gk%mats(k - j+1)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG,j,L) 
    enddo
  end subroutine neq_setup_initial_conditions










  

  !+-------------------------------------------------------------------+
  !PURPOSE: setup the Weiss Field G0 for the next time-step
  !+-------------------------------------------------------------------+
  subroutine neq_setup_weiss_field(g0,params)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,N,L
    complex(8),allocatable,dimension(:) :: SxG
    if(.not.g0%status)stop "init_g0: g0 is not allocated"
    if(.not.params%status)stop "init_g0: params is not allocated"
    !
    N = params%Itime
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
       do j=1,L
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
       do k=1,L
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











  !+-------------------------------------------------------------------+
  !PURPOSE  : Read equilibrium solution and initialize the corresponding
  ! function:
  ! - read_equilibrium_g0   : read and init Weiss FIeld G0
  ! - read_equilibrium_sigma: read and init Self-energy function Self
  !+-------------------------------------------------------------------+
  subroutine read_equilibrium_g0(g0,params,wband,Hk,Wtk)
    type(kb_contour_gf)     :: g0
    type(kb_contour_params) :: params
    real(8),optional        :: Hk(:)
    real(8),optional        :: Wtk(:)
    real(8),optional        :: wband
    real(8)                 :: wband_
    real(8)                 :: wm,res,ims
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8)                 :: mu,h1,h2,hDC
    real(8),dimension(0:3)  :: Gcoeff
    wband_=1d0;if(present(wband))wband_=wband
    if(.not.g0%status)stop "read_equilibrium_g0: g0 is not allocated"
    if(.not.params%status)stop "read_equilibrium_g0: params is not allocated"
    N = params%Itime
    L = params%Ntau
    Lf= params%Niw
    !
    !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
    inquire(file=trim(g0file),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial G0(iw) from file"//trim(g0file)
       unit = free_unit()
       open(unit,file=trim(g0file),status='old')
       i = file_length(trim(g0file))-1
       if(i/=Lf)then
          print*,"read_equilibrium_g0: Liw in "//reg(g0file)//"different from the input:",Lf
          stop
       endif
       read(unit,*)mu,h1,h2,hDC
       if(mu/=xmu)then
          print*,"read_equilibrium_g0: mu in "//reg(g0file)//"different from the input:",xmu
          stop
       endif
       write(*,"(A)")"Header of the file:",reg(txtfy(mu)),reg(txtfy(h1)),reg(txtfy(h2)),reg(txtfy(hDC))
       do i=1,Lf
          read(unit,*)wm,ims,res
          g0%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       write(*,"(A)")"Start from Non-interacting Bethe-model G0(iw)"
       mu=xmu ; h1=0d0 ; h2=0d0 ; hDC=0d0
       write(*,"(A)")"Header of the file:",reg(txtfy(mu)),reg(txtfy(h1)),reg(txtfy(h2)),reg(txtfy(hDC))
       if(present(Hk))then
          if(.not.present(Wtk))stop "read_equilibrium_g0: Wtk not present"
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu             
             g0%iw(i) = sum_overk_zeta(zeta,Hk,Wtk)
          enddo
       else
          do i=1,Lf
             wm    = pi/beta*(2*i-1)
             zeta  = xi*wm + xmu
             g0%iw(i) = gfbethe(wm,zeta,wband_)
          enddo
       endif
    endif
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    Gcoeff  = tail_coeff_g0(mu,h1,h2,hDC)
    call fft_gf_iw2tau(g0%iw,g0%mats,params%beta,Gcoeff)
    g0%less(1,1) = -xi*g0%mats(L)
    g0%ret(1,1)  = -xi
    forall(i=1:L)g0%lmix(1,i)=-xi*g0%mats(L-i+1)
  end subroutine read_equilibrium_g0



  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading self-energy
  !+-------------------------------------------------------------------+
  subroutine read_equilibrium_sigma(self,params)
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    real(8)                 :: wm,res,ims
    logical                 :: bool
    complex(8)              :: zeta
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8)                 :: u_,dens_
    logical                 :: hfmode_
    real(8),dimension(0:1)  :: Scoeff
    if(.not.self%status)stop "read_equilibrium_sigma: sigma is not allocated"
    if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
    N = params%Itime
    L = params%Ntau
    Lf= params%Niw
    !CHECK IF SIGMA(IW) IS AVAILABLE (START FROM THE EQUILIBRIUM SOLUTION)
    !IF NOT, START FROM NON-INTERACTING (SIGMA=0)
    inquire(file=trim(sigfile),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial Sigma(iw) from file "//reg(sigfile)
       unit = free_unit()
       open(unit,file=trim(sigfile),status='old')
       i = file_length(trim(sigfile)) - 1
       if(i/=Lf)then
          print*,"read_equilibrium_sigma: Liw in "//reg(sigfile)//"different from the input:",Lf
          stop
       endif
       read(unit,*)u_,dens_,hfmode_
       if(u_/=U0)then
          print*,"read_equilibrium_sigma: U_eq in "//reg(sigfile)//"different from the input:",u0
          stop
       endif
       write(*,"(A)")"Header of the file:",reg(txtfy(u_)),reg(txtfy(dens_)),reg(txtfy(hfmode_))
       do i=1,Lf
          read(unit,*)wm,res,ims
          self%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
       !NOW YOU NEED TO PROPAGATE THE SOLUTION to Sigma^{x=M,<,R,\lmix}
       Scoeff  = tail_coeff_sigma(u_,dens_,hfmode_)
       call fft_sigma_tau2iw(Self%iw,Self%mats,beta,Scoeff)
       self%less(1,1) = -xi*self%mats(L)
       self%ret(1,1) =  xi*(self%mats(1)+self%mats(L))
       forall(i=1:L)self%lmix(1,i)=-xi*self%mats(L-i+1)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Sigma(iw)"
       self=zero                !set self-energy to zero or whatever is the initial HF solution
    endif
  end subroutine read_equilibrium_sigma



end module NEQ_AUX_FUNX
