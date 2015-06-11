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
     module procedure neq_continue_equilibirum_bethe, &
                      neq_continue_equilibirum_, &
                      neq_continue_equilibirum__
  end interface neq_continue_equilibirum

  public  :: neq_continue_equilibirum     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.
  public  :: neq_setup_initial_conditions ! setup the initial conditions for the e/k dependent GF.


  ! public :: read_equilibrium_weiss
  ! public :: read_equilibrium_sigma
  ! public :: fft_extract_gtau

  real(8),public :: h1,h2,hDC,dens

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_(g0,gk,dgk,g,self,Hk,Wtk,params)
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
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    call read_equilibrium_weiss(g0,params,Hk=Hk,Wtk=Wtk)!<== get G0^{x=iw,tau,M,<,R,\lmix}
    ! Gcoeff = tail_coeff_glat(Ui,h1,h2,hDC)
    ! call fft_gf_iw2tau(g0%iw,g0%tau(0:),params%beta,Gcoeff)
    ! call extract_gtau_(g0%tau,g0%mats)
    ! g0%less(1,1) = -xi*g0%mats(L)
    ! g0%ret(1,1)  = -xi
    ! forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
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
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    G=zero
    do ik=1,Lk
       call neq_setup_initial_conditions(gk(ik),dgk(ik),self,hk(ik),params)
       G%mats(0:)  = G%mats(0:)  + wtk(ik)*gk(ik)%mats(0:)
       G%tau(0:)   = G%tau(0:)   + wtk(ik)*gk(ik)%tau(0:)
       G%iw(:)     = G%iw(:)     + wtk(ik)*gk(ik)%iw(:)
       G%ret(1,1)  = G%ret(1,1)  + wtk(ik)*gk(ik)%ret(1,1)
       G%less(1,1) = G%less(1,1) + wtk(ik)*gk(ik)%less(1,1)
       G%lmix(1,0:)= G%lmix(1,0:)+ wtk(ik)*gk(ik)%lmix(1,0:)
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
    !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
    call read_equilibrium_weiss(g0,params,wband=wband)!<== get G0^{x=iw,tau,M,<,R,\lmix}
    ! inquire(file=trim(g0file),exist=bool)
    ! if(bool)then
    !    write(*,"(A)")"Reading initial G0(iw) from file "//reg(g0file)
    !    unit = free_unit()
    !    open(unit,file=reg(g0file),status='old')
    !    i = file_length(reg(g0file))
    !    if(i/=Lf)then
    !       print*,"init_equilibrium_weiss_field: Lfreq in +g0file does not correspond",i
    !       stop
    !    endif
    !    do i=1,Lf
    !       read(unit,*)wm,ims,res
    !       g0%iw(i) = dcmplx(res,ims)
    !    enddo
    !    close(unit)
    ! else
    !    write(*,"(A)")"Start from Non-interacting G0(iw)"
    !    do i=1,Lf
    !       wm    = pi/beta*dble(2*i-1)
    !       zeta  = xi*wm
    !       g0%iw(i) = gfbethe(wm,zeta,wband)
    !    enddo
    ! endif
    ! !
    ! !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    ! Gcoeff = tail_coeff_glat(U,0.5d0,0d0,0d0)
    ! call fft_gf_iw2tau(g0%iw,g0%mats(0:),params%beta,Gcoeff)
    ! g0%less(1,1) = -xi*g0%mats(L)
    ! g0%ret(1,1)  = -xi
    ! forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call read_equilibrium_sigma(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
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
    real(8)                             :: nk
    complex(8)                          :: epsk
    complex(8),allocatable,dimension(:) :: SxG
    real(8),dimension(:),allocatable    :: ftau
    L = params%Ntau
    Lf= params%Niw
    do i=1,Lf
       Gk%iw(i) = one/(xi*params%wm(i) - Hk - Self%iw(i))
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
  end subroutine neq_setup_initial_conditions













  !+-------------------------------------------------------------------+
  !PURPOSE: Read equilibrium solution and initialize the corresponding
  ! function.
  ! - read_equilibrium_weiss: read and init Weiss FIeld G0
  ! - read_equilibrium_sigma: read and init Self-energy function Self
  !+-------------------------------------------------------------------+
  subroutine read_equilibrium_weiss(g0,params,wband,Hk,Wtk)
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
       write(*,"(A)")"Start from Non-interacting Bethe-model G0(iw)"
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
  end subroutine read_equilibrium_weiss

  subroutine read_equilibrium_sigma(self,params)
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
  end subroutine read_equilibrium_sigma





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



  subroutine neq_continue_equilibirum__(g0,gk,dgk,g,self,Hk,Wtk,params,Kerk)
    type(kb_contour_gf)                 :: g0(2,2)
    type(kb_contour_gf)                 :: gk(:,:),gk_aux(size(gk,1),2),Kerk(size(gk,1),4)
    type(kb_contour_dgf)                :: dgk(size(gk,1),2)
    type(kb_contour_gf)                 :: g(2,2)
    type(kb_contour_gf)                 :: self(2,2)
    type(kb_contour_params)             :: params
    real(8)                             :: Hk(size(gk)),Wtk(size(gk))
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: SxG
    Lk=size(gk,1)
!   if(.not.g0%status)stop "init_functions: g0 is not allocated"
!   if(.not.g%status)stop "init_functions: g is not allocated"
!   do ik=1,Lk
!      if(.not.gk(ik)%status)stop "init_functions: gk(ik) is not allocated"
!      if(.not.dgk(ik)%status)stop "init_functions: dgk(ik) is not allocated"
!   enddo
!   if(.not.self%status)stop "init_functions: self is not allocated"
!   if(.not.params%status)stop "init_functions: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
!   call read_equilibrium_weiss(g0,params,Hk=Hk,Wtk=Wtk)!<== get G0^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
!   call read_equilibrium_sigma(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
!   G=zero
    do ik=1,Lk
       call neq_setup_initial_conditions(gk_aux(ik,1),dgk(ik,1),self(1,1),hk(ik),params)
       call neq_setup_initial_conditions(gk_aux(ik,2),dgk(ik,2),self(1,1),-hk(ik),params)

       call convolute_kb_contour_gf(Kerk(ik,1),gk_aux(ik,1),self(1,2),params)
       call convolute_kb_contour_gf(Kerk(ik,2),Kerk(ik,1),gk_aux(ik,2),params)
       call convolute_kb_contour_gf(Kerk(ik,3),Kerk(ik,2),self(2,1),params)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Now the subroutine convolute, if N==0, does the convolution also in
       ! Matsubara frequencies (a mere multiplication) and in Matsubara
       ! imaginary times. In this way we can use the equation satisfied by G and
       ! F_bar also to initialise them.

       do i=1,Lf
          Gk(ik,1)%iw(i) = gk_aux(ik,1)%iw(i)/(1d0-Kerk(ik,3)%iw(i))
       enddo
       call fft_gf_iw2tau(Gk(ik,1)%iw,Gk(ik,1)%tau(0:),beta)
       call fft_extract_gtau(Gk(ik,1)%tau,Gk(ik,1)%mats)

       call convolute_kb_contour_gf(Kerk(ik,4),gk_aux(ik,2),self(2,1),params)
       call convolute_kb_contour_gf(Gk(ik,2),Kerk(ik,4),gk(ik,1),params)

       ! Alternatively, one can do the initialisation 'by hand'
       
       do i=i,Lf
          Gk(ik,1)%iw(i)=-(xi*params%wm(i)-hk(ik))/&
                          (params%wm(i)**2+(hk(ik)+self(1,1))**2+self(1,2)*self(2,1))
          Gk(ik,2)%iw(i)=-self(2,1)/(params%wm(i)**2+(hk(ik)+self(1,1))**2+self(1,2)*self(2,1))
       enddo
       call fft_gf_iw2tau(Gk(ik,1)%iw,Gk(ik,1)%tau(0:),beta)
       call fft_extract_gtau(Gk(ik,1)%tau,Gk(ik,1)%mats)
       call fft_gf_iw2tau(Gk(ik,2)%iw,Gk(ik,2)%tau(0:),beta)
       call fft_extract_gtau(Gk(ik,2)%tau,Gk(ik,2)%mats)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       G(1,1)%mats(0:)  = G(1,1)%mats(0:)  + wtk(ik)*gk(ik,1)%mats(0:)
       G(1,1)%tau(0:)   = G(1,1)%tau(0:)   + wtk(ik)*gk(ik,1)%tau(0:)
       G(1,1)%iw(:)     = G(1,1)%iw(:)     + wtk(ik)*gk(ik,1)%iw(:)
       G(1,1)%ret(1,1)  = G(1,1)%ret(1,1)  + wtk(ik)*gk(ik,1)%ret(1,1)
       G(1,1)%less(1,1) = G(1,1)%less(1,1) + wtk(ik)*gk(ik,1)%less(1,1)
       G(1,1)%lmix(1,0:)= G(1,1)%lmix(1,0:)+ wtk(ik)*gk(ik,1)%lmix(1,0:)

       G(2,1)%mats(0:)  = G(2,1)%mats(0:)  + wtk(ik)*gk(ik,2)%mats(0:)
       G(2,1)%tau(0:)   = G(2,1)%tau(0:)   + wtk(ik)*gk(ik,2)%tau(0:)
       G(2,1)%iw(:)     = G(2,1)%iw(:)     + wtk(ik)*gk(ik,2)%iw(:)
       G(2,1)%ret(1,1)  = G(2,1)%ret(1,1)  + wtk(ik)*gk(ik,2)%ret(1,1)
       G(2,1)%less(1,1) = G(2,1)%less(1,1) + wtk(ik)*gk(ik,2)%less(1,1)
       G(2,1)%lmix(1,0:)= G(2,1)%lmix(1,0:)+ wtk(ik)*gk(ik,2)%lmix(1,0:)

       call get_bar(G(2,2),G(1,1),params)
       call get_bar(G(1,2),G(2,1),params)
    enddo
    return
  end subroutine neq_continue_equilibirum__

end module NEQ_EQUILIBRIUM
