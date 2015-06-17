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
     module procedure neq_continue_equilibirum_normal
     module procedure neq_continue_equilibirum_bethe
  end interface neq_continue_equilibirum

  public  :: neq_continue_equilibirum     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.
  public  :: neq_setup_initial_conditions ! setup the initial conditions for the e/k dependent GF.


  real(8)                          :: h1,h2,hDC,dens
  real(8),allocatable,dimension(:) :: Rtau

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_normal(g0,gk,dgk,g,selfHF,selfReg,self,Hk,Wtk,params)
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
    call read_equilibrium_sigma(selfHF,selfReg,self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
    !
    !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
    G=zero
    do ik=1,Lk
       call neq_setup_initial_conditions(gk(ik),dgk(ik),self,hk(ik),params)
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
  end subroutine neq_continue_equilibirum_normal


  subroutine neq_continue_equilibirum_bethe(g0,dg0,g,selfHF,selfReg,self,params,wband)
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
    call read_equilibrium_sigma(selfHF,selfReg,self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
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
    !get d/dt G_k^R = -i e(k,0)G_k^R
    dGk%ret(1)  = -xi*Hk*Gk%ret(1,1)
    !get d/dt G_k^< = -i e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    do k=0,L
       SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k))
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
  ! We assume that the read self-energy function INCLUDES the Hartree-Fock
  ! term, so we are able to extract it and separate regular from local part.
  ! - read_equilibrium_sigma: read and init Self-energy function Self
  ! - read_equilibrium_weiss: read and init Weiss FIeld G0
  !+-------------------------------------------------------------------+
  subroutine read_equilibrium_sigma(selfHF,selfReg,self,params)
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
  end subroutine read_equilibrium_sigma




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






  ! subroutine read_equilibrium_weiss(g0,params,wband,Hk,Wtk)
  !   type(kb_contour_gf)     :: g0
  !   type(kb_contour_params) :: params
  !   real(8),optional        :: Hk(:)
  !   real(8),optional        :: Wtk(:)
  !   real(8),optional        :: wband
  !   real(8)                 :: wband_
  !   real(8)                 :: wm,res,ims,mu
  !   logical                 :: bool
  !   complex(8)              :: zeta
  !   integer                 :: i,j,k,ik,unit,len,N,L,Lf
  !   real(8)                 :: Gcoeff(4)
  !   wband_=1d0;if(present(wband))wband_=wband
  !   !
  !   if(.not.g0%status)stop "read_equilibrium_g0: g0 is not allocated"
  !   if(.not.params%status)stop "read_equilibrium_g0: params is not allocated"
  !   N = params%Nt
  !   L = params%Ntau
  !   Lf= params%Niw
  !   !
  !   !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
  !   inquire(file=trim(g0file),exist=bool)
  !   if(bool)then
  !      write(*,"(A)")"Reading initial G0(iw) from file"//trim(g0file)
  !      unit = free_unit()
  !      i = file_length(trim(g0file))-1
  !      open(unit,file=trim(g0file),status='old')
  !      if(i/=Lf)then
  !         print*,"read_equilibrium_weiss: Liw in "//reg(g0file)//" different from the input:",Lf
  !         print*,"read_equilibrium_weiss: check the header of the file xmu, <h>, <h**2>, <h_dc>"
  !         stop
  !      endif
  !      read(unit,*)mu,h1,h2,hDC
  !      if(mu/=xmu)then
  !         print*,"read_equilibrium_weiss: mu in "//reg(g0file)//" different from the input:",xmu
  !         stop
  !      endif
  !      write(*,"(A)")"Header of the file:",reg(txtfy(mu)),reg(txtfy(h1)),reg(txtfy(h2)),reg(txtfy(hDC))
  !      do i=1,Lf
  !         read(unit,*)wm,ims,res
  !         g0%iw(i) = dcmplx(res,ims)
  !      enddo
  !      close(unit)
  !   else
  !      write(*,"(A)")"Start from Non-interacting Bethe-model G0(iw)"
  !      mu=xmu ; h1=0d0 ; h2=0d0 ; hDC=0d0
  !      write(*,"(5A)")"Header of the file:",reg(txtfy(mu)),reg(txtfy(h1)),reg(txtfy(h2)),reg(txtfy(hDC))
  !      if(present(Hk))then
  !         if(.not.present(Wtk))stop "read_equilibrium_weiss: Wtk not present"
  !         do i=1,Lf
  !            wm    = pi/beta*(2*i-1)
  !            zeta  = xi*wm + xmu             
  !            g0%iw(i) = sum_overk_zeta(zeta,Hk,Wtk)
  !         enddo
  !      else
  !         do i=1,Lf
  !            wm    = pi/beta*(2*i-1)
  !            zeta  = xi*wm + xmu
  !            g0%iw(i) = gfbethe(wm,zeta,wband_)
  !         enddo
  !      endif
  !   endif
  !   !
  !   Gcoeff = tail_coeff_g0(Ui,h1,h2,hDC)
  !   call fft_gf_iw2tau(g0%iw,g0%tau(0:),params%beta,Gcoeff)
  !   call fft_extract_gtau(g0%tau,g0%mats)
  !   g0%less(1,1) = -xi*g0%mats(L)
  !   g0%ret(1,1)  = -xi
  !   forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
  ! end subroutine read_equilibrium_weiss


end module NEQ_EQUILIBRIUM
