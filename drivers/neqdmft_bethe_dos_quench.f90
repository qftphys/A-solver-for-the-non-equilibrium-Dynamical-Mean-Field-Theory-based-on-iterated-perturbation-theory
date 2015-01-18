program neqDMFT
  !USE NEQ_DMFT_IPT
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_VARS_GLOBAL

  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                               :: i,j,ik,itime,iloop,ix,iy,Lk
  logical                               :: converged
  real(8)                               :: wband
  character(len=16)                     :: finput
  type(kb_contour_gf)                   :: Sbath
  type(kb_contour_gf)                   :: Gloc,Sigma,Gwf
  type(kb_contour_gf)                   :: Ker
  type(kb_contour_gf),allocatable       :: Gk(:)
  type(kb_contour_dgf),allocatable      :: dGk(:),dGk_old(:)
  type(kb_contour_dgf)                  :: Gedge,G0edge
  complex(8),dimension(:,:),allocatable :: Ham
  complex(8)                            :: sigma_gtr
  !RESULTS:
  real(8),dimension(:),allocatable      :: epsik,wt
  real(8),dimension(:,:),allocatable    :: nk
  real(8)                               :: D,de,intwt

  !READ THE INPUT FILE (in vars_global):
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.conf')
  call parse_input_variable(wband,"wband",finput,default=2d0,comment="half-bandwidth W=2t")
  call parse_input_variable(Lk,"Lk",finput,default=100,comment="Number of energy levels for Bethe DOS")
  call read_input_init(trim(finput))


  !BUILD THE LATTICE STRUCTURE (in lib/square_lattice):
  allocate(epsik(Lk),wt(Lk))
  epsik = linspace(-wband,wband,Lk,mesh=de)
  do i=1,Lk
     wt(i) = dens_bethe(epsik(i),wband)
  enddo
  call splot("DOSbethe.ipt",epsik,wt)
  wt=wt/sum(wt)


  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  !=====================================================================
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,Lfreq)
  call setup_kb_contour_params(cc_params,dt,beta)
  ! !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  ! call allocate_kb_contour_params(cc_params,Ntime,Ntau,dt,beta,Lfreq)
  ! cc_params%t = linspace(0.d0,cc_params%tmax,cc_params%Ntime,mesh=dt)
  ! cc_params%tau(0:) = linspace(0.d0,cc_params%beta,cc_params%Ntau+1,mesh=dtau)
  ! cc_params%wm  = pi/cc_params%beta*dble(2*arange(1,cc_params%Lf)-1)
  ! print*,"dt=",cc_params%dt
  ! print*,"dtau=",cc_params%dtau


  !ALLOCATE ALL THE FUNCTIONS INVOLVED IN THE CALCULATION:
  call allocate_kb_contour_gf(Sigma,cc_params) !Self-Energy function
  call allocate_kb_contour_gf(Gloc,cc_params)  !Local Green's function
  call allocate_kb_contour_gf(Gwf,cc_params)   !Local Weiss-Field function
  call allocate_kb_contour_gf(Ker,cc_params)
  call allocate_kb_contour_dgf(Gedge,cc_params)
  allocate(Gk(Lk),dGk(Lk),dGk_old(Lk))
  do ik=1,Lk
     call allocate_kb_contour_gf(Gk(ik),cc_params)
     call allocate_kb_contour_dgf(dGk(ik),cc_params)
     call allocate_kb_contour_dgf(dGk_old(ik),cc_params)
     Gk(ik)=zero
     dGk(ik)=zero
     dGk_old(ik)=zero
  end do
  allocate(ham(cc_params%Ntime,Lk))
  allocate(nk(cc_params%Ntime,Lk))



  !READ OR GUESS THE INITIAL WEISS FIELD G0 (t=t'=0)
  cc_params%Nt=1
  Gloc = zero
  call init_equilibrium_functions(Gwf,Gk,dGk,Gloc,Sigma,cc_params)
  call measure_observables(Gloc,Sigma,cc_params)
  do ik=1,Lk
     nk(1,ik)=dimag(Gk(ik)%less(1,1))
  enddo
  do i=1,cc_params%Ntime
     do ik=1,Lk
        ham(i,ik)=epsik(ik)
     enddo
  enddo



  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Nt=itime
     !prepare the weiss-field at this actual time_step for DMFT:
     call setup_weiss_field(Gwf,cc_params)
     do ik=1,Lk
        dGk_old(ik) = dGk(ik)
     enddo

     !DMFT LOOP:
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        write(*,"(A,I4,A1)",advance='no')"dmft loop=",iloop," "

        !IMPURITY SOLVER: IPT.
        !GET SIGMA FROM THE WEISS FIELD
        call neq_solve_ipt(Gwf,Sigma,cc_params)

        !PERFORM THE SELF_CONSISTENCY: get local GF + update Weiss-Field
        Gedge%ret(1:itime)=zero
        Gedge%less(1:itime)=zero
        Gedge%lmix(0:)=zero
        do ik=1,Lk
           call vide_kb_contour_gf(Ham(:,ik),Sigma,Gk(ik),dGk_old(ik),dGk(ik),cc_params)
           Gedge%ret(1:itime) = Gedge%ret(1:itime)  + wt(ik)*Gk(ik)%ret(itime,1:itime)
           Gedge%less(1:itime)= Gedge%less(1:itime) + wt(ik)*Gk(ik)%less(itime,1:itime)
           Gedge%lmix(0:)     = Gedge%lmix(0:)      + wt(ik)*Gk(ik)%lmix(itime,0:)
        enddo
        Gloc%ret(itime,1:itime)   = Gedge%ret(1:itime)
        Gloc%less(itime,1:itime)  = Gedge%less(1:itime)
        Gloc%lmix(itime,0:)       = Gedge%lmix(0:)
        Gloc%less(1:itime-1,itime)=-conjg(Gedge%less(1:itime-1))


        !update the weiss field by solving the integral equation:
        ! G0 + K*G0 = Q , with K = G*\Sigma and Q = G
        call convolute_kb_contour_gf(Gloc,Sigma,Ker,cc_params,dcoeff=-1.d0)
        call vie_kb_contour_gf(Gloc,Ker,Gwf,cc_params)
        !
        !CHECK CONVERGENCE
        converged = convergence_check(Gwf,cc_params)
     enddo



     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
     call measure_observables(Gloc,Sigma,cc_params)
     forall(ik=1:Lk)nk(itime,ik)=dimag(Gk(ik)%less(itime,itime))
  enddo


  !EVALUATE AND PRINT OTHER RESULTS OF THE CALCULATION
  call splot3d("nkVSepsikVStime.ipt",cc_params%t,epsik,nk)
  call plot_kb_contour_gf("Sigma.ipt",Sigma,cc_params)
  call plot_kb_contour_gf("Gloc.ipt",Gloc,cc_params)
  call plot_kb_contour_gf("G0.ipt",Gwf,cc_params)
  print*,"BRAVO"


contains




  subroutine init_equilibrium_functions(g0,gk,dgk,g,self,params)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_gf)                 :: gk(:)
    type(kb_contour_dgf)                :: dgk(size(gk))
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: SxG
    real(8)                             :: Scoeff(2),Gcoeff(4)
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
          zeta  = dcmplx(0.d0,wm)
          g0%iw(i) = sum_overk_zeta(zeta,epsik,wt)
       enddo
    endif
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    Gcoeff = tail_coeff_glat(U,0.5d0,0d0,0d0)
    call fft_gf_iw2tau(g0%iw,g0%mats(0:),params%beta,Gcoeff)
    ! call fftgf_iw2tau(g0%iw,g0%mats(0:),params%beta)
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
    ! call fftgf_tau2iw(Self%mats(0:),Self%iw,beta)
    Self%iw = xi*dimag(self%iw) !!ACTHUNG: imposing half-filling symmetry
    Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
    Self_gtr      =(xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
    do j=0,L
       Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
    end do
    Self%ret(1,1) = Self_gtr - Self%less(1,1)
    !   
    do ik=1,Lk 
       call setup_initial_conditions(self,gk(ik),dgk(ik),epsik(ik),params)
       G%mats(0:)  = G%mats(0:)  + wt(ik)*gk(ik)%mats(0:)
       G%iw(:)     = G%iw(:)     + wt(ik)*gk(ik)%iw(:)
       G%ret(1,1)  = G%ret(1,1)  + wt(ik)*gk(ik)%ret(1,1)
       G%less(1,1) = G%less(1,1) + wt(ik)*gk(ik)%less(1,1)
       G%lmix(1,0:)= G%lmix(1,0:)+ wt(ik)*gk(ik)%lmix(1,0:)
    enddo
    return
  end subroutine init_equilibrium_functions



  subroutine setup_initial_conditions(Self,Gk,dGk,hk,params)
    type(kb_contour_gf)                 :: Gk,Self
    type(kb_contour_dgf)                :: dGk
    real(8)                             :: hk
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,Ltau
    real(8)                             :: nktmp
    complex(8),allocatable,dimension(:) :: SxG
    Ltau  = params%Ntau
    Gk%iw = one/(xi*params%wm - hk - self%iw)          !get G_k(iw) 
    call fft_gf_iw2tau(Gk%iw,Gk%mats(0:),beta)        !get G_k(tau)
    nktmp = -Gk%mats(Ltau)                           !n(k,t=0)=-G^M_k(beta)=G^M_k(0-)
    Gk%less(1,1) =  xi*nktmp                         !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk%ret(1,1)  = -xi                               !get G^R_k(0,0)=-xi
    forall(i=0:Ltau)Gk%lmix(1,i)=-xi*Gk%mats(Ltau-i) !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !
    !Derivatives
    allocate(SxG(0:Ltau))
    !get d/dt G_k^R = -i*e(k,0)G_k^R
    dGk%ret(1)  = -xi*hk*Gk%ret(1,1)            
    !get d/dt G_k^< = -i*e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    do k=0,Ltau
       SxG(k)=Self%lmix(1,k)*conjg(Gk%lmix(1,Ltau-k))
    end do
    dGk%less(1) = -xi*hk*Gk%less(1,1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,Ltau) 
    !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(0:)= -xi*hk*Gk%lmix(1,0:)
    do j=0,Ltau
       do k=0,j
          SxG(k)=Self%lmix(1,k)*Gk%mats(Ltau+k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       do k=j,Ltau
          SxG(k)=Self%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,Ltau) 
    enddo
  end subroutine setup_initial_conditions






  subroutine setup_weiss_field(g0,params)
    type(kb_contour_gf)                   :: g0
    type(kb_contour_params)               :: params
    integer                               :: i,j,k,N,L
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
  end subroutine setup_weiss_field



  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve with the 2^nd IPT sigma functions
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt(G0,Sigma,params)
    type(kb_contour_gf)                   :: G0
    type(kb_contour_gf)                   :: Sigma
    type(kb_contour_params)               :: params
    integer                               :: N,L
    complex(8),dimension(:,:),allocatable :: G0_gtr,Sigma_gtr,G0_rmix
    integer                               :: i,j,itau
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    allocate(G0_gtr(N,N),Sigma_gtr(N,N))
    G0_gtr(N,1:N)=G0%less(N,1:N)+ G0%ret(N,1:N)
    G0_gtr(1:N-1,N)=G0%less(1:N-1,n)-conjg(G0%ret(N,1:N-1))
    !
    !Vertical edge
    do j=1,N
       Sigma%less(N,j)= U*U*G0%less(N,j)*G0_gtr(j,N)*G0%less(N,j)
       Sigma_gtr(N,j) = U*U*G0_gtr(N,j)*G0%less(j,N)*G0_gtr(N,j)
    end do
    !Horizontal edge
    do i=1,N-1
       Sigma%less(i,N)= U*U*G0%less(i,N)*G0_gtr(N,i)*G0%less(i,N)
       Sigma_gtr(i,N) = U*U*G0_gtr(i,N)*G0%less(N,i)*G0_gtr(i,N)
    end do
    !Imaginary time edge:
    forall(i=0:L)Sigma%lmix(N,i)  = U*Ui*G0%lmix(N,i)*(conjg(G0%lmix(N,L-i)))*G0%lmix(N,i)
    forall(j=1:N)Sigma%ret(N,j) = Sigma_gtr(N,j) - Sigma%less(N,j)

  end subroutine neq_solve_ipt



  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check(G,params) result(converged)
    type(kb_contour_gf)                 :: G
    type(kb_contour_params)             :: params
    logical                             :: converged
    complex(8),dimension(:),allocatable :: test_func
    integer                             :: i,N,L,Ntot
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    Ntot=2*N+L+1
    allocate(test_func(Ntot))
    test_func=zero
    do i=1,N
       test_func(i)  = G%ret(N,i)
       test_func(N+i)= G%less(N,i)
    enddo
    do i=0,L
       test_func(2*N+i+1)=G%lmix(N,i)
    enddo
    converged=check_convergence(test_func,dmft_error,Nsuccess,nloop)
    deallocate(test_func)
    !if(isnan(err))stop "Aborted convergence: error=NaN"
  end function convergence_check





  !+-------------------------------------------------------------------+
  !PURPOSE: measure some observables and print them
  !+-------------------------------------------------------------------+
  subroutine measure_observables(g,self,params)
    type(kb_contour_gf)     :: g
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    integer                 :: unit,itime
    real(8)                 :: dens,docc,ekin,epot,etot
    itime = params%Nt
    unit = free_unit()
    open(unit,file="observables.plot",position="append")
    dens = measure_dens(g,self,params)
    docc = measure_docc(g,self,params)
    ekin = measure_ekin(g,self,params)
    epot = measure_epot(g,self,params)
    etot = ekin + epot
    write(unit,"(6F20.12)")params%t(itime),dens,docc,ekin,epot,etot
    close(unit)
  end subroutine measure_observables



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the density at a given istant of time
  ! n(t)=-xi*G^<(t,t)
  !+-------------------------------------------------------------------+
  function measure_dens(g,self,params) result(dens)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: dens
    integer                             :: N
    N = params%Nt
    dens = dimag(G%less(N,N))
  end function measure_dens


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the double occupancy at a given istant of time
  ! d(t)=n_up(t)*n_do(t)-1/U0*[Self^M*G^M]
  !      n_up(t)*n_do(t)-i/U*[Self^R*G^< + Self^<*G^A + Self^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_docc(g,self,params) result(docc)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: docc
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: SxG
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    nt   = dimag(G%less(N,N))
    allocate(SxG(0:max(N,L)))
    docc = nt**2
    if(N==1)then
       if(ui/=0.d0)then
          do k=0,L
             SxG(k)=Self%mats(L-k)*G%mats(k)
          end do
          docc=docc-1.d0/Ui*params%dtau*kb_trapz(SxG(0:),0,L)
       endif
    else
       if(u/=0.d0)then
          do k=0,L
             SxG(k)=Self%lmix(N,k)*conjg(G%lmix(N,L-k))
          end do
          docc=docc + 1.d0/U*params%dtau*dimag( (-xi)*kb_trapz(SxG(0:),0,L) )
          do k=1,N
             SxG(k)=Self%ret(N,k)*G%less(k,N)
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
          do k=1,N
             SxG(k)=Self%less(N,k)*conjg(G%ret(N,k))
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
       endif
    endif
    deallocate(SxG)
  end function measure_docc


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_ekin(g,self,params) result(ekin)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: ekin
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: Ker
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    allocate(Ker(0:max(N,L)))
    if(N==1)then
       do k=0,L
          Ker(k)=G%mats(L-k)*G%mats(k)
       end do
       ekin = -2.d0*params%dtau*kb_trapz(Ker(0:),0,L)
    else
       do k=0,L
          Ker(k)=G%lmix(N,k)*conjg(G%lmix(N,L-k))
       end do
       ekin=2.d0*params%dtau*dimag( (-xi)*kb_trapz(Ker(0:),0,L) )
       do k=1,N
          Ker(k)=G%ret(N,k)*G%less(k,N)
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
       do k=1,N
          Ker(k)=G%less(N,k)*conjg(G%ret(N,k))
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
    endif
    deallocate(Ker)
  end function measure_ekin



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! U(t)= U*docc(t) - n(t) + 1/4
  !+-------------------------------------------------------------------+
  function measure_epot(g,self,params) result(epot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: epot,docc,nt
    integer                             :: i,k,j,N,L
    N = params%Nt
    L = params%Ntau
    !
    if(N==1)then
       nt   = measure_dens(g,self,params)
       docc = measure_docc(g,self,params)
       epot = Ui*(docc - nt + 0.25d0)
    else
       nt   = measure_dens(g,self,params)
       docc = measure_docc(g,self,params)
       epot = U*(docc - nt + 0.25d0)
    endif
  end function measure_epot



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_etot(g,self,params) result(etot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: etot,ekin,epot
    ekin = measure_ekin(g,self,params)
    epot = measure_epot(g,self,params)
    etot = ekin + epot
  end function measure_etot



end PROGRAM neqDMFT
