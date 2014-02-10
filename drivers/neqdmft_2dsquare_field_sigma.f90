!###################################################################
!PURPOSE  : Solve conduction band electrons driven 
! by electric field interacting with a resevoir of free 
! electrons at temperature T
!AUTHORS  : Adriano Amaricci 
!###################################################################
program neqDMFT
  USE SQUARE_LATTICE
  !
  USE NEQ_VARS_GLOBAL   !global variables, calls to 3rd library 
  USE ELECTRIC_FIELD    !contains electric field init && routines
  USE NEQ_THERMOSTAT    !contains bath inizialization
  USE NEQ_IPT           !performs the non-eq. IPT. Write Sigma
  !USE NEQ_RESULTS
  implicit none
  integer           :: i,j,ik,itime,iloop,ix,iy
  logical           :: converged
  character(len=16) :: finput
  type(kb_contour_gf)                 :: Sbath
  type(kb_contour_gf)                 :: Gloc,Sigma,Gwf
  type(kb_contour_gf)                 :: Ker
  type(kb_contour_gf),allocatable     :: Gk(:)
  type(kb_contour_dgf),allocatable    :: dGk(:),dGk_old(:)
  type(kb_contour_dgf)                :: Gedge
  complex(8),dimension(:),allocatable :: Ham
  complex(8)                          :: sigma_gtr
  !RESULTS:
  real(8),dimension(:,:,:),allocatable  :: nDens
  real(8),dimension(:),allocatable      :: n_t
  real(8),dimension(:,:),allocatable    :: nk
  type(vect2D),dimension(:),allocatable :: Jloc
  type(vect2D)                          :: kt,Jk

  !READ THE INPUT FILE (in vars_global):
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.in')
  call read_input_init(trim(finput))


  !BUILD THE LATTICE STRUCTURE (in lib/square_lattice):
  Lk   = square_lattice_dimension(Nx,Nx)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)
  call get_free_dos(epsik,wt)


  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,dt,beta,Lfreq)
  cc_params%t = linspace(0.d0,cc_params%tmax,cc_params%Ntime,mesh=dt)
  cc_params%tau(0:) = linspace(0.d0,cc_params%beta,cc_params%Ntau+1,mesh=dtau)
  cc_params%wm  = pi/cc_params%beta*dble(2*arange(1,cc_params%Lf)-1)
  print*,"dt  =",cc_params%dt
  print*,"dtau=",cc_params%dtau


  !SET THE ELECTRIC FIELD (in electric_field):
  call set_efield_vector()
  call print_field(cc_params%t)


  !SET THE THERMOSTAT FUNCTION (in neq_thermostat):
  call allocate_kb_contour_gf(Sbath,cc_params)
  call get_thermostat_bath(Sbath,cc_params)


  !ALLOCATE ALL THE FUNCTIONS INVOLVED IN THE CALCULATION:
  call allocate_kb_contour_gf(Sigma,cc_params) !Self-Energy function
  call allocate_kb_contour_gf(Gloc,cc_params)  !Local Green's function
  call allocate_kb_contour_gf(Gwf,cc_params)   !Local Weiss-Field function
  call allocate_kb_contour_gf(Ker,cc_params)
  allocate(Gk(Lk),dGk(Lk),dGk_old(Lk))
  do ik=1,Lk
     call allocate_kb_contour_gf(Gk(ik),cc_params)
     call allocate_kb_contour_dgf(dGk(ik),cc_params)
     call allocate_kb_contour_dgf(dGk_old(ik),cc_params)
     Gk(ik)=zero
  end do
  call allocate_kb_contour_dgf(Gedge,cc_params)
  allocate(ham(cc_params%Ntime))
  allocate(nk(cc_params%Ntime,Lk))



  !READ OR GUESS THE INITIAL SELF-ENERGY
  call init_equilibrium_guess(Sigma,cc_params)


  !GET THE STARTING CONDITIONS: 
  !BUILD THE T=0(itime=1) FUNCTION FROM NON-INTERACTING SOLUTION
  Gloc = zero
  do ik=1,Lk 
     call setup_initial_conditions(Sigma,Gk(ik),dGk(ik),ik,cc_params)
     Gloc%mats(0: ) = Gloc%mats(0:)  + wt(ik)*Gk(ik)%mats(0:)
     Gloc%iw(:)     = Gloc%iw(:)     + wt(ik)*Gk(ik)%iw(:)
     Gloc%ret(1,1)  = Gloc%ret(1,1)  + wt(ik)*Gk(ik)%ret(1,1)
     Gloc%less(1,1) = Gloc%less(1,1) + wt(ik)*Gk(ik)%less(1,1)
     Gloc%lmix(1,0:)= Gloc%lmix(1,0:)+ wt(ik)*Gk(ik)%lmix(1,0:)
     nk(1,ik)=-xi*Gk(ik)%less(1,1)
  enddo


  !================================================================
  !UP TO HERE YOU HAVE ALL G_loc,G_k,G_0,Sigma in the square [1,1]
  !================================================================



  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Nt=itime

     call setup_sigma(Sigma,cc_params)
     do ik=1,Lk
        dGk_old(ik) = dGk(ik)
     enddo
     !
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        write(*,"(A,I2,A1)",advance='no')"dmft loop=",iloop," "

        !IMPURITY SOLVER: IPT.
        !GET SIGMA FROM THE WEISS FIELD
        call neq_solve_ipt(Gwf,Sigma,cc_params)
        call plot_kb_contour_gf("Sigma",Sigma,cc_params)
        !PERFORM THE SELF_CONSISTENCY: get local GF + update Weiss-Field
        !
        !prepare the kernel for evolution: Ker=Sigma+S_bath
        call add_kb_contour_gf(Sbath,Sigma,Ker,cc_params)
        !
        !sum over all k-points
        Gedge=zero
        do ik=1,Lk
           !call setup_initial_conditions(Ker,Gk(ik),dGk(ik),ik,cc_params)
           ham(:)=hamkt(ik,cc_params%ntime,cc_params)
           call vide_kb_contour_gf(Ham,Ker,Gk(ik),dGk_old(ik),dGk(ik),cc_params)
           nk(itime,ik)=-xi*Gk(ik)%less(itime,itime)
           Gedge%ret =Gedge%ret  + wt(ik)*Gk(ik)%ret(itime,:)
           Gedge%less=Gedge%less + wt(ik)*Gk(ik)%less(itime,:)
           Gedge%lmix(0:)=Gedge%lmix(0:) + wt(ik)*Gk(ik)%lmix(itime,0:)
        enddo
        Gloc%ret(itime,:)=Gedge%ret
        Gloc%less(itime,:)=Gedge%less
        Gloc%less(1:itime-1,itime)=-conjg(Gedge%less(1:itime-1))
        Gloc%lmix(itime,0:)=Gedge%lmix(0:)
        call plot_kb_contour_gf("Gloc",Gloc,cc_params)

        !
        !update the weiss field by solving the integral equation:
        ! G0 + K*G0 = Q, with K = \Sigma*G and Q = G
        call convolute_kb_contour_gf(Sigma,Gloc,Ker,cc_params)
        call vie_kb_contour_gf(Gloc,Ker,Gwf,cc_params)
        call plot_kb_contour_gf("G0",Gwf,cc_params)
        !
        !CHECK CONVERGENCE
        converged = convergence_check(Gwf,cc_params)
     enddo
  enddo


  !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
  allocate(n_t(cc_params%Ntime))
  allocate(Jloc(cc_params%Ntime))
  allocate(ndens(0:Nx,0:Nx,cc_params%Ntime))
  forall(i=1:cc_params%Ntime)n_t(i) = -xi*Gloc%less(i,i)
  Jloc=Vzero
  do ik=1,Lk
     ix=ik2ix(ik);iy=ik2iy(ik)
     do i=1,cc_params%Ntime
        Ak      = Afield(cc_params%t(i),Ek)
        kt      = kgrid(ix,iy)-Ak
        Jk      = nk(i,ik)*square_lattice_velocity(kt)
        Jloc(i) = Jloc(i) +  wt(ik)*Jk
     enddo
  enddo

  forall(i=1:cc_params%Ntime,ik=1:Lk)nDens(ik2ix(ik),ik2iy(ik),i)=nk(i,ik)
  call splot3d("3dFSVSpiVSt.plot",kgrid(0:Nx,0)%x,kgrid(0,0:Nx)%y,nDens(0:Nx,0:Nx,:))
  call splot3d("nkVSepsikVStime.plot",cc_params%t,epsik,nk)
  call splot("JlocVStime.plot",cc_params%t,Jloc%x,Jloc%y)
  call splot("ntVStime.plot",cc_params%t,n_t)


  print*,"BRAVO"



contains




  subroutine init_equilibrium_guess(self,params)
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    real(8)                 :: wm,res,ims,C0,C1,n0
    logical                 :: bool
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    complex(8)              :: zeta,hamk(1)
    if(.not.g0%status)stop "init_equilibrium_weiss_field: g0 is not allocated"
    if(.not.params%status)stop "init_equilibrium_weiss_field: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Lf
    !
    !CHECK IF G0(IW) IS AVAILABLE
    !IF NOT, START FROM NON-INTERACTING (SIGMA=0)
    inquire(file=trim(g0file),exist=bool)
    if(bool)then
       write(*,"(A)")"Reading initial Sigma(iw) from file "//reg(g0file)
       unit = free_unit()
       open(unit,file=reg(g0file),status='old')
       i = file_length(reg(g0file))
       if(i/=Lf)then
          print*,"init_equilibrium_weiss_field: Lfreq in +g0file does not correspond",i
          stop
       endif
       do i=1,Lf
          read(unit,*)wm,ims,res
          self%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       stop "Need a non-interacting Sigma(iw) to start... go get one!"
    endif
    !NOW YOU NEED TO PROPAGATE THE EQUILIBRIUM SOLUTION to SIGMA^{x=M,<,R,\lmix}
    n0=1.d0/2.d0
    C0=Uloc(1)*(n0-0.5d0)
    C1=Uloc(1)**2*n0*(1.d0-n0)
    call fftgf_iw2tau(self%iw -C1/(xi*params%wm)-C0,self%mats(0:),beta,notail=.true.)
    self%mats=self%mats-C1*0.5d0
    self%less(1,1) = -xi*self%mats(L)
    self%ret(1,1)  = xi*(self%mats(0)+self%mats(L))
    forall(i=0:L)self%lmix(1,i)=-xi*self%mats(L-i)
    !
  end subroutine init_equilibrium_guess





  subroutine setup_initial_conditions(Self,Gk,dGk,ik,params)
    type(kb_contour_gf)                 :: Gk,Self
    type(kb_contour_dgf)                :: dGk
    integer                             :: ik
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,Ltau
    real(8)                             :: nktmp
    complex(8)                          :: epsk(1)
    complex(8),allocatable,dimension(:) :: SxG
    Ltau  = params%Ntau
    epsk  = Hamkt(ik,1,params)
    Gk%iw = one/(xi*params%wm - epsk(1) - self%iw)          !get G_k(iw) 
    !
    call fftgf_iw2tau(Gk%iw,Gk%mats(0:),beta)        !get G_k(tau)
    nktmp = -Gk%mats(Ltau)                           !n(k,t=0)=-G^M_k(beta)=G^M_k(0-)
    Gk%less(1,1) =  xi*nktmp                         !get G^<_k(0,0)= xi*G^M_k(0-)
    Gk%ret(1,1)  = -xi                               !get G^R_k(0,0)=-xi
    forall(i=0:Ltau)Gk%lmix(1,i)=-xi*Gk%mats(Ltau-i) !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
    !Derivatives
    allocate(SxG(0:Ltau))
    !
    !get d/dt G_k^R = -i e(k,0)G_k^R
    dGk%ret(1)  = -xi*epsk(1)*Gk%ret(1,1)            
    !
    !get d/dt G_k^< = -i e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    do k=0,Ltau
       SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,Ltau-k))
    end do
    dGk%less(1) = -xi*epsk(1)*Gk%less(1,1)-&
         xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,Ltau) 
    !
    !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(0:)= -xi*epsk(1)*Gk%lmix(1,0:)
    do j=0,Ltau
       do k=0,j
          SxG(k)=self%lmix(1,k)*Gk%mats(k+Ltau-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       do k=j,Ltau
          SxG(k)=self%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,Ltau) 
    enddo
  end subroutine setup_initial_conditions




  subroutine setup_sigma(self,params)
    type(kb_contour_gf)                   :: self
    type(kb_contour_params)               :: params
    integer                               :: i,j,k,N,L
    if(.not.self%status)stop "setup_sigma: self is not allocated"
    if(.not.params%status)stop "setup_sigma: params is not allocated"
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
          self%ret(N,j) =self%ret(1,1)
          self%less(N,j)=self%less(1,1)
       end do
       do i=1,N-1
          self%less(i,N)=self%less(1,1)
       end do
       do j=0,L
          self%lmix(N,j)=self%lmix(1,j)
       end do

    case default
       !EXTEND SELF FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
       !USING QUADRATIC EXTRAPOLATION
       do k=1,N-1
          self%less(N,k)=2.d0*self%less(N-1,k)-self%less(N-2,k)
          self%less(k,N)=2.d0*self%less(k,N-1)-self%less(k,N-2)
       end do
       self%less(N,N)=2.d0*self%less(N-1,N-1)-self%less(N-2,N-2)
       !
       do k=0,L
          self%lmix(N,k)=2.d0*self%lmix(N-1,k)-self%lmix(N-2,k)
       end do
       !
       self%ret(N,N)=-xi
       do k=1,N-2
          self%ret(N,k)=2.d0*self%ret(N-1,k)-self%ret(N-2,k)
       end do
       self%ret(N,N-1)=0.5d0*(self%ret(N,N)+self%ret(N,N-2))
    end select
  end subroutine setup_sigma




  function Hamkt(ik,N,params) result(hk)
    integer                 :: ik,N
    type(kb_contour_params) :: params
    integer                 :: ix,iy,it
    type(vect2D)            :: kt,Ak
    complex(8)              :: hk(N)
    ix  = ik2ix(ik)
    iy  = ik2iy(ik)
    do it=1,N
       Ak = Afield(params%t(it),Ek)
       kt = kgrid(ix,iy) - Ak
       hk(it)= one*square_lattice_dispersion(kt)
    enddo
  end function Hamkt










  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check(G,params) result(converged)
    type(kb_contour_gf)                 :: G
    type(kb_contour_params)             :: params
    logical                             :: converged
    complex(8),dimension(:),allocatable :: test_func
    integer :: N,L,Ntot
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



end PROGRAM neqDMFT
