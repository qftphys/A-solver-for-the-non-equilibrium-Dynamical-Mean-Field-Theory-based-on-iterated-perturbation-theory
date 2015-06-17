program neqDMFT
  USE NEQ_DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                               :: i,j,ik,itime,iloop,ix,iy,iz,Lk,Nx
  logical                               :: converged
  real(8)                               :: ts,kx,ky,time
  character(len=16)                     :: finput
  type(kb_contour_gf)                   :: Sbath
  type(kb_contour_gf)                   :: Gloc
  type(kb_contour_gf)                   :: Gwf
  type(kb_contour_gf)                   :: Sigma
  type(kb_contour_gf)                   :: SigmaReg
  complex(8),dimension(:),allocatable   :: SigmaHF
  !
  type(kb_contour_gf)                   :: Ker
  type(kb_contour_gf),allocatable       :: Gk(:)
  type(kb_contour_dgf),allocatable      :: dGk(:),dGk_old(:)
  type(kb_contour_dgf)                  :: Gedge
  complex(8),dimension(:),allocatable   :: Ham
  !RESULTS:
  complex(8),dimension(:,:),allocatable :: Hk
  real(8),dimension(:,:,:),allocatable  :: Vk
  real(8),dimension(:),allocatable      :: kxgrid,kygrid,Wt,Epsik
  real(8),dimension(:,:,:),allocatable  :: nDens
  real(8),dimension(:,:),allocatable    :: nk


  !READ THE INPUT FILE (in vars_global):
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.conf')
  call parse_input_variable(ts,"ts",finput,default=1d0,comment="hopping")
  call parse_input_variable(Nx,"Nx",finput,default=21,comment="Number of k-points")
  call read_input_init(trim(finput))


  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,Niw)
  call setup_kb_contour_params(cc_params,dt,beta)

  !SET THE ELECTRIC FIELD (in electric_field):
  call set_efield_vector(cc_params%t)


  !BUILD THE LATTICE STRUCTURE (use tight_binding):
  Lk = Nx*Nx
  allocate(Epsik(Lk),Wt(Lk))
  allocate(Hk(Ntime,Lk),Vk(Ntime,Lk,2))
  allocate(kxgrid(Nx),kygrid(Nx))
  write(*,*) "Using Nk_total="//txtfy(Lk)
  kxgrid = kgrid(Nx)
  kygrid = kgrid(Nx)
  Epsik  = build_hk_model(hk_model,kxgrid,kygrid,[0d0])
  Wt     = 1d0/Lk
  call write_hk_w90("Hk2d.dat",1,1,0,1,dcmplx(Epsik,0d0),kxgrid,kygrid,[0d0])
  call get_free_dos(Epsik,Wt)
  do i=1,cc_params%Ntime
     do ik=1,Lk
        call indx2coord(ik,ix,iy,iz,[Nx,Nx,1])
        kx=kxgrid(ix)
        ky=kygrid(iy)
        Hk(i,ik)   = hkt_model([kx,ky],cc_params%t(i))
        Vk(i,ik,:) = vk_model([kx,ky],cc_params%t(i))
     enddo
  enddo


  !SET THE THERMOSTAT FUNCTION (in neq_thermostat):
  call allocate_kb_contour_gf(Sbath,cc_params)
  call get_thermostat_bath(Sbath,cc_params)


  !ALLOCATE ALL THE FUNCTIONS INVOLVED IN THE CALCULATION:
  allocate(SigmaHF(cc_params%Ntime)) !Self-Energy function (Hartree-Fock term)
  call allocate_kb_contour_gf(SigmaReg,cc_params) !Self-Energy function (Regular or 2nd-order term)
  call allocate_kb_contour_gf(Sigma,cc_params) !Self-Energy function
  call allocate_kb_contour_gf(Gloc,cc_params)  !Local Green's function
  call allocate_kb_contour_gf(Gwf,cc_params)   !Local Weiss-Field function
  call allocate_kb_contour_gf(Ker,cc_params)
  allocate(Gk(Lk),dGk(Lk),dGk_old(Lk))
  do ik=1,Lk
     call allocate_kb_contour_gf(Gk(ik),cc_params)
     call allocate_kb_contour_dgf(dGk(ik),cc_params)
     call allocate_kb_contour_dgf(dGk_old(ik),cc_params)
  end do
  allocate(ham(cc_params%Ntime))
  allocate(nk(cc_params%Ntime,Lk))


  !READ OR GUESS THE INITIAL WEISS FIELD G0 (t=t'=0)
  cc_params%Nt=1
  Gloc = zero
  call neq_continue_equilibirum(Gwf,Gk,dGk,Gloc,SigmaHF,SigmaReg,Sigma,epsik,wt,cc_params)
  call measure_observables(Gloc,Sigma,Gk,Hk,Wt,cc_params)
  call measure_current(Gk,Vk,Wt,cc_params)
  do ik=1,Lk
     nk(1,ik)=dimag(Gk(ik)%less(1,1))
  enddo


  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Nt=itime
     !prepare the weiss-field at this actual time_step for DMFT:
     call extrapolate_kb_contour_gf(Gwf,cc_params)
     call extrapolate_kb_contour_gf(Gloc,cc_params)
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
        sigmaHF(itime) = zero
        call neq_solve_ipt(Gwf,SigmaReg,cc_params)
        call sum_kb_contour_gf(SigmaReg,1d0,SigmaHF,1d0,Sigma,cc_params)

        !PERFORM THE SELF_CONSISTENCY: get local GF + update Weiss-Field
        !prepare the kernel for evolution: Ker=Sigma+S_bath
        call sum_kb_contour_gf(Sbath,1d0,Sigma,1d0,Ker,cc_params)
        do ik=1,Lk
           Ham = Hk(:,ik)
           call vide_kb_contour_gf(Ham+sigmaHF,Ker,Gk(ik),dGk_old(ik),dGk(ik),cc_params)
        enddo
        call sum_kb_contour_gf(Gk(:),wt(:),Gloc,cc_params)

        !update the weiss field by solving the integral equation:
        ! G0 + K*G0 = Q , with K = G*\Sigma and Q = G
        call convolute_kb_contour_gf(Gloc,SigmaReg,Ker,cc_params,dcoeff=-1.d0)
        call vie_kb_contour_gf(Gloc,Ker,Gwf,cc_params)
        !
        !CHECK CONVERGENCE
        converged = convergence_check(Gwf,cc_params)
     enddo

     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
     call measure_observables(Gloc,Sigma,Gk,Hk,Wt,cc_params)
     call measure_current(Gk,Vk,Wt,cc_params)
     forall(ik=1:Lk)nk(itime,ik)=dimag(Gk(ik)%less(itime,itime))
  enddo


  !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
  allocate(ndens(1:Nx,1:Nx,cc_params%Ntime))
  do i=1,cc_params%Ntime
     do ik=1,Lk
        ix = indx2ix(ik,[Nx,Nx,1])
        iy = indx2iy(ik,[Nx,Nx,1])
        nDens(ix,iy,i)=nk(i,ik)
     enddo
  enddo
  call splot3d("3dFSVSpiVSt.nipt",kxgrid,kygrid,nDens(:,:,:))
  call splot3d("nkVSepsikVStime.nipt",cc_params%t,epsik,nk)
  call plot_kb_contour_gf("Sigma.nipt",Sigma,cc_params)
  call plot_kb_contour_gf("Gloc.nipt",Gloc,cc_params)
  call plot_kb_contour_gf("G0.nipt",Gwf,cc_params)

  print*,"BRAVO"



contains


  function hk_model(kpoint) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky
    complex(8)           :: hk
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = -one*2d0*ts*(cos(kx)+cos(ky))
  end function hk_model

  function hkt_model(kpoint,time) result(hk)
    real(8),dimension(:) :: kpoint
    real(8),dimension(3) :: Ak
    real(8)              :: time
    integer              :: ik
    real(8)              :: kx,ky
    complex(8)           :: hk
    Ak = Afield(time)
    kx=kpoint(1) - Ak(1)
    ky=kpoint(2) - Ak(2)
    hk = hk_model([kx,ky])
  end function Hkt_Model

  function Vk_model(kpoint,time) result(vel)
    real(8),dimension(:)               :: kpoint
    real(8),dimension(3)               :: Ak
    real(8)                            :: time
    integer                            :: ik
    real(8)                            :: kx,ky
    complex(8),dimension(size(kpoint)) :: vel
    Ak = Afield(time)
    kx=kpoint(1) - Ak(1)
    ky=kpoint(2) - Ak(2)
    vel(1) = 2d0*ts*sin(kx)
    vel(2) = 2d0*ts*sin(ky)
  end function Vk_model






  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check(G,params) result(converged)
    type(kb_contour_gf)                 :: G
    type(kb_contour_params)             :: params
    logical                             :: converged
    complex(8),dimension(:),allocatable :: test_func
    integer                             :: N,L,Ntot
    real(8)                             :: err
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
    converged=check_convergence(test_func,dmft_error,Nsuccess,nloop,oerr=err)
    deallocate(test_func)
    if(isnan(err))stop "neqdmft_2dsquare_field: convergence gives NaN"
    !if(isnan(err))stop "Aborted convergence: error=NaN"
  end function convergence_check



  ! subroutine init_equilibrium_functions(g0,gk,dgk,g,self,params)
  !   type(kb_contour_gf)                 :: g0
  !   type(kb_contour_gf)                 :: gk(:)
  !   type(kb_contour_dgf)                :: dgk(size(gk))
  !   type(kb_contour_gf)                 :: g
  !   type(kb_contour_gf)                 :: self
  !   type(kb_contour_params)             :: params
  !   real(8)                             :: wm,res,ims
  !   logical                             :: bool
  !   integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
  !   complex(8)                          :: zeta
  !   complex(8)                          :: Self_gtr
  !   complex(8),allocatable,dimension(:) :: SxG
  !   real(8),dimension(:),allocatable    :: ftau,stau
  !   !
  !   Lk=size(gk)
  !   if(.not.g0%status)stop "init_functions: g0 is not allocated"
  !   if(.not.g%status)stop "init_functions: g is not allocated"
  !   do ik=1,Lk
  !      if(.not.gk(ik)%status)stop "init_functions: gk(ik) is not allocated"
  !      if(.not.dgk(ik)%status)stop "init_functions: dgk(ik) is not allocated"
  !   enddo
  !   if(.not.self%status)stop "init_functions: self is not allocated"
  !   if(.not.params%status)stop "init_functions: params is not allocated"
  !   !
  !   N = params%Nt
  !   L = params%Ntau
  !   Lf= params%Lf
  !   !
  !   !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
  !   inquire(file=trim(g0file),exist=bool)
  !   if(bool)then
  !      write(*,"(A)")"Reading initial G0(iw) from file "//reg(g0file)
  !      i = file_length(reg(g0file))
  !      if(i/=Lf)then
  !         print*,"init_equilibrium_weiss_field: Lfreq in +g0file does not correspond",i
  !         stop
  !      endif
  !      unit = free_unit()
  !      open(unit,file=reg(g0file),status='old')
  !      do i=1,Lf
  !         read(unit,*)wm,ims,res
  !         g0%iw(i) = dcmplx(res,ims)
  !      enddo
  !      close(unit)
  !   else
  !      write(*,"(A)")"Start from Non-interacting G0(iw)"
  !      do i=1,Lf
  !         wm    = pi/beta*dble(2*i-1)
  !         zeta  = dcmplx(0.d0,wm)
  !         g0%iw(i) = sum_overk_zeta(zeta,epsik,wt)
  !      enddo
  !   endif
  !   !
  !   !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
  !   allocate(ftau(0:Lf),stau(0:Lf))
  !   call fft_gf_iw2tau(g0%iw,ftau(0:),params%beta)
  !   call extract_gtau_(ftau,g0%mats)
  !   g0%less(1,1) = -xi*g0%mats(L)
  !   g0%ret(1,1)  = -xi
  !   forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
  !   !
  !   !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
  !   !(this step depends on the imp. solv.)
  !   ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
  !   ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
  !   ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
  !   ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
  !   ! self^R(0,0) = self^> - self^<
  !   do j=0,Lf
  !      stau(j) = Ui*Ui*ftau(j)*ftau(Lf-j)*ftau(j)
  !   end do
  !   call extract_gtau_(stau,Self%mats)
  !   call fft_gf_tau2iw(Self%iw,stau,beta)
  !   Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
  !   Self_gtr      =(xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
  !   do j=0,L
  !      Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
  !   end do
  !   Self%ret(1,1) = Self_gtr - Self%less(1,1)
  !   deallocate(ftau,stau)
  !   !   
  !   G%mats=0.d0
  !   G%iw = zero
  !   do ik=1,Lk 
  !      call setup_initial_conditions(self,gk(ik),dgk(ik),ik,params)
  !      G%mats(0:)  = G%mats(0:)  + wt(ik)*gk(ik)%mats(0:)
  !      G%iw(:)     = G%iw(:)     + wt(ik)*gk(ik)%iw(:)
  !      G%ret(1,1)  = G%ret(1,1)  + wt(ik)*gk(ik)%ret(1,1)
  !      G%less(1,1) = G%less(1,1) + wt(ik)*gk(ik)%less(1,1)
  !      G%lmix(1,0:)= G%lmix(1,0:)+ wt(ik)*gk(ik)%lmix(1,0:)
  !   enddo
  !   return
  ! end subroutine init_equilibrium_functions



  ! subroutine setup_initial_conditions(Self,Gk,dGk,ik,params)
  !   type(kb_contour_gf)                 :: Gk,Self
  !   type(kb_contour_dgf)                :: dGk
  !   integer                             :: ik
  !   type(kb_contour_params)             :: params
  !   integer                             :: i,j,k,Ltau,Lf
  !   complex(8),allocatable,dimension(:) :: SxG
  !   real(8),dimension(:),allocatable    :: ftau
  !   Ltau  = params%Ntau
  !   Lf    = params%Lf
  !   Gk%iw = one/(xi*params%wm - epsik(ik) - self%iw)          !get G_k(iw) 
  !   !
  !   allocate(ftau(0:Lf))
  !   call fft_gf_iw2tau(Gk%iw,ftau(0:),beta)           !get G_k(tau)
  !   call extract_gtau_(ftau,Gk%mats)
  !   Gk%less(1,1) = -xi*Gk%mats(Ltau)                 !get G^<_k(0,0)= xi*G^M_k(0-)
  !   Gk%ret(1,1)  = -xi                               !get G^R_k(0,0)=-xi
  !   forall(i=0:Ltau)Gk%lmix(1,i)=-xi*Gk%mats(Ltau-i) !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
  !   !Derivatives
  !   allocate(SxG(0:Ltau))
  !   !get d/dt G_k^R = -i*e(k,0)G_k^R
  !   dGk%ret(1)  = -xi*epsik(ik)*Gk%ret(1,1)            
  !   !get d/dt G_k^< = -i*e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
  !   do k=0,Ltau
  !      SxG(k)=Self%lmix(1,k)*conjg(Gk%lmix(1,Ltau-k))
  !   end do
  !   dGk%less(1) = -xi*epsik(ik)*Gk%less(1,1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,Ltau) 
  !   !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
  !   dGk%lmix(0:)= -xi*epsik(ik)*Gk%lmix(1,0:)
  !   do j=0,Ltau
  !      do k=0,j
  !         SxG(k)=Self%lmix(1,k)*Gk%mats(Ltau+k-j)
  !      end do
  !      dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
  !      do k=j,Ltau
  !         SxG(k)=Self%lmix(1,k)*Gk%mats(k-j)
  !      end do
  !      dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,Ltau) 
  !   enddo
  !   deallocate(SxG,ftau)
  ! end subroutine setup_initial_conditions




  ! subroutine setup_weiss_field(g0,params)
  !   type(kb_contour_gf)                   :: g0
  !   type(kb_contour_params)               :: params
  !   integer                               :: i,j,k,N,L
  !   if(.not.g0%status)stop "init_g0: g0 is not allocated"
  !   if(.not.params%status)stop "init_g0: params is not allocated"
  !   !
  !   N = params%Nt
  !   L = params%Ntau
  !   !
  !   select case(N)
  !   case(1)
  !      return

  !   case(2)
  !      !GUESS G0 AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
  !      do j=1,N
  !         g0%ret(N,j) =g0%ret(1,1)
  !         g0%less(N,j)=g0%less(1,1)
  !      end do
  !      do i=1,N-1
  !         g0%less(i,N)=g0%less(1,1)
  !      end do
  !      do j=0,L
  !         g0%lmix(N,j)=g0%lmix(1,j)
  !      end do

  !   case default
  !      !EXTEND G0 FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
  !      !USING QUADRATIC EXTRAPOLATION
  !      do k=1,N-1
  !         g0%less(N,k)=2.d0*g0%less(N-1,k)-g0%less(N-2,k)
  !         g0%less(k,N)=2.d0*g0%less(k,N-1)-g0%less(k,N-2)
  !      end do
  !      g0%less(N,N)=2.d0*g0%less(N-1,N-1)-g0%less(N-2,N-2)
  !      !
  !      do k=0,L
  !         g0%lmix(N,k)=2.d0*g0%lmix(N-1,k)-g0%lmix(N-2,k)
  !      end do
  !      !
  !      g0%ret(N,N)=-xi
  !      do k=1,N-2
  !         g0%ret(N,k)=2.d0*g0%ret(N-1,k)-g0%ret(N-2,k)
  !      end do
  !      g0%ret(N,N-1)=0.5d0*(g0%ret(N,N)+g0%ret(N,N-2))
  !   end select
  ! end subroutine setup_weiss_field




  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : Solve with the 2^nd IPT sigma functions
  ! !+-------------------------------------------------------------------+
  ! subroutine neq_solve_ipt(G0,Sigma,params)
  !   type(kb_contour_gf)                   :: G0
  !   type(kb_contour_gf)                   :: Sigma
  !   type(kb_contour_params)               :: params
  !   integer                               :: N,L
  !   complex(8),dimension(:,:),allocatable :: G0_gtr,Sigma_gtr,G0_rmix
  !   integer                               :: i,j,itau
  !   !
  !   N   = params%Nt                 !<== work with the ACTUAL size of the contour
  !   L   = params%Ntau
  !   allocate(G0_gtr(N,N),Sigma_gtr(N,N))
  !   do j=1,N
  !      G0_gtr(N,j)=G0%less(N,j)+ G0%ret(N,j)
  !   end do
  !   do i=1,N-1
  !      G0_gtr(i,N)=G0%less(i,n)-conjg(G0%ret(N,i))
  !   end do
  !   !
  !   !Vertical edge
  !   do j=1,N
  !      Sigma%less(N,j)= U*U*G0%less(N,j)*G0_gtr(j,N)*G0%less(N,j)
  !      Sigma_gtr(N,j) = U*U*G0_gtr(N,j)*G0%less(j,N)*G0_gtr(N,j)
  !   end do
  !   !Horizontal edge
  !   do i=1,N-1
  !      Sigma%less(i,N)= U*U*G0%less(i,N)*G0_gtr(N,i)*G0%less(i,N)
  !      Sigma_gtr(i,N) = U*U*G0_gtr(i,N)*G0%less(N,i)*G0_gtr(i,N)
  !   end do
  !   !Imaginary time edge:
  !   forall(i=0:L)Sigma%lmix(N,i)  = U*Ui*G0%lmix(N,i)*(conjg(G0%lmix(N,L-i)))*G0%lmix(N,i)
  !   forall(j=1:N)Sigma%ret(N,j) = Sigma_gtr(N,j) - Sigma%less(N,j)

  ! end subroutine neq_solve_ipt



















  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: measure some observables and print them
  ! !+-------------------------------------------------------------------+
  ! subroutine measure_current(gk,Vkt_,Wtk_,params)
  !   type(kb_contour_gf)                   :: gk(:)
  !   real(8),dimension(:,:,:)       :: Vkt_
  !   real(8),dimension(:)           :: wtk_
  !   type(kb_contour_params)               :: params
  !   integer                               :: unit,itime,Lk,ix,iy,ik
  !   ! type(vect2D)                          :: Ak
  !   real(8)                               :: nk,Jloc(2)!,Ak(3)
  !   Lk=size(gk)
  !   itime = params%Nt
  !   unit = free_unit()
  !   open(unit,file="current.info")
  !   write(unit,"(8A20)")"time","Jx","Jy","Jz"
  !   close(unit)
  !   open(unit,file="current.ipt",position="append")
  !   Jloc=0d0
  !   do ik=1,Lk
  !      nk = dimag(Gk(ik)%less(itime,itime))
  !      ! Ak   = Afield(cc_params%t(itime))
  !      ! call indx2coord(ik,ix,iy,iz,[Nx,Nx,1])
  !      ! kx=kxgrid(ix) - Ak(1)!%x
  !      ! ky=kygrid(iy) - Ak(2)!%y
  !      Jloc = Jloc +  wt(ik)*nk*Vk(itime,ik,:)
  !   enddo
  !   write(unit,"(4F20.12)")params%t(itime),Jloc(1),Jloc(2)
  !   close(unit)
  ! end subroutine measure_current








  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: measure some observables and print them
  ! !+-------------------------------------------------------------------+
  ! subroutine measure_observables(g,gk,self,params)
  !   type(kb_contour_gf)                   :: g
  !   type(kb_contour_gf)                   :: gk(:)
  !   type(kb_contour_gf)                   :: self
  !   type(kb_contour_params)               :: params
  !   integer                               :: unit,itime,Lk,ix,iy,ik
  !   type(vect2D)                          :: kt,Ak,Jloc
  !   real(8)                               :: dens,docc,ekin,epot,etot,nk
  !   Lk=size(gk)
  !   itime = params%Nt
  !   unit = free_unit()
  !   open(unit,file="columns.plot")
  !   write(unit,"(8A20)")"time","Jx","Docc","Jy","n","Ekin","Epot","Etot"
  !   close(unit)
  !   unit = free_unit()
  !   open(unit,file="observables.plot",position="append")
  !   dens = measure_dens(g,self,params)
  !   docc = measure_docc(g,self,params)
  !   ekin = measure_ekin(g,self,params)
  !   epot = measure_epot(g,self,params)
  !   etot = ekin + epot
  !   Jloc=Vzero
  !   do ik=1,Lk
  !      ix=ik2ix(ik)
  !      iy=ik2iy(ik)
  !      nk = dimag(Gk(ik)%less(itime,itime))
  !      Ak   = Afield(cc_params%t(itime),Ek)
  !      kt   = kgrid(ix,iy) - Ak
  !      Jloc = Jloc +  wt(ik)*nk*square_lattice_velocity(kt)
  !   enddo
  !   write(unit,"(8F20.12)")params%t(itime),Jloc%x,docc,Jloc%y,dens,ekin,epot,etot
  !   close(unit)
  ! end subroutine measure_observables



  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: return the value of the density at a given istant of time
  ! ! n(t)=-xi*G^<(t,t)
  ! !+-------------------------------------------------------------------+
  ! function measure_dens(g,self,params) result(dens)
  !   type(kb_contour_gf)                 :: g
  !   type(kb_contour_gf)                 :: self
  !   type(kb_contour_params)             :: params
  !   real(8)                             :: dens
  !   integer                             :: N
  !   N = params%Nt
  !   dens = dimag(G%less(N,N))
  ! end function measure_dens


  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: return the value of the double occupancy at a given istant of time
  ! ! d(t)=n_up(t)*n_do(t)-1/U0*[Self^M*G^M]
  ! !      n_up(t)*n_do(t)-i/U*[Self^R*G^< + Self^<*G^A + Self^\lmix*G^\rmix](t,t)
  ! !+-------------------------------------------------------------------+
  ! function measure_docc(g,self,params) result(docc)
  !   type(kb_contour_gf)                 :: g
  !   type(kb_contour_gf)                 :: self
  !   type(kb_contour_params)             :: params
  !   real(8)                             :: docc,epot
  !   integer                             :: i,k,j,N,L,Lf
  !   complex(8),dimension(:),allocatable :: SxG
  !   real(8)                             :: nt
  !   N = params%Nt
  !   L = params%Ntau
  !   Lf= params%Lf
  !   !
  !   nt   = dimag(G%less(N,N))
  !   allocate(SxG(0:max(N,L)))
  !   docc = nt**2
  !   if(N==1)then
  !      if(ui/=0.d0)then
  !         ! epot=0.d0
  !         ! do i=1,L
  !         !    epot=epot+dreal(self%iw(i)*g%iw(i))
  !         ! enddo
  !         ! epot=2.d0*epot/beta
  !         ! docc=epot/Ui + 0.5d0*nt*2.d0 - 0.25d0
  !         do k=0,L
  !            SxG(k)=Self%mats(L-k)*G%mats(k)
  !         end do
  !         docc=docc-1.d0/Ui*params%dtau*kb_trapz(SxG(0:),0,L)
  !      endif
  !   else
  !      if(u/=0.d0)then
  !         do k=0,L
  !            SxG(k)=Self%lmix(N,k)*conjg(G%lmix(N,L-k))
  !         end do
  !         docc=docc + 1.d0/U*params%dtau*dimag( (-xi)*kb_trapz(SxG(0:),0,L) )
  !         do k=1,N
  !            SxG(k)=Self%ret(N,k)*G%less(k,N)
  !         end do
  !         docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
  !         do k=1,N
  !            SxG(k)=Self%less(N,k)*conjg(G%ret(N,k))
  !         end do
  !         docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
  !      endif
  !   endif
  !   deallocate(SxG)
  ! end function measure_docc


  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  ! !+-------------------------------------------------------------------+
  ! function measure_ekin(g,self,params) result(ekin)
  !   type(kb_contour_gf)                 :: g
  !   type(kb_contour_gf)                 :: self
  !   type(kb_contour_params)             :: params
  !   real(8)                             :: ekin
  !   integer                             :: i,k,j,N,L
  !   complex(8),dimension(:),allocatable :: Ker
  !   real(8)                             :: nt
  !   N = params%Nt
  !   L = params%Ntau
  !   !
  !   allocate(Ker(0:max(N,L)))
  !   if(N==1)then
  !      do k=0,L
  !         Ker(k)=G%mats(L-k)*G%mats(k)
  !      end do
  !      ekin = -2.d0*params%dtau*kb_trapz(Ker(0:),0,L)
  !   else
  !      do k=0,L
  !         Ker(k)=G%lmix(N,k)*conjg(G%lmix(N,L-k))
  !      end do
  !      ekin=2.d0*params%dtau*dimag( (-xi)*kb_trapz(Ker(0:),0,L) )
  !      do k=1,N
  !         Ker(k)=G%ret(N,k)*G%less(k,N)
  !      end do
  !      ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
  !      do k=1,N
  !         Ker(k)=G%less(N,k)*conjg(G%ret(N,k))
  !      end do
  !      ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
  !   endif
  !   deallocate(Ker)
  ! end function measure_ekin



  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! ! U(t)= U*docc(t) - n(t) + 1/4
  ! !+-------------------------------------------------------------------+
  ! function measure_epot(g,self,params) result(epot)
  !   type(kb_contour_gf)                 :: g
  !   type(kb_contour_gf)                 :: self
  !   type(kb_contour_params)             :: params
  !   real(8)                             :: epot,docc,nt
  !   integer                             :: i,k,j,N,L
  !   N = params%Nt
  !   L = params%Ntau
  !   !
  !   if(N==1)then
  !      nt   = measure_dens(g,self,params)
  !      docc = measure_docc(g,self,params)
  !      epot = Ui*(docc - nt + 0.25d0)
  !   else
  !      nt   = measure_dens(g,self,params)
  !      docc = measure_docc(g,self,params)
  !      epot = U*(docc - nt + 0.25d0)
  !   endif
  ! end function measure_epot



  ! subroutine extract_gtau_(g,gred)
  !   real(8),dimension(0:) :: g
  !   real(8),dimension(0:) :: gred
  !   integer               :: N,Nred
  !   integer               :: i,ip
  !   real(8)               :: p,mismatch
  !   N   =size(g)-1
  !   Nred=size(gred)-1
  !   gred(0)=g(0)
  !   ! if(g(0) > 0.d0) gred(Nred)=1.d0-g(0)
  !   ! if(g(0) < 0.d0) gred(Nred)=-(g(0)+1.d0)
  !   gred(Nred)=g(N)
  !   mismatch=dble(N)/dble(Nred)
  !   do i=1,Nred-1
  !      p=dble(i)*mismatch
  !      ip=int(p)
  !      gred(i)=g(ip)
  !   enddo
  ! end subroutine extract_gtau_


end PROGRAM neqDMFT
