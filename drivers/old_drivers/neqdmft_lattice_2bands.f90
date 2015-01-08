!###################################################################
!PURPOSE  : Solve quantum quench dynamic in the 2-band lattice (2Dsquare)
!AUTHORS  : Adriano Amaricci 
!###################################################################
program neqDMFT
  USE ERROR
  USE FFTGF
  USE MATRIX
  USE NEQ_VARS_GLOBAL   !global variables, calls to 3rd library 
  USE NEQ_THERMOSTAT    !contains bath inizialization
  USE SQUARE_LATTICE
  implicit none
  !Number of the bands: STATIC!
  integer,parameter :: Nb=2
  !
  integer                             :: i,j,it,ib,jb,ik,itime,iloop,ix,iy
  logical                             :: converged
  character(len=16)                   :: finput
  type(kb_contour_gf),allocatable     :: Gwf(:)
  type(kb_contour_gf),allocatable     :: Sigma(:)
  type(kb_contour_gf),allocatable     :: Gloc(:)
  type(kb_contour_gf)                 :: Ker
  type(kb_contour_dgf),allocatable    :: Sedge(:,:),Gedge(:)
  type(kb_contour_gf),allocatable     :: Gk(:,:,:),Fk(:,:)
  type(kb_contour_dgf),allocatable    :: dGk(:,:,:),dGk_old(:,:,:)
  type(kb_contour_dgf),allocatable    :: dFk(:,:),dFk_old(:,:)
  complex(8),dimension(:),allocatable :: Ham
  real(8),dimension(:,:),allocatable  :: nk
  !INTERNAL VARIABLES:
  real(8)                             :: tpd
  !
  !READ THE INPUT FILE (in vars_global):
  !=====================================================================
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.in')
  call parse_cmd_variable(ts,"TS",default=0.5d0)
  call parse_cmd_variable(tpd,"TPD",default=0.d0)
  call read_input_init(trim(finput))

  !BUILD THE LATTICE STRUCTURE (in lib/square_lattice):
  Lk   = square_lattice_dimension(Nx,Nx)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)
  call get_free_dos(epsik,wt)


  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  !=====================================================================
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,dt,beta,Lfreq)
  cc_params%t = linspace(0.d0,cc_params%tmax,cc_params%Ntime,mesh=dt)
  cc_params%tau(0:) = linspace(0.d0,cc_params%beta,cc_params%Ntau+1,mesh=dtau)
  cc_params%wm  = pi/cc_params%beta*dble(2*arange(1,cc_params%Lf)-1)
  print*,"dt=",dt,cc_params%dt
  print*,"dtau=",dtau,cc_params%dtau


  !ALLOCATE ALL THE FUNCTIONS INVOLVED IN THE CALCULATION:
  !=====================================================================
  allocate(Sigma(Nb),Gloc(Nb),Gwf(Nb))
  allocate(Sedge(Nb,2),Gedge(Nb))
  allocate(Gk(Lk,Nb,Nb),dGk(Lk,Nb,Nb),dGk_old(Lk,Nb,Nb))
  allocate(Fk(Lk,Nb),dFk(Lk,Nb),dFk_old(Lk,Nb))

  do ib=1,Nb
     call allocate_kb_contour_gf(Sigma(ib),cc_params) !Self-Energy function
     call allocate_kb_contour_dgf(Sedge(ib,1),cc_params,wgtr=.true.)
     call allocate_kb_contour_dgf(Sedge(ib,2),cc_params,wgtr=.true.)
     call allocate_kb_contour_gf(Gwf(ib),cc_params)   !Local Weiss-Field function
     call allocate_kb_contour_dgf(Gedge(ib),cc_params)
     call allocate_kb_contour_gf(Gloc(ib),cc_params)  !Local Green's function
     Gloc(ib) = zero
  enddo
  do ik=1,Lk
     do ib=1,Nb
        do jb=1,Nb
           call allocate_kb_contour_gf(Gk(ik,ib,jb),cc_params)
           call allocate_kb_contour_dgf(dGk(ik,ib,jb),cc_params)
           call allocate_kb_contour_dgf(dGk_old(ik,ib,jb),cc_params)
           Gk(ik,ib,jb)=zero
           dGk(ik,ib,jb)=zero
           dGk_old(ik,ib,jb)=zero
        enddo
        call allocate_kb_contour_gf(Fk(ik,ib),cc_params)
        call allocate_kb_contour_dgf(dFk(ik,ib),cc_params)
        call allocate_kb_contour_dgf(dFk_old(ik,ib),cc_params)
        Fk(ik,ib)=zero
        dFk(ik,ib)=zero
        dFk_old(ik,ib)=zero
     enddo
  enddo
  call allocate_kb_contour_gf(Ker,cc_params)
  allocate(ham(cc_params%Ntime))


  !READ OR GUESS THE INITIAL WEISS FIELD G0 (t=t'=0)
  !=====================================================================
  cc_params%Nt=1
  call init_equilibrium_functions(Gwf,Gk,dGk,Fk,dFk,Gloc,Sigma,cc_params)
  call measure_observables(Gloc,Sigma,Sedge,cc_params)


  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  !=====================================================================
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Nt=itime
     !
     call setup_weiss_field(Gwf,cc_params)
     do ik=1,Lk
        do ib=1,Nb
           do jb=1,Nb
              dGk_old(ik,ib,jb) = dGk(ik,ib,jb)
           enddo
           dFk_old(ik,ib) = dFk(ik,ib)
        enddo
     enddo
     !
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        write(*,"(A,I2,A1)",advance='no')"dmft loop=",iloop," "
        !
        !IMPURITY SOLVER: IPT.
        !GET SIGMA FROM THE WEISS FIELD
        call neq_solve_ipt(Gwf,Sigma,cc_params)
        !PERFORM THE SELF_CONSISTENCY: get local GF + update Weiss-Field
        do ib=1,Nb
           Gedge(ib)=zero
        enddo
        do ik=1,Lk
           !1) get f_{a,k}:
           do ib=1,Nb
              do it=1,cc_params%Nt
                 ham(it)=Hamk(ik,ib,ib,cc_params)
              enddo
              call vide_kb_contour_gf(Ham(:),Sigma(ib),Fk(ik,ib),dFk_old(ik,ib),dFk(ik,ib),cc_params)
           enddo
           !2) get K_a = f_{a,k}*f_{abar,k}
           !3) get G_{aa,k} = K_a*G_{aa,k} + f_{a,k}
           do ib=1,Nb
              call convolute_kb_contour_gf(Fk(ik,ib),Fk(ik,3-ib),Ker,cc_params,dcoeff=dreal(Hamk(ik,ib,3-ib,cc_params))**2)
              call vie_kb_contour_gf(Fk(ik,ib),Ker,Gk(ik,ib,ib),cc_params)
           enddo
           !
           do ib=1,Nb
              Gedge(ib)%ret(1:itime) = Gedge(ib)%ret(1:itime)  + wt(ik)*Gk(ik,ib,ib)%ret(itime,1:itime)
              Gedge(ib)%less(1:itime)= Gedge(ib)%less(1:itime) + wt(ik)*Gk(ik,ib,ib)%less(itime,1:itime)
              Gedge(ib)%lmix(0:)     = Gedge(ib)%lmix(0:)      + wt(ik)*Gk(ik,ib,ib)%lmix(itime,0:)
           enddo
        enddo
        do ib=1,Nb
           Gloc(ib)%ret(itime,1:itime)   = Gedge(ib)%ret(1:itime)
           Gloc(ib)%less(itime,1:itime)  = Gedge(ib)%less(1:itime)
           Gloc(ib)%lmix(itime,0:)       = Gedge(ib)%lmix(0:)
           Gloc(ib)%less(1:itime-1,itime)=-conjg(Gedge(ib)%less(1:itime-1))
        enddo

        !update the weiss field by solving the integral equation:
        ! G0 + K*G0 = Q , with K = G*\Sigma and Q = G
        do ib=1,2
           call convolute_kb_contour_gf(Gloc(ib),Sigma(ib),Ker,cc_params,dcoeff=-1.d0)
           call vie_kb_contour_gf(Gloc(ib),Ker,Gwf(ib),cc_params)
        enddo
        !
        !CHECK CONVERGENCE
        converged = convergence_check(Gwf,cc_params)
        !
     enddo
     do ib=1,Nb
        call plot_kb_contour_gf("Sigma_l"//txtfy(ib),Sigma(ib),cc_params)
        call plot_kb_contour_gf("G0_l"//txtfy(ib),Gwf(ib),cc_params)
        call plot_kb_contour_gf("Gloc_l"//txtfy(ib),Gloc(ib),cc_params)
     enddo
     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
     call measure_observables(Gloc,Sigma,Sedge,cc_params)
  enddo

  print*,"BRAVO"

contains



  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine init_equilibrium_functions(g0,gk,dgk,fk,dfk,g,self,params)
    !subroutine init_equilibrium_functions(g0,dg0,g,self,params)
    type(kb_contour_gf)                 :: g0(:)
    type(kb_contour_gf)                 :: gk(:,:,:),fk(:,:)
    type(kb_contour_dgf)                :: dgk(:,:,:),dfk(:,:)
    type(kb_contour_gf)                 :: g(:)
    type(kb_contour_gf)                 :: self(:)
    type(kb_contour_params)             :: params
    integer                             :: Nb
    real(8)                             :: wm,res(2),ims(2)
    logical                             :: bool
    integer                             :: i,j,ib,jb,k,ik,unit,len,N,L,Lf
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: SxG
    Lk=size(gk,1)
    if(.not.g0(1)%status)stop "init_functions: g0 is not allocated"
    do ik=1,Lk
       if(.not.gk(ik,1,1)%status)stop "init_functions: gk(ik) is not allocated"
       if(.not.dgk(ik,1,1)%status)stop "init_functions: dgk(ik) is not allocated"
       if(.not.fk(ik,1)%status)stop "init_functions: gk(ik) is not allocated"
       if(.not.dfk(ik,1)%status)stop "init_functions: dgk(ik) is not allocated"
    enddo
    if(.not.g(1)%status)stop "init_functions: g is not allocated"
    if(.not.self(1)%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Lf
    Nb= size(g0);if(Nb>2)stop "Nb>2 is not yet implemented"
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
          read(unit,*)wm,ims(1),res(1),ims(2),res(2)
          g0(1)%iw(i) = dcmplx(res(1),ims(1))
          g0(2)%iw(i) = dcmplx(res(2),ims(2))
       enddo
       close(unit)
    else
       write(*,"(A)")"Start from Non-interacting G0(iw)"
       do i=1,Lf
          wm    = pi/beta*dble(2*i-1)
          zeta  = dcmplx(0.d0,wm)
          g0(1)%iw(i) = sum_overk_zeta(zeta,epsik,wt)
          g0(2)%iw(i) = sum_overk_zeta(zeta,epsik,wt)
       enddo
    endif
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    do ib=1,Nb
       call fftgf_iw2tau(g0(ib)%iw,g0(ib)%mats(0:),params%beta)
       g0(ib)%less(1,1) = -xi*g0(ib)%mats(L)
       g0(ib)%ret(1,1)  = -xi
       forall(i=0:L)g0(ib)%lmix(1,i)=-xi*g0(ib)%mats(L-i)
    enddo
    !
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    !(this step depends on the imp. solv.)
    ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
    ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
    ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
    ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
    ! self^R(0,0) = self^> - self^<
    do ib=1,Nb
       do j=0,L
          Self(ib)%mats(j) = Ui*Ui*g0(ib)%mats(j)*g0(ib)%mats(L-j)*g0(ib)%mats(j)+&
               2.d0*Usti*Usti*g0(ib)%mats(j)*g0(3-ib)%mats(L-j)*g0(3-ib)%mats(j)
       end do
       call fftgf_tau2iw(Self(ib)%mats(0:),Self(ib)%iw,beta)
       Self(ib)%iw = xi*dimag(Self(ib)%iw) !!ACTHUNG: imposing half-filling symmetry
       !
       Self(ib)%less(1,1)=(xi**3)*U*U*g0(ib)%mats(L)*g0(ib)%mats(0)*g0(ib)%mats(L)+&
            2.d0*(xi**3)*Ust*Ust*g0(ib)%mats(L)*g0(3-ib)%mats(0)*g0(3-ib)%mats(L)
       !
       Self_gtr          =(xi**3)*U*U*g0(ib)%mats(0)*g0(ib)%mats(L)*g0(ib)%mats(0)+&
            2.d0*(xi**3)*Ust*Ust*g0(ib)%mats(0)*g0(3-ib)%mats(L)*g0(3-ib)%mats(0)
       !
       do j=0,L
          Self(ib)%lmix(1,j)=(xi**3)*U*Ui*g0(ib)%mats(L-j)*g0(ib)%mats(j)*g0(ib)%mats(L-j)+&
               2.d0*(xi**3)*Ust*Usti*g0(ib)%mats(L-j)*g0(3-ib)%mats(j)*g0(3-ib)%mats(L-j)
       end do
       Self(ib)%ret(1,1) = Self_gtr - Self(ib)%less(1,1)
    enddo
    !
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    ! 
    do ik=1,Lk 
       call setup_initial_conditions(self,gk(ik,:,:),dgk(ik,:,:),fk(ik,:),dfk(ik,:),ik,params)
       if(mod(ik,Lk/10)==0)print*,"ik=",ik,"/",Lk
       do ib=1,Nb
          G(ib)%mats(0:)  = G(ib)%mats(0:)  + wt(ik)*gk(ik,ib,ib)%mats(0:)
          G(ib)%iw(:)     = G(ib)%iw(:)     + wt(ik)*gk(ik,ib,ib)%iw(:)
          G(ib)%ret(1,1)  = G(ib)%ret(1,1)  + wt(ik)*gk(ik,ib,ib)%ret(1,1)
          G(ib)%less(1,1) = G(ib)%less(1,1) + wt(ik)*gk(ik,ib,ib)%less(1,1)
          G(ib)%lmix(1,0:)= G(ib)%lmix(1,0:)+ wt(ik)*gk(ik,ib,ib)%lmix(1,0:)
       enddo
    enddo
    return
  end subroutine init_equilibrium_functions




  subroutine setup_initial_conditions(Self,Gk,dGk,Fk,dFk,ik,params)
    type(kb_contour_gf)                 :: Self(:)
    type(kb_contour_gf)                 :: Gk(:,:),Fk(:)
    type(kb_contour_dgf)                :: dGk(:,:),dFk(:)
    integer                             :: ik
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,Ltau,ib,jb,Nb,Lf
    real(8)                             :: nktmp
    complex(8)                          :: iw,Gkmat(size(self),size(self))
    complex(8),allocatable,dimension(:) :: SxG
    Nb=size(self);if(Nb>2)stop "setup_initial_conditions: Nb>2 not yet implemented..."
    Ltau  = params%Ntau
    Lf    = params%Lf
    do i=1,Lf                   !get G_k(iw) 
       do ib=1,Nb
          Gkmat(ib,ib) = xi*params%wm(i) - self(ib)%iw(i) - Hamk(ik,ib,ib,params)
          do jb=ib+1,Nb
             Gkmat(ib,jb) = -Hamk(ik,ib,jb,params)
             Gkmat(jb,ib) = -Hamk(ik,jb,ib,params)
          enddo
       enddo
       call matrix_inverse(Gkmat)
       do ib=1,Nb
          do jb=1,Nb
             Gk(ib,jb)%iw(i) = Gkmat(ib,jb)
          enddo
          Fk(ib)%iw(i) = one/(xi*params%wm(i) - Hamk(ik,ib,ib,params) - self(ib)%iw(i))
       enddo
    enddo
    !
    do ib=1,Nb
       do jb=1,Nb
          Gk(ib,jb)%ret(1,1)  = zero  !ACTHUNG!!
       enddo
    enddo
    do ib=1,Nb
       do jb=1,Nb
          call fftgf_iw2tau(Gk(ib,jb)%iw(:),Gk(ib,jb)%mats(0:),beta)
          !<DEBUG
          if(sum(Gk(ib,jb)%iw(:))==zero)Gk(ib,jb)%mats(0:)=0.d0
          !>DEBUG
          Gk(ib,jb)%less(1,1) = -xi*Gk(ib,jb)%mats(Ltau)                    
          do i=0,Ltau
             Gk(ib,jb)%lmix(1,i)=-xi*Gk(ib,jb)%mats(Ltau-i)
          enddo
       enddo
       Gk(ib,ib)%ret(1,1)  = -xi      
       !
       call fftgf_iw2tau(Fk(ib)%iw(:),Fk(ib)%mats(0:),beta)
       Fk(ib)%less(1,1) = -xi*Fk(ib)%mats(Ltau)                    
       do i=0,Ltau
          Fk(ib)%lmix(1,i)=-xi*Fk(ib)%mats(Ltau-i)
       enddo
       Fk(ib)%ret(1,1)  = -xi      
    enddo


    !Derivatives
    allocate(SxG(0:Ltau))
    !get d/dt G^{ab}_k^R = -i*sum_c Hk^{ac}(k,0)G^{cb}_k^R
    do ib=1,Nb
       do jb=1,Nb
          dGk(ib,jb)%ret(1)  = -xi*(Hamk(ik,ib,1,params)*Gk(1,jb)%ret(1,1)+Hamk(ik,ib,2,params)*Gk(2,jb)%ret(1,1))
       enddo
    enddo
    !
    !get d/dt G^{ab}_k^< = -i*sum_c Hk^{ac}(k,0)G^{cb}_k^<  -xi(-xi)int_0^beta S^{aa}^\lmix*G^{ab}_k^\rmix
    do ib=1,Nb
       do jb=1,Nb
          dGk(ib,jb)%less(1) = -xi*(Hamk(ik,ib,1,params)*Gk(1,jb)%less(1,1)+Hamk(ik,ib,2,params)*Gk(2,jb)%less(1,1))
          SxG=zero
          do k=0,Ltau
             SxG(k)=Self(ib)%lmix(1,k)*conjg(Gk(ib,jb)%lmix(1,Ltau-k))
          end do
          dGk(ib,jb)%less(1) = dGk(ib,jb)%less(1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,Ltau) 
       enddo
    enddo
    !
    !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    do ib=1,Nb
       do jb=1,Nb
          do j=0,Ltau
             dGk(ib,jb)%lmix(j)= -xi*(Hamk(ik,ib,1,params)*Gk(1,jb)%lmix(1,j)+Hamk(ik,ib,2,params)*Gk(2,jb)%lmix(1,j))
             SxG=zero
             do k=0,j
                SxG(k)=Self(ib)%lmix(1,k)*Gk(ib,jb)%mats(Ltau+k-j)
             end do
             dGk(ib,jb)%lmix(j)=dGk(ib,jb)%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
             do k=j,Ltau
                SxG(k)=Self(ib)%lmix(1,k)*Gk(ib,jb)%mats(k-j)
             end do
             dGk(ib,jb)%lmix(j)=dGk(ib,jb)%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,Ltau) 
          enddo
       enddo
    enddo
    !
    !
    !Get dFk:
    !get d/dt F_k^R = -i*e(k,0)F_k^R
    do ib=1,Nb
       dFk(ib)%ret(1)  = -xi*Hamk(ik,ib,ib,params)*Fk(ib)%ret(1,1)
       !get d/dt G_k^< = -i*e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
       dFk(ib)%less(1) = -xi*Hamk(ik,ib,ib,params)*Fk(ib)%less(1,1)
       do k=0,Ltau
          SxG(k)=Self(ib)%lmix(1,k)*conjg(Fk(ib)%lmix(1,Ltau-k))
       end do
       dFk(ib)%less(1)=dFk(ib)%less(1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,Ltau) 
       !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
       dFk(ib)%lmix(0:)= -xi*Hamk(ik,ib,ib,params)*Fk(ib)%lmix(1,0:)
       do j=0,Ltau
          do k=0,j
             SxG(k)=Self(ib)%lmix(1,k)*Fk(ib)%mats(Ltau+k-j)
          end do
          dFk(ib)%lmix(j)=dFk(ib)%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
          do k=j,Ltau
             SxG(k)=Self(ib)%lmix(1,k)*Fk(ib)%mats(k-j)
          end do
          dFk(ib)%lmix(j)=dFk(ib)%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,Ltau) 
       enddo
    enddo
  end subroutine setup_initial_conditions






  !+-------------------------------------------------------------------+
  !PURPOSE: setup the Weiss Field G0 for the next time-step
  !+-------------------------------------------------------------------+
  subroutine setup_weiss_field(g0,params)
    type(kb_contour_gf)                   :: g0(:)
    type(kb_contour_params)               :: params
    integer                               :: i,ib,Nb,j,k,N,L
    if(.not.g0(1)%status)stop "init_g0: g0 is not allocated"
    if(.not.params%status)stop "init_g0: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Nb=size(g0);if(Nb>2)stop "Nb>2 not yet implemented... setup_weiss_field" 
    !
    select case(N)
    case(1)
       return
    case(2)
       !GUESS G0 AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
       do ib=1,Nb
          do j=1,N
             g0(ib)%ret(N,j) =g0(ib)%ret(1,1)
             g0(ib)%less(N,j)=g0(ib)%less(1,1)
          end do
          do i=1,N-1
             g0(ib)%less(i,N)=g0(ib)%less(1,1)
          end do
          do j=0,L
             g0(ib)%lmix(N,j)=g0(ib)%lmix(1,j)
          end do
       enddo
    case default
       !EXTEND G0 FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
       !USING QUADRATIC EXTRAPOLATION
       do ib=1,Nb
          do k=1,N-1
             g0(ib)%less(N,k)=2.d0*g0(ib)%less(N-1,k)-g0(ib)%less(N-2,k)
             g0(ib)%less(k,N)=2.d0*g0(ib)%less(k,N-1)-g0(ib)%less(k,N-2)
          end do
          g0(ib)%less(N,N)=2.d0*g0(ib)%less(N-1,N-1)-g0(ib)%less(N-2,N-2)
          !
          do k=0,L
             g0(ib)%lmix(N,k)=2.d0*g0(ib)%lmix(N-1,k)-g0(ib)%lmix(N-2,k)
          end do
          !
          g0(ib)%ret(N,N)=-xi
          do k=1,N-2
             g0(ib)%ret(N,k)=2.d0*g0(ib)%ret(N-1,k)-g0(ib)%ret(N-2,k)
          end do
          g0(ib)%ret(N,N-1)=0.5d0*(g0(ib)%ret(N,N)+g0(ib)%ret(N,N-2))
       enddo
    end select
  end subroutine setup_weiss_field



  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve with the 2^nd IPT sigma functions
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt(G0,Sigma,params)
    type(kb_contour_gf)                     :: G0(:)
    type(kb_contour_gf)                     :: Sigma(size(G0))
    type(kb_contour_params)                 :: params
    integer                                 :: N,L,Nb
    complex(8),dimension(:,:,:),allocatable :: G0_gtr,Sigma_gtr,G0_rmix
    integer                                 :: i,j,k,itau,ib
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    Nb  = size(G0);if(Nb>2)stop "neq_solve_ipt: Nb>2 not yet implemented..."
    !
    allocate(G0_gtr(2,N,N),Sigma_gtr(2,N,N),G0_rmix(2,0:L,N))
    do ib=1,Nb
       do j=1,N
          G0_gtr(ib,N,j)=G0(ib)%less(N,j)+G0(ib)%ret(N,j)
       end do
       do i=1,N-1
          G0_gtr(ib,i,N)=G0(ib)%less(i,n)-conjg(G0(ib)%ret(N,i))
       end do
       do j=0,L
          G0_rmix(ib,j,N)  = conjg(G0(ib)%lmix(N,L-j))
       enddo
    enddo
    !

    !Vertical edge
    do i=1,Nb
       do j=1,Nb
          Sedge(i,j)=zero
          Sedge(i,j)%gtr=zero
       enddo
    enddo

    do ib=1,Nb
       do j=1,N
          Sedge(ib,1)%less(j)=U*U*G0(ib)%less(N,j)*G0_gtr(ib,j,N)*G0(ib)%less(N,j)
          Sedge(ib,2)%less(j)=2.d0*Ust*Ust*G0(ib)%less(N,j)*G0_gtr(3-ib,j,N)*G0(3-ib)%less(N,j)
          !
          Sedge(ib,1)%gtr(j) =U*U*G0_gtr(ib,N,j)*G0(ib)%less(j,N)*G0_gtr(ib,N,j)
          Sedge(ib,2)%gtr(j) =2.d0*Ust*Ust*G0_gtr(ib,N,j)*G0(3-ib)%less(j,N)*G0_gtr(3-ib,N,j)
          !
          Sigma(ib)%less(N,j)= U*U*G0(ib)%less(N,j)*G0_gtr(ib,j,N)*G0(ib)%less(N,j)+&
               2.d0*Ust*Ust*G0(ib)%less(N,j)*G0_gtr(3-ib,j,N)*G0(3-ib)%less(N,j)
          !
          Sigma_gtr(ib,N,j) = U*U*G0_gtr(ib,N,j)*G0(ib)%less(j,N)*G0_gtr(ib,N,j)+&
               2.d0*Ust*Ust*G0_gtr(ib,N,j)*G0(3-ib)%less(j,N)*G0_gtr(3-ib,N,j)
       end do
       !Horizontal edge
       do i=1,N-1
          Sigma(ib)%less(i,N)= U*U*G0(ib)%less(i,N)*G0_gtr(ib,N,i)*G0(ib)%less(i,N)+&
               2.d0*Ust*Ust*G0(ib)%less(i,N)*G0_gtr(3-ib,N,i)*G0(3-ib)%less(i,N)
          !
          Sigma_gtr(ib,i,N) = U*U*G0_gtr(ib,i,N)*G0(ib)%less(N,i)*G0_gtr(ib,i,N)+&
               2.d0*Ust*Ust*G0_gtr(ib,i,N)*G0(3-ib)%less(N,i)*G0_gtr(3-ib,i,N)
       end do
       !Imaginary time edge:
       do i=0,L
          Sedge(ib,1)%lmix(i)=U*Ui*G0(ib)%lmix(N,i)*G0_rmix(ib,i,N)*G0(ib)%lmix(N,i)
          Sedge(ib,2)%lmix(i)=2.d0*Ust*Usti*G0(ib)%lmix(N,i)*G0_rmix(3-ib,i,N)*G0(3-ib)%lmix(N,i)
          !
          Sigma(ib)%lmix(N,i)  = U*Ui*G0(ib)%lmix(N,i)*G0_rmix(ib,i,N)*G0(ib)%lmix(N,i)+&
               2.d0*Ust*Usti*G0(ib)%lmix(N,i)*G0_rmix(3-ib,i,N)*G0(3-ib)%lmix(N,i)
       enddo
       !Retarded
       forall(j=1:N)
          Sigma(ib)%ret(N,j) = Sigma_gtr(ib,N,j) - Sigma(ib)%less(N,j)
          Sedge(ib,1)%ret(j) = Sedge(ib,1)%gtr(j)-Sedge(ib,1)%less(j)
          Sedge(ib,2)%ret(j) = Sedge(ib,2)%gtr(j)-Sedge(ib,2)%less(j)
       end forall
    enddo

  end subroutine neq_solve_ipt






  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check(G,params) result(converged)
    type(kb_contour_gf)                 :: G(:)
    type(kb_contour_params)             :: params
    logical                             :: converged
    complex(8),dimension(:,:),allocatable :: test_func_tmp
    complex(8),dimension(:),allocatable :: test_func
    integer :: i,N,L,Ntot,ib,Nb
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    Nb  = size(G);if(Nb>2)stop "Nb>2 is not yet implemented... convergence_check"
    !
    Ntot=2*N+L+1
    allocate(test_func(2*Ntot),test_func_tmp(Nb,Ntot))
    test_func=zero
    test_func_tmp=zero
    do ib=1,Nb
       do i=1,N
          test_func_tmp(ib,i)  = G(ib)%ret(N,i)
          test_func_tmp(ib,N+i)= G(ib)%less(N,i)
       enddo
       do i=0,L
          test_func_tmp(ib,2*N+i+1)=G(ib)%lmix(N,i)
       enddo
    enddo
    test_func(1:Ntot)=test_func_tmp(1,:)
    test_func(Ntot+1:2*Ntot)=test_func_tmp(2,:)
    converged=check_convergence(test_func,dmft_error,Nsuccess,nloop)
    deallocate(test_func)
    !if(isnan(err))stop "Aborted convergence: error=NaN"
  end function convergence_check



  !+-------------------------------------------------------------------+
  !PURPOSE: measure some observables and print them
  !+-------------------------------------------------------------------+
  subroutine measure_observables(g,self,sedge,params)
    type(kb_contour_gf)     :: g(:)
    type(kb_contour_gf)     :: self(size(g))
    type(kb_contour_dgf)    :: sedge(2,size(g))    
    type(kb_contour_params) :: params
    integer                 :: unit,itime,ib,Nb,ic
    real(8)                 :: dens(size(g)),docc(size(g),size(g))
    real(8)                 :: ekin,epot,etot
    itime = params%Nt
    unit = free_unit()
    open(unit,file="observables.plot",position="append")
    dens = measure_dens(g,params)
    docc = measure_docc(g,sedge,params)
    ekin = measure_ekin(g,self,params)
    epot = measure_epot(g,sedge,params)
    etot = ekin + epot
    write(unit,"(20F18.9)")params%t(itime),dens(1),docc(1,1),dens(2),docc(2,1),docc(1,2),docc(2,2),ekin,epot,etot
    close(unit)
  end subroutine measure_observables



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the density at a given istant of time
  ! n(t)=-xi*G^<(t,t)
  !+-------------------------------------------------------------------+
  function measure_dens(g,params) result(dens)
    type(kb_contour_gf)     :: g(:)
    type(kb_contour_params) :: params
    real(8)                 :: dens(size(g))
    integer                 :: N,ib
    N = params%Nt
    do ib=1,size(g)
       dens(ib) = dimag(G(ib)%less(N,N))
    enddo
  end function measure_dens


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the double occupancy at a given istant of time
  ! d(t)=n_up(t)*n_do(t)-1/U0*[Self^M*G^M]
  !      n_up(t)*n_do(t)-i/U*[Self^R*G^< + Self^<*G^A + Self^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_docc(g,sedge,params) result(docc)
    type(kb_contour_gf)                   :: g(:)
    type(kb_contour_dgf)                  :: sedge(2,size(g))
    type(kb_contour_params)               :: params
    real(8)                               :: docc(2,size(g))
    integer                               :: i,k,j,N,L,ib,Nb
    complex(8),dimension(:,:),allocatable :: SxG
    real(8)                               :: nt(size(g))
    N = params%Nt
    L = params%Ntau
    Nb= size(g)
    !
    do ib=1,Nb
       nt(ib)    = dimag(G(ib)%less(N,N))
    enddo
    do ib=1,Nb
       docc(ib,1)= nt(ib)*nt(ib)
       docc(ib,2)= nt(ib)*nt(3-ib)
    enddo
    !
    !evaluate the temporary value of the double occupancy:
    allocate(SxG(2,0:max(L,N)))
    do ib=1,Nb
       do k=0,L
          SxG(1,k)=Sedge(ib,1)%lmix(k)*conjg(G(ib)%lmix(N,L-k))
          SxG(2,k)=Sedge(ib,2)%lmix(k)*conjg(G(ib)%lmix(N,L-k))
       end do
       if(U/=0)docc(ib,1)=docc(ib,1) + 1.d0/U*params%dtau*dimag( (-xi)*kb_trapz(SxG(1,0:),0,L) )
       if(Ust/=0)docc(ib,2)=docc(ib,2) + 0.5d0/Ust*params%dtau*dimag( (-xi)*kb_trapz(SxG(2,0:),0,L) )
       do k=1,N
          SxG(1,k)=Sedge(ib,1)%ret(k)*G(ib)%less(k,N)
          SxG(2,k)=Sedge(ib,2)%ret(k)*G(ib)%less(k,N)
       end do
       if(U/=0)docc(ib,1)=docc(ib,1) + 1.d0/U*params%dt*dimag(kb_trapz(SxG(1,0:),1,N))
       if(Ust/=0)docc(ib,2)=docc(ib,2) + 0.5d0/Ust*params%dt*dimag(kb_trapz(SxG(2,0:),1,N))
       do k=1,N
          SxG(1,k)=Sedge(ib,1)%less(k)*conjg(G(ib)%ret(N,k))
          SxG(2,k)=Sedge(ib,2)%less(k)*conjg(G(ib)%ret(N,k))
       end do
       if(U/=0)docc(ib,1)=docc(ib,1) + 1.d0/U*params%dt*dimag(kb_trapz(SxG(1,0:),1,N))
       if(Ust/=0)docc(ib,2)=docc(ib,2) + 0.5d0/Ust*params%dt*dimag(kb_trapz(SxG(2,0:),1,N))
    enddo
    deallocate(SxG)
  end function measure_docc


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_ekin(g,self,params) result(ekin)
    type(kb_contour_gf)                 :: g(:)
    type(kb_contour_gf)                 :: self(size(g))
    type(kb_contour_params)             :: params
    real(8)                             :: ekin_(size(g)),ekin
    integer                             :: i,k,j,N,L,ib,Nb
    complex(8),dimension(:),allocatable :: Ker
    N = params%Nt
    L = params%Ntau
    Nb=size(g)
    !
    allocate(Ker(0:max(N,L)))
    ekin_=0.d0
    if(N==1)then
       do ib=1,Nb
          do k=0,L
             Ker(k)=G(ib)%mats(L-k)*G(ib)%mats(k)
          end do
          ekin_(ib) = -2.d0*params%dtau*kb_trapz(Ker(0:),0,L)
       enddo
    else
       do ib=1,Nb
          do k=0,L
             Ker(k)=G(ib)%lmix(N,k)*conjg(G(ib)%lmix(N,L-k))
          end do
          ekin_(ib)=2.d0*params%dtau*dimag( (-xi)*kb_trapz(Ker(0:),0,L) )
          do k=1,N
             Ker(k)=G(ib)%ret(N,k)*G(ib)%less(k,N)
          end do
          ekin_(ib)=ekin_(ib) + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
          do k=1,N
             Ker(k)=G(ib)%less(N,k)*conjg(G(ib)%ret(N,k))
          end do
          ekin_(ib)=ekin_(ib) + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
       enddo
    endif
    ekin=sum(ekin_)
    deallocate(Ker)
  end function measure_ekin



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! U(t)= U*docc(t) - n(t) + 1/4
  !+-------------------------------------------------------------------+
  function measure_epot(g,sedge,params) result(epot)
    type(kb_contour_gf)     :: g(:)
    type(kb_contour_dgf)     :: sedge(2,size(g))
    type(kb_contour_params)             :: params
    real(8)                             :: epot,docc(size(g),size(g)),nt(size(g))
    integer                             :: i,k,j,N,L,ib,Nb
    N = params%Nt
    L = params%Ntau
    !
    if(N==1)then
       nt   = measure_dens(g,params)
       docc = measure_docc(g,sedge,params)
       epot = Ui*(docc(1,1)+docc(2,2) - (nt(1)+nt(2)) + 0.5d0) + &
            Usti*(2.d0*docc(2,1)+2.d0*docc(1,2) - 2.d0*(nt(1)+nt(2)) + 0.5d0)
    else
       nt   = measure_dens(g,params)
       docc = measure_docc(g,sedge,params)
       epot = U*(docc(1,1)+docc(2,2) - (nt(1)+nt(2)) + 0.5d0) + &
            Ust*(2.d0*docc(2,1)+2.d0*docc(1,2) - 2.d0*(nt(1)+nt(2)) + 0.5d0)
    endif
  end function measure_epot





  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the 1-b hamiltonian H(k) [2d-square] 
  !+-------------------------------------------------------------------+
  function hamk(ik,i,j,params) result(hk)
    integer                 :: ik,i,j
    type(kb_contour_params) :: params
    integer                 :: ix,iy
    type(vect2D)            :: kt
    complex(8)              :: hk
    ix  = ik2ix(ik)
    iy  = ik2iy(ik)
    kt  = kgrid(ix,iy)
    if(i==j)then
       hk  = -2.d0*one*ts*(cos(kt%x)+cos(kt%y))
    else
       hk  = -4.d0*one*tpd*(sin(kt%x)*sin(kt%y))
    endif
  end function Hamk



end PROGRAM neqDMFT


