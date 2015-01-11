program neqDMFT
  USE NEQ_DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                   :: i,j,ik,itime,iloop,ix,iy
  logical                                   :: converged
  character(len=16)                         :: finput
  real(8)                                   :: ts
  integer                                   :: Nx,Ny,Lk
  type(kb_contour_gf)                       :: Sbath
  type(kb_contour_gf)                       :: Gloc,Sigma,Gwf
  type(kb_contour_gf)                       :: Ker
  type(kb_contour_gf),allocatable           :: Gk(:)
  type(kb_contour_dgf),allocatable          :: dGk(:),dGk_old(:)
  type(kb_contour_dgf)                      :: Gedge
  complex(8),dimension(:,:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable          :: kxgrid,kygrid,Wtk
  !RESULTS:
  real(8),dimension(:,:,:),allocatable      :: nDens
  real(8),dimension(:),allocatable          :: n_t
  real(8),dimension(:,:),allocatable        :: nk

  !READ THE INPUT FILE (in vars_global):
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.in')
  call parse_input_variable(ts,"ts",finput,default=1d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=20,comment="Number of points along x-axis of the BZ")
  call parse_input_variable(Ny,"Ny",finput,default=20,comment="Number of points along y-axis of the BZ")
  call read_input_init(trim(finput))



  !BUILD THE LATTICE STRUCTURE (in lib/square_lattice):
  Lk = Nx*Ny
  allocate(Hk(Ntime,Lk,1,1),Wtk(Lk))
  allocate(kxgrid(Nx),kygrid(Ny))
  write(*,*) "Using Nk_total="//txtfy(Lk)
  kxgrid = kgrid(Nx)
  kygrid = kgrid(Ny)
  Hk(1,:,:,:) = build_hk_model(Lk,1,hk_model,kxgrid,kygrid,[0d0])
  Wtk    = 1d0/Lk
  call write_hk_w90("Hk2d.dat",1,2,1,hk_model,kxgrid,kygrid,[0d0])


  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,Niw)
  call setup_kb_contour_params(cc_params,dt,beta)



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
  allocate(nk(cc_params%Ntime,Lk))



  !READ OR GUESS THE INITIAL WEISS FIELD G0 (t=t'=0)
  cc_params%Itime=1
  Gloc = zero
  call neq_continue_equilibirum(Gwf,Gk,dGk,Gloc,Sigma,dreal(Hk(1,:,1,1)),Wtk,cc_params)
  call neq_measure_observables(Gloc,Sigma,cc_params)
  do ik=1,Lk
     nk(1,ik)=dimag(Gk(ik)%less(1,1))
  enddo
  do i=2,cc_params%Ntime
     Hk(i,:,:,:)=Hk(1,:,:,:)
  enddo

  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Itime=itime
     !prepare the weiss-field at this actual time_step for DMFT:
     call neq_setup_weiss_field(Gwf,cc_params)
     do ik=1,Lk
        dGk_old(ik) = dGk(ik)
     enddo
     !
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        write(*,"(A,I4,A1)",advance='no')"dmft loop=",iloop," "

        !IMPURITY SOLVER: IPT.
        call neq_solve_ipt(Gwf,Sigma,cc_params)

        !PERFORM THE SELF_CONSISTENCY: get local GF + update Weiss-Field
        Gedge=zero
        do ik=1,Lk
           call vide_kb_contour_gf(Hk(:,ik,1,1),Sigma,Gk(ik),dGk_old(ik),dGk(ik),cc_params)
           Gedge%ret(1:itime) = Gedge%ret(1:itime)  + Wtk(ik)*Gk(ik)%ret(itime,1:itime)
           Gedge%less(1:itime)= Gedge%less(1:itime) + Wtk(ik)*Gk(ik)%less(itime,1:itime)
           Gedge%lmix(:)      = Gedge%lmix(:)       + Wtk(ik)*Gk(ik)%lmix(itime,:)
        enddo
        Gloc%ret(itime,1:itime)   = Gedge%ret(1:itime)
        Gloc%less(itime,1:itime)  = Gedge%less(1:itime)
        Gloc%lmix(itime,:)        = Gedge%lmix(:)
        Gloc%less(1:itime-1,itime)=-conjg(Gedge%less(1:itime-1))

        !update the weiss field by solving the integral equation:
        ! G0 + K*G0 = Q , with K = G*\Sigma and Q = G
        call convolute_kb_contour_gf(Gloc,Sigma,Ker,cc_params,dcoeff=-1.d0)
        call vie_kb_contour_gf(Gloc,Ker,Gwf,cc_params)
        !
        !CHECK CONVERGENCE
        converged = convergence_check(Gwf,cc_params)
     enddo
     call plot_kb_contour_gf("Sigma",Sigma,cc_params)
     call plot_kb_contour_gf("Gloc",Gloc,cc_params)
     call plot_kb_contour_gf("G0",Gwf,cc_params)

     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
     call neq_measure_observables(Gloc,Sigma,cc_params)
     forall(ik=1:Lk)nk(itime,ik)=dimag(Gk(ik)%less(itime,itime))
  enddo


  !EVALUATE AND PRINT OTHER RESULTS OF THE CALCULATION
  allocate(ndens(1:Nx,1:Ny,cc_params%Ntime))
  do i=1,cc_params%Ntime
     do ik=1,Lk
        ix = indx2ix(ik,[Nx,Ny,1])
        iy = indx2iy(ik,[Nx,Ny,1])
        nDens(ix,iy,i)=nk(i,ik)
     enddo
  enddo
  call splot3d("3dFSVSpiVSt.plot",kxgrid,kygrid,nDens(1:Nx,1:Nx,:))
  call splot3d("nkVSepsikVStime.plot",cc_params%t,dreal(Hk(1,:,1,1)),nk)

  print*,"BRAVO"

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the 1-b hamiltonian H(k) [2d-square] 
  !+-------------------------------------------------------------------+
  ! function hamk(ik,params) result(hk)
  !   integer                 :: ik
  !   type(kb_contour_params) :: params
  !   integer                 :: ix,iy
  !   type(vect2D)            :: kt
  !   complex(8)              :: hk
  !   ix  = ik2ix(ik)
  !   iy  = ik2iy(ik)
  !   hk  = one*square_lattice_dispersion(kgrid(ix,iy))
  ! end function Hamk
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky
    real(8)                   :: ck
    complex(8),dimension(N,N) :: hk
    kx=kpoint(1)
    ky=kpoint(2)
    !
    ck=cos(kx)+cos(ky)
    !
    Hk      = zero
    Hk(1,1) = -2d0*ts*ck
  end function hk_model













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
    N   = params%itime                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    Ntot=2*N+L
    allocate(test_func(Ntot))
    test_func=zero
    do i=1,N
       test_func(i)  = G%ret(N,i)
       test_func(N+i)= G%less(N,i)
    enddo
    do i=1,L
       test_func(2*N+i)=G%lmix(N,i)
    enddo
    converged=check_convergence(test_func,dmft_error,Nsuccess,nloop)
    deallocate(test_func)
    !if(isnan(err))stop "Aborted convergence: error=NaN"
  end function convergence_check






end PROGRAM neqDMFT
