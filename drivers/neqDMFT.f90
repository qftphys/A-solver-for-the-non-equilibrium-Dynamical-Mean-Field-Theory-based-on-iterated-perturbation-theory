!###################################################################
!PURPOSE  : Solve conduction band electrons driven 
! by electric field interacting with a resevoir of free 
! electrons at temperature T
!AUTHORS  : Adriano Amaricci && Cedric Weber
!###################################################################
program neqDMFT
  USE NEQ_VARS_GLOBAL   !global variables, calls to 3rd library 
  USE ELECTRIC_FIELD    !contains electric field init && routines
  USE NEQ_BATH          !contains bath inizialization
  USE NEQ_IPT           !performs the non-eq. IPT. Write Sigma
  USE NEQ_UPDATE_WF     !contains routines for WF update and printing.
  USE NEQ_KADANOFF_BAYM !solves KB equations numerically to get k-sum
  USE NEQ_RESULTS
  implicit none

  integer :: i
  logical :: converged
  real(8),allocatable,dimension(:) :: trel,wrel
  real(8) :: dtrel,dwrel
  integer :: Nrel

  !READ THE INPUT FILE (in vars_global):
  call read_input_init("inputFILE.in")

  !BUILD THE TIME,FREQUENCY GRIDS:
  allocate(time(Nstep))
  tmax = dt*real(Nstep-1,8)
  time = linspace(0.d0,tmax,nstep)

  Nrel=Nstep-1
  allocate(trel(-Nrel:Nrel))
  allocate(wrel(2*Nrel))
  wmax=pi/dt
  trel=linspace(-dt*real(Nrel,8),dt*real(Nrel,8),2*Nrel+1,mesh=dtrel)
  wrel=linspace(-wmax,wmax,2*Nrel,mesh=dwrel)
  ! allocate(wr(2*nstep),t(-nstep:nstep))
  ! allocate(wm(L),tau(0:Ltau))
  ! wmax  = pi/dt
  ! t     = linspace(-dt*real(nstep,8),dt*real(nstep,8),2*nstep+1)
  ! wr    = linspace(-wmax,wmax,2*nstep,mesh=fmesh)-fmesh/2.d0
  ! !t     = wr(-nstep:nstep)/fmesh*dt
  ! wm    = pi/beta*real(2*arange(1,L)-1,8)
  ! tau   = linspace(0.d0,beta,Ltau+1,mesh=dtau)

  write(*,'(A,F12.6)')"dt   =",dt
  write(*,'(A,F12.6)')"dtrel=",dtrel
  write(*,'(A,F12.6)')"dwrel=",dwrel
  write(*,'(A,F12.6)')"tmax =",tmax


  !BUILD THE 2D-SQUARE LATTICE STRUCTURE (in lib/square_lattice):
  Lk   = square_lattice_dimension(Nx,Nx)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)
  call get_free_dos(epsik,wt)

  !SET THE ELECTRIC FIELD (in electric_field):
  call set_efield_vector()
  call print_field(time)

  !ALLOCATE FUNCTIONS IN THE MEMORY (in vars_global):
  call global_memory_allocation

  !BUILD THE  DISSIPATIVE BATH FUNCTIONS (in bath):
  call get_thermostat_bath()

  ! !SOLVE THE EQUILIBRIUM PROBLEM WITH IPT (in equilibrium):
  ! if(solve_eq)call solve_equilibrium_ipt()

  !START DMFT LOOP SEQUENCE:
  !==============================================================
  !initialize the run by guessing/reading the self-energy functions (in IPT_NEQ.f90):
  call read_initial_self

  iloop=0;converged=.false.
  do while(.not.converged);iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop",unit=6)
     !
     call neq_get_localgf        !-|(in kadanoff-baym)

     call neq_update_weiss_field !-|SELF-CONSISTENCY (in funx_neq)
     !
     call neq_solve_ipt          !-|IMPURITY SOLVER (in ipt_neq)
     !
     converged = convergence_check()
     !
     call end_loop()
  enddo

  !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
  call plot_results

  call msg("BRAVO")


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check() result(converged)
    logical                       :: converged
    integer                       :: i,ik,ix,iy
    type(vect2D)                  :: Jk,Ak
    type(vect2D),dimension(nstep) :: Jloc                   !local Current 
    real(8),dimension(nstep)      :: test_func
    integer                       :: selector
    real(8)                       :: err
    if(Efield/=0.d0)then
       Jloc=Vzero
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)
          do i=1,nstep
             Ak= Afield(time(i),Ek)
             Jk= nk(i,ik)*square_lattice_velocity(kgrid(ix,iy) - Ak)
             Jloc(i) = Jloc(i) +  wt(ik)*Jk
          enddo
       enddo
       test_func=modulo(Jloc)
    else
       forall(i=1:nstep)test_func(i)=-xi*locG%less(i,i)
    endif
    converged=check_convergence(test_func,eps_error,Nsuccess,nloop,oerr=err)
    if(isnan(err))call error("Aborted convergence: error=NaN")
  end function convergence_check



end PROGRAM neqDMFT
