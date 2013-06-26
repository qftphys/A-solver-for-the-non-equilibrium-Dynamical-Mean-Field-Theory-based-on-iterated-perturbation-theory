!###################################################################
!PURPOSE  : Solve conduction band electrons driven 
! by electric field interacting with a resevoir of free 
! electrons at temperature T
!AUTHORS  : Adriano Amaricc & G.Mazza & C.Weber
!###################################################################
program neqDMFT
  USE NEQ_VARS_GLOBAL   !global variables, calls to 3rd library 
  USE ELECTRIC_FIELD    !contains electric field init && routines
  USE BATH              !contains bath inizialization
  USE EQ_IPT            !solves the equilibrium problem w/ IPT
  USE NEQ_IPT           !performs the non-eq. IPT. Write Sigma
  USE NEQ_UPDATE_WF     !contains routines for WF update and printing.
  USE NEQ_KADANOFF_BAYM !solves KB equations numerically to get k-sum
  USE RESULTS           !get results and print observables
  implicit none

  integer :: i
  logical :: converged

  !READ THE INPUT FILE (in vars_global):
  call read_input_init("inputFILE.in")

  !BUILD THE TIME,FREQUENCY GRIDS:
  allocate(t(1:nstep))
  tmax = dt*real(nstep-1,8)
  t    = linspace(0.d0,tmax,Nstep)
  Nfit=2*Nstep-1
  dtfit= dt/2.d0
  allocate(tfit(Nfit))
  tfit = linspace(0.d0,tmax,Nfit)
  write(*,'(A,2F12.6)')"dt   =",dt,dtfit
  write(*,'(A,F12.6)')"tmax =",tmax

  !BUILD THE 2D-SQUARE LATTICE STRUCTURE (in lib/square_lattice):
  Lk   = square_lattice_dimension(Nx,Ny)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Ny)
  epsik= square_lattice_dispersion_array(Lk,ts)
  call get_free_dos(epsik,wt)

  !SET THE ELECTRIC FIELD (in electric_field):
  call set_efield_vector()

  !ALLOCATE FUNCTIONS IN THE MEMORY (in vars_global):
  call global_memory_allocation

  !BUILD THE  DISSIPATIVE BATH FUNCTIONS (in bath):
  call get_thermostat_bath()

  !SOLVE THE EQUILIBRIUM PROBLEM WITH IPT (in equilibrium):
  if(solve_eq)call solve_equilibrium_ipt()

  !START DMFT LOOP SEQUENCE:
  !==============================================================
  !initialize the run by guessing/reading the self-energy functions (in IPT_NEQ.f90):
  call neq_init_run

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
             Ak= Afield(t(i),Ek)
             Jk= nk(i,ik)*square_lattice_velocity(kgrid(ix,iy) - Ak)
             Jloc(i) = Jloc(i) +  wt(ik)*Jk
          enddo
       enddo
       test_func(:)=modulo(Jloc(:))
    else
       do i=1,nstep
          test_func(i)=-xi*pack_less_tri(i,i,locG)
       enddo
    endif
    converged=check_convergence(test_func(:),eps_error,Nsuccess,nloop,oerr=err)
    if(isnan(err))call error("Aborted convergence: error=NaN")
  end function convergence_check



end PROGRAM neqDMFT
