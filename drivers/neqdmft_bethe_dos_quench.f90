program neqDMFT
  USE NEQ_DMFT_IPT
  ! USE NEQ_CONTOUR
  ! USE NEQ_CONTOUR_GF
  ! USE NEQ_VARS_GLOBAL
  ! USE NEQ_AUX_FUNX
  ! USE NEQ_MEASURE
  ! USE NEQ_IPT
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


  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  !=====================================================================
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,Niw)
  call setup_kb_contour_params(cc_params,dt,beta)


  !BUILD THE LATTICE STRUCTURE (in lib/square_lattice):
  allocate(epsik(Lk),wt(Lk))
  epsik = linspace(-wband,wband,Lk,mesh=de)
  do i=1,Lk
     wt(i) = dens_bethe(epsik(i),wband)
  enddo
  call splot("DOSbethe.ipt",epsik,wt)
  wt=wt/sum(wt)



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
     dGk(ik)=zero
     dGk_old(ik)=zero
  end do
  allocate(ham(cc_params%Ntime,Lk))
  allocate(nk(cc_params%Ntime,Lk))



  !READ OR GUESS THE INITIAL WEISS FIELD G0 (t=t'=0)
  cc_params%Nt=1
  Gloc = zero
  call neq_continue_equilibirum(Gwf,Gk,dGk,Gloc,Sigma,epsik,wt,cc_params)
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
     call  extrapolate_kb_contour_gf(Gwf,cc_params)
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
        call del_kb_contour_gf(Gloc,cc_params)
        do ik=1,Lk
           call vide_kb_contour_gf(Ham(:,ik),Sigma,Gk(ik),dGk_old(ik),dGk(ik),cc_params)
        enddo
        call sum_kb_contour_gf(Gk(:),wt(:),Gloc,cc_params)

        !update the weiss field by solving the integral equation:
        ! G0  = Q + K*G0 , with K = -G*\Sigma and Q = G
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
  call splot3d("nkVSepsikVStime.nipt",cc_params%t,epsik,nk)
  call plot_kb_contour_gf("Sigma.nipt",Sigma,cc_params)
  call plot_kb_contour_gf("Gloc.nipt",Gloc,cc_params)
  call plot_kb_contour_gf("G0.nipt",Gwf,cc_params)
  print*,"BRAVO"


contains


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



end PROGRAM neqDMFT
