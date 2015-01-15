program neqDMFT
  USE NEQ_DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                               :: i,j,ik,itime,iloop,ix,iy
  logical                               :: converged,t_converged
  character(len=16)                     :: finput
  real(8)                               :: wband
  integer                               :: Lk
  type(kb_contour_gf)                   :: Sbath
  type(kb_contour_gf)                   :: Gloc,Sigma,Gwf
  type(kb_contour_gf)                   :: Ker
  type(kb_contour_gf),allocatable       :: Gk(:)
  type(kb_contour_dgf),allocatable      :: dGk(:),dGk_old(:)
  type(kb_contour_dgf)                  :: Gedge,G0edge
  complex(8),dimension(:,:),allocatable :: Hk
  real(8),dimension(:),allocatable      :: Epsk,Wtk
  !RESULTS:
  real(8),dimension(:,:),allocatable    :: nk
  real(8)                               :: de

  !READ THE INPUT FILE (in vars_global):
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.conf')
  call parse_input_variable(wband,"wband",finput,default=1d0,comment="half-bandwidth W=2t")
  call parse_input_variable(Lk,"Lk",finput,default=100,comment="Number of energy levels for Bethe DOS")
  call read_input_init(trim(finput))


  !BUILD THE LATTICE STRUCTURE (in lib/square_lattice):
  allocate(Hk(Ntime,Lk),Epsk(Lk),Wtk(Lk))
  Epsk = linspace(-wband,wband,Lk,mesh=de)
  do i=1,Lk
     Wtk(i) = dens_bethe(Epsk(i),wband)
  enddo
  call splot("DOSbethe.ipt",Epsk,Wtk)
  Wtk=Wtk/sum(Wtk)


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
  itime=1
  cc_params%Itime=itime
  Gloc = zero
  call neq_continue_equilibirum(Gwf,Gk,dGk,Gloc,Sigma,Epsk,Wtk,cc_params)
  call neq_measure_observables(Gloc,Sigma,cc_params)
  do ik=1,Lk
     nk(1,ik)=dimag(Gk(ik)%less(1,1))
  enddo
  Hk=zero
  do i=1,cc_params%Ntime
     Hk(i,:)=Epsk
  enddo


  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  ! t_converged=.false.
  ! do while(.not.t_converged.AND.itime<cc_params%Ntime)
  !    itime=itime+1
  do itime=2,cc_params%Ntime
     cc_params%Itime=itime
     print*,""
     print*,"time step=",itime
     !prepare the weiss-field at this actual time_step for DMFT:
     call neq_setup_weiss_field(Gwf,cc_params)
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
        Gedge=zero
        do ik=1,Lk
           call vide_kb_contour_gf(Hk(:,ik),Sigma,Gk(ik),dGk_old(ik),dGk(ik),cc_params)
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


     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
     call neq_measure_observables(Gloc,Sigma,cc_params)
     forall(ik=1:Lk)nk(itime,ik)=dimag(Gk(ik)%less(itime,itime))
  enddo


  !EVALUATE AND PRINT OTHER RESULTS OF THE CALCULATION
  call splot3d("nkVSepsikVStime.ipt",cc_params%t,dreal(Hk(1,:)),nk)
  call plot_kb_contour_gf("Sigma.ipt",Sigma,cc_params)
  call plot_kb_contour_gf("Gloc.ipt",Gloc,cc_params)
  call plot_kb_contour_gf("G0.ipt",Gwf,cc_params)
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
  end function convergence_check


end PROGRAM neqDMFT
