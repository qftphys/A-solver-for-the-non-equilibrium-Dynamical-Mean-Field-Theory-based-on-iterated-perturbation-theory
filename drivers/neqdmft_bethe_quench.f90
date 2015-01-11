!PURPOSE  : Solve quantum quench dynamic in the 1-band Bethe lattice 
!THIS VERSION IS USED FOR TEST REFERENCE
program neqDMFT
  USE NEQ_DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                             :: i,j,ik,itime,iloop,ix,iy
  logical                             :: converged
  real(8)                             :: wband
  integer                             :: Lk
  character(len=16)                   :: finput
  type(kb_contour_gf)                 :: Sbath
  type(kb_contour_gf)                 :: Gwf
  type(kb_contour_gf)                 :: Sigma
  type(kb_contour_gf)                 :: Gloc
  type(kb_contour_gf)                 :: Ker
  type(kb_contour_dgf)                :: dGwf,dGwf_old
  complex(8),dimension(:),allocatable :: Ham
  !RESULTS:
  real(8),dimension(:),allocatable    :: Hk,Wtk
  real(8),dimension(:,:),allocatable  :: nk

  !READ THE INPUT FILE (in vars_global):
  !=====================================================================
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.in')
  call parse_input_variable(wband,"wband",finput,default=1d0,comment="half-bandwidth W=2t")
  call parse_input_variable(Lk,"Lk",finput,default=100,comment="Number of energy levels for Bethe DOS")
  call read_input_init(trim(finput))


  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  !=====================================================================
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,Niw)
  call setup_kb_contour_params(cc_params,dt,beta)


  !SET THE THERMOSTAT FUNCTION (in neq_thermostat):
  !=====================================================================
  call allocate_kb_contour_gf(Sbath,cc_params)
  call get_thermostat_bath(Sbath,cc_params)


  !ALLOCATE ALL THE FUNCTIONS INVOLVED IN THE CALCULATION:
  !=====================================================================
  call allocate_kb_contour_gf(Sigma,cc_params) !Self-Energy function
  call allocate_kb_contour_gf(Gloc,cc_params)  !Local Green's function
  call allocate_kb_contour_gf(Gwf,cc_params)   !Local Weiss-Field function
  call allocate_kb_contour_dgf(dGwf,cc_params) !d/dt Local Weiss-Field
  call allocate_kb_contour_dgf(dGwf_old,cc_params)
  call allocate_kb_contour_gf(Ker,cc_params)
  allocate(ham(cc_params%Ntime))


  !GET THE STARTING CONDITIONS: 
  !build the equilibrium t=0 functions
  !=====================================================================
  cc_params%Itime=1
  call neq_continue_equilibirum(Gwf,dGwf,Gloc,Sigma,cc_params,wband)
  call neq_measure_observables(Gloc,Sigma,cc_params)

  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  !=====================================================================
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Itime=itime
     !
     call neq_setup_weiss_field(Gwf,cc_params)
     dGwf_old = dGwf
     !
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        write(*,"(A,I2,A1)",advance='no')"dmft loop=",iloop," "
        !
        !IMPURITY SOLVER: IPT
        call neq_solve_ipt(Gwf,Sigma,cc_params)
        !Gimp = G0 + K*G, K=G0*Sigma, Q = G0
        call convolute_kb_contour_gf(Gwf,Sigma,Ker,cc_params)
        call vie_kb_contour_gf(Gwf,Ker,Gloc,cc_params)
        !
        !PERFORM THE SELF_CONSISTENCY 
        !update Weiss-Field: id/dt Gwf(t,t')-Gloc*Gwf(t,t')=delta(t,t')
        !this is the analog of \Delta = t^2Gloc
        Ham = zero
        call vide_kb_contour_gf(Ham,Gloc,Gwf,dGwf_old,dGwf,cc_params)
        !
        !CHECK CONVERGENCE
        converged = convergence_check(Gwf,cc_params)
        !
     enddo
     call plot_kb_contour_gf("Sigma",Sigma,cc_params)
     call plot_kb_contour_gf("Gloc",Gloc,cc_params)
     call plot_kb_contour_gf("G0",Gwf,cc_params)

     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
     call neq_measure_observables(Gloc,Sigma,cc_params)
  enddo


  print*,"Getting n(e,t):"
  allocate(Hk(Lk),Wtk(Lk))
  allocate(nk(cc_params%Ntime,Lk))
  call bethe_lattice(Wtk,Hk,Lk,Wband)
  !Gwf,dGwf <--- G_e(t,t'),d/dtG_e(t,t')
  do ik=1,Lk
     ham = Hk(ik)
     call neq_setup_initial_conditions(Gwf,dGwf_old,Sigma,Hk(ik),cc_params)
     do itime=2,cc_params%Ntime
        cc_params%Itime = itime
        call vide_kb_contour_gf(Ham,Sigma,Gwf,dGwf_old,dGwf,cc_params)
        dGwf_old = dGwf
     enddo
     forall(itime=1:cc_params%Ntime)nk(itime,ik)=dimag(Gwf%less(itime,itime))
  enddo
  call splot3d("nkVSepsikVStime.plot",cc_params%t,Hk,nk,wlines=.true.)


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check(G,params) result(converged)
    type(kb_contour_gf)                 :: G
    type(kb_contour_params)             :: params
    logical                             :: converged
    complex(8),dimension(:),allocatable :: test_func
    integer :: i,N,L,Ntot
    !
    N   = params%Itime                 !<== work with the ACTUAL size of the contour
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



! subroutine test_(G,Self,params)
!   type(kb_contour_gf)                   :: G,Self,GxS
!   type(kb_contour_params)               :: params
!   integer                               :: N,L,i,j,Ntot
!   complex(8),dimension(:,:),allocatable :: Gmat,Smat,GxSmat
!   real(8),dimension(:),allocatable      :: time
!   N=params%Ntime
!   L=params%Ntau
!   Ntot=2*N+L+1
!   allocate(time(Ntot))
!   allocate(Gmat(Ntot,Ntot),Smat(Ntot,Ntot),GxSmat(Ntot,Ntot))
!   time(1:N)=params%t(1:N)
!   time(N+1:2*N)=time(N) + params%t(1:N)
!   time(2*N+1:2*N+L+1)=time(2*N)+params%tau(0:L)/beta
!   call kb_contour_gf2kb_matrix(G,N,L,Gmat)
!   call splot3d("Gmat_test.plot",time,time,Gmat)
!   call splot3d("Gmat_13.plot",params%t(1:N),params%tau(0:L),Gmat(1:N,2*N+1:2*N+1+L))
!   call splot3d("Gmat_12.plot",params%t(1:N),params%t(1:N),Gmat(1:N,N+1:2*N))
!   G=zero
!   call kb_matrix2kb_contour_gf(Gmat,N,L,G)
!   call plot_kb_contour_gf("G_test",G,params)
!   !
!   !get convolution in ordinary way:
!   call allocate_kb_contour_gf(GxS,params)
!   do i=1,params%Ntime
!      params%Nt=i
!      call convolute_kb_contour_gf(G,Self,GxS,params)
!   enddo
!   call plot_kb_contour_gf("GxS",GxS,params)
!   !
!   !map convolution to a matrix:
!   call kb_contour_gf2kb_matrix(GxS,N,L,Smat)
!   call splot3d("GxSmat_test1.plot",time,time,Smat)
!   !
!   !Evaluate convolution using Keldysh Matrices
!   call kb_contour_gf2kb_matrix(Self,N,L,Smat)
!   call convolute_kb_matrix_gf(Gmat,Smat,N,L,params,GxSmat)
!   call splot3d("GxSmat_test2.plot",time,time,GxSmat)
!   call splot3d("GxS_13.plot",params%t(1:N),params%tau(0:L),GxSmat(1:N,2*N+1:2*N+1+L))
!   call splot3d("GxS_31.plot",params%tau(0:L),params%t(1:N),GxSmat(2*N+1:2*N+1+L,1:N))
!   GxS=zero
!   call kb_matrix2kb_contour_gf(GxSmat,N,L,GxS)
!   call plot_kb_contour_gf("GxS_test",GxS,params)
! end subroutine test_




! !!Sigma^(4c):
! print*,"Get 4c"
! Smat=zero
! do i=1,Ntot
!    do j=1,Ntot
!       int_ij=zero

!       inta=zero
!       do ia=1,N
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia)*G0mat(ia,ib)*G0mat(ib,j)*G2mat(i,ib)*G2mat(ia,j)
!          enddo
!          inta(ia) = inta(ia) + params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia)*G0mat(ia,ib+N)*G0mat(ib+N,j)*G2mat(i,ib+N)*G2mat(ia,j)
!          enddo
!          inta(ia) = inta(ia) - params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=0,L
!             intb(ib) = G0mat(i,ia)*G0mat(ia,2*N+1+ib)*G0mat(2*N+1+ib,j)*G2mat(i,2*N+1+ib)*G2mat(ia,j)
!          enddo
!          inta(ia) = inta(ia) -xi*params%dtau*kb_trapz(intb(0:),0,L)
!       enddo
!       int_ij = int_ij + params%dt*kb_trapz(inta(0:),1,N)

!       inta=zero
!       do ia=1,N
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia+N)*G0mat(ia+N,ib)*G0mat(ib,j)*G2mat(i,ib)*G2mat(ia+N,j)
!          enddo
!          inta(ia) = inta(ia) + params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia+N)*G0mat(ia+N,ib+N)*G0mat(ib+N,j)*G2mat(i,ib+N)*G2mat(ia+N,j)
!          enddo
!          inta(ia) = inta(ia) - params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=0,L
!             intb(ib) = G0mat(i,ia+N)*G0mat(ia+N,2*N+1+ib)*G0mat(2*N+1+ib,j)*G2mat(i,2*N+1+ib)*G2mat(ia+N,j)
!          enddo
!          inta(ia) = inta(ia) -xi*params%dtau*kb_trapz(intb(0:),0,L)
!       enddo
!       int_ij = int_ij - params%dt*kb_trapz(inta(0:),1,N)

!       inta=zero
!       do ia=0,L
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,2*N+1+ia)*G0mat(2*N+1+ia,ib)*G0mat(ib,j)*G2mat(i,ib)*G2mat(2*N+1+ia,j)
!          enddo
!          inta(ia) = inta(ia) + params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,2*N+1+ia)*G0mat(2*N+1+ia,ib+N)*G0mat(ib+N,j)*G2mat(i,ib+N)*G2mat(2*N+1+ia,j)
!          enddo
!          inta(ia) = inta(ia) - params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=0,L
!             intb(ib) = G0mat(i,2*N+1+ia)*G0mat(2*N+1+ia,2*N+1+ib)*G0mat(2*N+1+ib,j)*G2mat(i,2*N+1+ib)*G2mat(2*N+1+ia,j)
!          enddo
!          inta(ia) = inta(ia) -xi*params%dtau*kb_trapz(intb(0:),0,L)
!       enddo
!       int_ij = int_ij -xi*params%dtau*kb_trapz(inta(0:),0,L)

!       !
!       Smat(i,j) = int_ij
!    enddo
! enddo
! print*,"extract Sigma^c:"
! call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(3))



! !!Sigma^(4d):
! print*,"Get 4d"
! Smat=zero
! do i=1,Ntot
!    do j=1,Ntot
!       int_ij=zero

!       inta=zero
!       do ia=1,N
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia)*G0mat(ia,j)*G0mat(i,ib)*G0mat(ib,j)*G2mat(ia,ib)
!          enddo
!          inta(ia) = inta(ia) + params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia)*G0mat(ia,j)*G0mat(i,N+ib)*G0mat(N+ib,j)*G2mat(ia,N+ib)
!          enddo
!          inta(ia) = inta(ia) - params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=0,L
!             intb(ib) = G0mat(i,ia)*G0mat(ia,j)*G0mat(i,2*N+i+ib)*G0mat(2*N+i+ib,j)*G2mat(ia,2*N+i+ib)
!          enddo
!          inta(ia) = inta(ia) -xi*params%dtau*kb_trapz(intb(0:),0,L)
!       enddo
!       int_ij = int_ij + params%dt*kb_trapz(inta(0:),1,N)

!       inta=zero
!       do ia=1,N
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia+N)*G0mat(ia+N,j)*G0mat(i,ib)*G0mat(ib,j)*G2mat(ia+N,ib)
!          enddo
!          inta(ia) = inta(ia) + params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,ia+N)*G0mat(ia+N,j)*G0mat(i,N+ib)*G0mat(N+ib,j)*G2mat(ia+N,N+ib)
!          enddo
!          inta(ia) = inta(ia) - params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=0,L
!             intb(ib) = G0mat(i,ia+N)*G0mat(ia+N,j)*G0mat(i,2*N+i+ib)*G0mat(2*N+i+ib,j)*G2mat(ia+N,2*N+i+ib)
!          enddo
!          inta(ia) = inta(ia) -xi*params%dtau*kb_trapz(intb(0:),0,L)
!       enddo
!       int_ij = int_ij - params%dt*kb_trapz(inta(0:),1,N)


!       inta=zero
!       do ia=0,L
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,2*N+1+ia)*G0mat(2*N+1+ia,j)*G0mat(i,ib)*G0mat(ib,j)*G2mat(2*N+1+ia,ib)
!          enddo
!          inta(ia) = inta(ia) + params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=1,N
!             intb(ib) = G0mat(i,2*N+1+ia)*G0mat(2*N+1+ia,j)*G0mat(i,N+ib)*G0mat(N+ib,j)*G2mat(2*N+1+ia,N+ib)
!          enddo
!          inta(ia) = inta(ia) - params%dt*kb_trapz(intb(0:),1,N)
!          intb(0:)=zero
!          do ib=0,L
!             intb(ib) = G0mat(i,2*N+1+ia)*G0mat(2*N+1+ia,j)*G0mat(i,2*N+i+ib)*G0mat(2*N+i+ib,j)*G2mat(2*N+1+ia,2*N+i+ib)
!          enddo
!          inta(ia) = inta(ia) -xi*params%dtau*kb_trapz(intb(0:),0,L)
!       enddo
!       int_ij = int_ij - xi*params%dtau*kb_trapz(inta(0:),0,L)


!       !
!       Smat(i,j) = G0mat(i,j)*int_ij
!    enddo
! enddo
! print*,"extract Sigma^c:"
! call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(4))

