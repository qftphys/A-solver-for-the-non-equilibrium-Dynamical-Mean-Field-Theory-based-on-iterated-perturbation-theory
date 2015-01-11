!###################################################################
!PURPOSE  : Solve quantum quench dynamic in the 1-band Bethe lattice 
!AUTHORS  : Adriano Amaricci 
!###################################################################
program neqDMFT
  USE ERROR
  USE FFTGF
  USE NEQ_VARS_GLOBAL   !global variables, calls to 3rd library 
  USE NEQ_THERMOSTAT    !contains bath inizialization
  implicit none
  integer                             :: i,j,ib,ik,itime,iloop,ix,iy
  logical                             :: converged
  character(len=16)                   :: finput
  type(kb_contour_gf),allocatable     :: Gwf(:)
  type(kb_contour_gf),allocatable     :: Sigma(:)
  type(kb_contour_gf),allocatable     :: Gloc(:)
  type(kb_contour_gf)                 :: Ker
  type(kb_contour_dgf),allocatable    :: dGwf(:),dGwf_old(:),Sedge(:,:)
  complex(8),dimension(:),allocatable :: Ham
  !RESULTS:
  real(8),dimension(:),allocatable    :: dos,epsi
  real(8),dimension(:,:),allocatable  :: nk

  !READ THE INPUT FILE (in vars_global):
  !=====================================================================
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.in')
  call read_input_init(trim(finput))
  ts=1.d0

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
  allocate(Sigma(2),Gloc(2),Gwf(2),dGwf(2),dGwf_old(2),Sedge(2,2))
  do ib=1,2
     call allocate_kb_contour_gf(Sigma(ib),cc_params) !Self-Energy function
     call allocate_kb_contour_gf(Gloc(ib),cc_params)  !Local Green's function
     call allocate_kb_contour_gf(Gwf(ib),cc_params)   !Local Weiss-Field function
     call allocate_kb_contour_dgf(dGwf(ib),cc_params) !d/dt Local Weiss-Field
     call allocate_kb_contour_dgf(dGwf_old(ib),cc_params)
     call allocate_kb_contour_dgf(Sedge(ib,1),cc_params,wgtr=.true.)
     call allocate_kb_contour_dgf(Sedge(ib,2),cc_params,wgtr=.true.)
  enddo
  call allocate_kb_contour_gf(Ker,cc_params)
  allocate(ham(cc_params%Ntime))


  !GET THE STARTING CONDITIONS: 
  !build the equilibrium t=0 functions
  !=====================================================================
  cc_params%Nt=1
  call init_equilibrium_functions(Gwf,dGwf,Gloc,Sigma,cc_params)
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
     do ib=1,2
        dGwf_old(ib) = dGwf(ib)
     enddo
     !
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        write(*,"(A,I2,A1)",advance='no')"dmft loop=",iloop," "
        !
        !IMPURITY SOLVER: IPT
        call neq_solve_ipt(Gwf,Sigma,cc_params)
        !get the impurity green's function:
        !G = G0 + K*G, K=G0*Sigma, Q = G0
        do ib=1,2
           call convolute_kb_contour_gf(Gwf(ib),Sigma(ib),Ker,cc_params)
           call vie_kb_contour_gf(Gwf(ib),Ker,Gloc(ib),cc_params)
        enddo
        !
        !PERFORM THE SELF_CONSISTENCY 
        !update Weiss-Field: id/dt Gwf(t,t')-Gloc*Gwf(t,t')=delta(t,t')
        !this is the analog of \Delta = t^2Gloc
        Ham = zero
        do ib=1,2
           call vide_kb_contour_gf(Ham,Gloc(ib),Gwf(ib),dGwf_old(ib),dGwf(ib),cc_params)
        enddo
        !
        !CHECK CONVERGENCE
        converged = convergence_check(Gwf,cc_params)
        !
     enddo
     do ib=1,2
        call plot_kb_contour_gf("Sigma_"//txtfy(ib),Sigma(ib),cc_params)
        call plot_kb_contour_gf("Gloc_"//txtfy(ib),Gloc(ib),cc_params)
        call plot_kb_contour_gf("G0_"//txtfy(ib),Gwf(ib),cc_params)
     enddo
     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION

     call measure_observables(Gloc,Sigma,Sedge,cc_params)
  enddo

  print*,"BRAVO"

contains



  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine init_equilibrium_functions(g0,dg0,g,self,params)
    type(kb_contour_gf)                 :: g0(:)
    type(kb_contour_dgf)                :: dg0(:)
    type(kb_contour_gf)                 :: g(:)
    type(kb_contour_gf)                 :: self(:)
    type(kb_contour_params)             :: params
    integer                             :: Nb
    real(8)                             :: wm,res(2),ims(2)
    logical                             :: bool
    integer                             :: i,j,ib,k,ik,unit,len,N,L,Lf
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: GxG0
    if(.not.g0(1)%status)stop "init_functions: g0 is not allocated"
    if(.not.dg0(1)%status)stop "init_functions: dg0 is not allocated"
    if(.not.g(1)%status)stop "init_functions: g is not allocated"
    if(.not.self(1)%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
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
          zeta  = xi*wm
          g0(1)%iw(i) = gfbethe(wm,zeta,2.d0*ts)
          g0(2)%iw(i) = gfbethe(wm,zeta,2.d0*ts)
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
    do ib=1,Nb
       do i=1,Lf
          wm    = pi/beta*dble(2*i-1)
          zeta  = xi*wm - Self(ib)%iw(i)
          G(ib)%iw(i) = gfbethe(wm,zeta,2.d0*ts)
       enddo
       call fftgf_iw2tau(G(ib)%iw,G(ib)%mats(0:),beta)         !get G(tau)
       G(ib)%ret(1,1)  = -xi                               !get G^R(0,0)=-xi
       G(ib)%less(1,1) = -xi*G(ib)%mats(L)                  !get G^<(0,0)= xi*G^M(0-)
       forall(i=0:L)G(ib)%lmix(1,i)=-xi*G(ib)%mats(L-i) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
    enddo
    !
    !DERIVATIVES
    allocate(GxG0(0:L))
    do ib=1,Nb
       !get d/dt G0^R = 0.d0
       dG0(ib)%ret(1)  = zero
       !
       !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
       do k=0,L
          GxG0(k)=G(ib)%lmix(1,k)*conjg(G0(ib)%lmix(1,L-k))
       end do
       dG0(ib)%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L) 
       !
       !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
       dG0(ib)%lmix(0:)= zero
       do j=0,L
          do k=0,j
             GxG0(k)=G(ib)%lmix(1,k)*G0(ib)%mats(k+L-j)
          end do
          dG0(ib)%lmix(j)=dG0(ib)%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
          do k=j,L
             GxG0(k)=G(ib)%lmix(1,k)*G0(ib)%mats(k-j)
          end do
          dG0(ib)%lmix(j)=dG0(ib)%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
       enddo
    enddo
    deallocate(GxG0)
    return
  end subroutine init_equilibrium_functions



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



end PROGRAM neqDMFT


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

