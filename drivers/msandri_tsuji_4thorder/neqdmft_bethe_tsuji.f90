!###################################################################
!PURPOSE  : Solve quantum quench dynamic in the 1-band Bethe lattice 
!AUTHORS  : Adriano Amaricci 
!###################################################################
program neqDMFT
  !USE ERROR
  !USE FFTGF
  USE NEQ_VARS_GLOBAL   !global variables, calls to 3rd library 
  USE NEQ_THERMOSTAT    !contains bath inizialization
  implicit none
  integer                             :: i,j,ik,itime,iloop,ix,iy
  logical                             :: converged
  character(len=16)                   :: finput
  type(kb_contour_gf)                 :: Sbath
  type(kb_contour_gf)                 :: Gwf
  type(kb_contour_gf)                 :: Sigma
  type(kb_contour_gf)                 :: Gloc
  type(kb_contour_gf)                 :: Ker
  type(kb_contour_dgf)                :: dGwf,dGwf_old
  complex(8),dimension(:),allocatable :: Ham
  !RESULTS:
  real(8),dimension(:),allocatable    :: dos,epsi
  real(8),dimension(:,:),allocatable  :: nk

  !READ THE INPUT FILE (in vars_global):
  !=====================================================================
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.in')
  call parse_cmd_variable(KMAT,"KMAT",default=.false.)
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
  cc_params%Nt=1
  call init_equilibrium_functions(Gwf,dGwf,Gloc,Sigma,cc_params)
  call measure_observables(Gloc,Sigma,cc_params)


  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  !=====================================================================
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Nt=itime
     !
     call setup_weiss_field(Gwf,cc_params)
     dGwf_old = dGwf
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
     call measure_observables(Gloc,Sigma,cc_params)
  enddo




  print*,"Getting n(e,t):"
  allocate(epsi(Nx),dos(Nx))
  allocate(nk(cc_params%Ntime,Nx))
  call bethe_lattice(dos,epsi,Nx,2.d0*ts)
  !Gwf,dGwf <--- G_e(t,t'),d/dtG_e(t,t')
  do ik=1,Nx
     ham = epsi(ik)
     call setup_initial_conditions(Gwf,dGwf_old,Sigma,ik,cc_params)
     do itime=2,cc_params%Ntime
        cc_params%Nt = itime
        call vide_kb_contour_gf(Ham,Sigma,Gwf,dGwf_old,dGwf,cc_params)
        dGwf_old = dGwf
     enddo
     forall(itime=1:cc_params%Ntime)nk(itime,ik)=dimag(Gwf%less(itime,itime))
  enddo
  call splot3d("nkVSepsikVStime.plot",cc_params%t,epsi,nk,wlines=.true.)


  call test_(Gloc,Sigma,cc_params)

  print*,"BRAVO"

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !+-------------------------------------------------------------------+
  subroutine setup_initial_conditions(Gk,dGk,Self,ik,params)
    type(kb_contour_gf)                 :: Gk,Self
    type(kb_contour_dgf)                :: dGk
    integer                             :: ik
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    complex(8),dimension(Lfreq)         :: Gkiw
    real(8),dimension(Lfreq)            :: wm
    real(8)                             :: nk
    complex(8)                          :: epsk
    complex(8),allocatable,dimension(:) :: SxG
    L = params%Ntau
    Lf= params%Lf
    do i=1,Lf
       Gk%iw(i) = one/(xi*params%wm(i) - epsi(ik) - Self%iw(i))
    enddo
    call fftgf_iw2tau(Gk%iw,Gk%mats(0:),beta)        
    Gk%less(1,1) = -xi*Gk%mats(L)
    Gk%ret(1,1)  = -xi
    forall(i=0:L)Gk%lmix(1,i)=-xi*Gk%mats(L-i)
    !
    !Derivatives
    allocate(SxG(0:L))
    !
    !get d/dt G_k^R = -i e(k,0)G_k^R
    dGk%ret(1)  = -xi*epsi(ik)*Gk%ret(1,1)
    !
    !get d/dt G_k^< = -i e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    do k=0,L
       SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k))
    end do
    dGk%less(1) = -xi*epsi(ik)*Gk%less(1,1)-&
         xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
    !
    !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(0:)= -xi*epsi(ik)*Gk%lmix(1,0:)
    do j=0,L
       do k=0,j
          SxG(k)=self%lmix(1,k)*Gk%mats(k+L-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       do k=j,L
          SxG(k)=self%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
    enddo
  end subroutine setup_initial_conditions




  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine init_equilibrium_functions(g0,dg0,g,self,params)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_dgf)                :: dg0
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Lf
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: GxG0
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.dg0%status)stop "init_functions: dg0 is not allocated"
    if(.not.g%status)stop "init_functions: g is not allocated"
    if(.not.self%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    N = params%Nt
    L = params%Ntau
    Lf= params%Lf
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
          read(unit,*)wm,ims,res
          g0%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
    else
       write(*,"(A)")"Start from Non-interacting G0(iw)"
       do i=1,Lf
          wm    = pi/beta*dble(2*i-1)
          zeta  = xi*wm
          g0%iw(i) = gfbethe(wm,zeta,2.d0*ts)
       enddo
    endif
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    call fftgf_iw2tau(g0%iw,g0%mats(0:),params%beta)
    g0%less(1,1) = -xi*g0%mats(L)
    g0%ret(1,1)  = -xi
    forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    !(this step depends on the imp. solv.)
    ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
    ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
    ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
    ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
    ! self^R(0,0) = self^> - self^<
    do j=0,L
       Self%mats(j) = Ui*Ui*g0%mats(j)*g0%mats(L-j)*g0%mats(j)
    end do
    call fftgf_tau2iw(Self%mats(0:),Self%iw,beta)
    Self%iw = xi*dimag(self%iw) !!ACTHUNG: imposing half-filling symmetry
    Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(1)*g0%mats(L)
    Self_gtr      =(xi**3)*U*U*g0%mats(1)*g0%mats(L)*g0%mats(1)
    do j=0,L
       Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
    end do
    Self%ret(1,1) = Self_gtr - Self%less(1,1)
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    do i=1,Lf
       wm    = pi/beta*dble(2*i-1)
       zeta  = xi*wm - Self%iw(i)
       G%iw(i) = gfbethe(wm,zeta,2.d0*ts)
    enddo
    call fftgf_iw2tau(G%iw,G%mats(0:),beta)         !get G(tau)
    G%ret(1,1)  = -xi                               !get G^R(0,0)=-xi
    G%less(1,1) = -xi*G%mats(L)                  !get G^<(0,0)= xi*G^M(0-)
    forall(i=0:L)G%lmix(1,i)=-xi*G%mats(L-i) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
    !Derivatives
    allocate(GxG0(0:L))
    !
    !get d/dt G0^R = 0.d0
    dG0%ret(1)  = zero
    !
    !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
    do k=0,L
       GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k))
    end do
    dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L) 
    !
    !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
    dG0%lmix(0:)= zero
    do j=0,L
       do k=0,j
          GxG0(k)=G%lmix(1,k)*G0%mats(k+L-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
       do k=j,L
          GxG0(k)=G%lmix(1,k)*G0%mats(k-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
    enddo
    deallocate(GxG0)
    return
  end subroutine init_equilibrium_functions



  !+-------------------------------------------------------------------+
  !PURPOSE: setup the Weiss Field G0 for the next time-step
  !+-------------------------------------------------------------------+
  subroutine setup_weiss_field(g0,params)
    type(kb_contour_gf)                   :: g0
    type(kb_contour_params)               :: params
    integer                               :: i,j,k,N,L
    complex(8),allocatable,dimension(:) :: SxG
    if(.not.g0%status)stop "init_g0: g0 is not allocated"
    if(.not.params%status)stop "init_g0: params is not allocated"
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
          g0%ret(N,j) =g0%ret(1,1)
          g0%less(N,j)=g0%less(1,1)
       end do
       do i=1,N-1
          g0%less(i,N)=g0%less(1,1)
       end do
       do j=0,L
          g0%lmix(N,j)=g0%lmix(1,j)
       end do
    case default
       !EXTEND G0 FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
       !USING QUADRATIC EXTRAPOLATION
       do k=1,N-1
          g0%less(N,k)=2.d0*g0%less(N-1,k)-g0%less(N-2,k)
          g0%less(k,N)=2.d0*g0%less(k,N-1)-g0%less(k,N-2)
       end do
       g0%less(N,N)=2.d0*g0%less(N-1,N-1)-g0%less(N-2,N-2)
       !
       do k=0,L
          g0%lmix(N,k)=2.d0*g0%lmix(N-1,k)-g0%lmix(N-2,k)
       end do
       !
       g0%ret(N,N)=-xi
       do k=1,N-2
          g0%ret(N,k)=2.d0*g0%ret(N-1,k)-g0%ret(N-2,k)
       end do
       g0%ret(N,N-1)=0.5d0*(g0%ret(N,N)+g0%ret(N,N-2))
    end select
  end subroutine setup_weiss_field



  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve with the 2^nd IPT sigma functions
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt(G0,Sigma,params)
    type(kb_contour_gf)                   :: G0
    type(kb_contour_gf)                   :: Sigma,Sigma4(4)
    type(kb_contour_params)               :: params
    integer                               :: N,L
    complex(8),dimension(:,:),allocatable :: G0_gtr,Sigma_gtr,G0_rmix
    integer                               :: i,j,itau
    integer                               :: Ntot
    real(8),dimension(:),allocatable      :: Utime
    integer                               :: ia,ib,ic,Cont(3)
    complex(8)                            :: int_ij,intb
    complex(8),dimension(:,:),allocatable :: Smat,G0mat,G2mat,G3mat,Sfoo,Chimat
    complex(8),dimension(:),allocatable   :: tdiff

    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    Cont= (/N,2*N,L/)
    !
    allocate(G0_gtr(N,N),Sigma_gtr(N,N),G0_rmix(0:L,N))
    do j=1,N
       G0_gtr(N,j)=G0%less(N,j)+G0%ret(N,j)
    end do
    do i=1,N-1
       G0_gtr(i,N)=G0%less(i,n)-conjg(G0%ret(N,i))
    end do
    do j=0,L
       G0_rmix(j,N)  = conjg(G0%lmix(N,L-j))
    enddo
    !
    !Vertical edge
    do j=1,N
       Sigma%less(N,j)= U*U*G0%less(N,j)*G0_gtr(j,N)*G0%less(N,j)
       Sigma_gtr(N,j) = U*U*G0_gtr(N,j)*G0%less(j,N)*G0_gtr(N,j)
    end do
    !Horizontal edge
    do i=1,N-1
       Sigma%less(i,N)= U*U*G0%less(i,N)*G0_gtr(N,i)*G0%less(i,N)
       Sigma_gtr(i,N) = U*U*G0_gtr(i,N)*G0%less(N,i)*G0_gtr(i,N)
    end do
    !Imaginary time edge:
    forall(i=0:L)Sigma%lmix(N,i)  = U*Ui*G0%lmix(N,i)*G0_rmix(i,N)*G0%lmix(N,i)
    forall(j=1:N)Sigma%ret(N,j) = Sigma_gtr(N,j) - Sigma%less(N,j)

    !#################################
    

    if(ifourth)then
       print*,"Entering 4th order..."
       do i=1,4
          call allocate_kb_contour_gf(sigma4(i),params)
          sigma4(i) = zero
       enddo
       Ntot = 2*N+L+1

       allocate(G0mat(Ntot,Ntot),G2mat(Ntot,Ntot),G3mat(Ntot,Ntot))
       allocate(Chimat(Ntot,Ntot),Smat(Ntot,Ntot),Sfoo(Ntot,Ntot))
       allocate(tdiff(Ntot),Utime(Ntot))
       tdiff(1:N)           =  params%dt
       tdiff(N+1:2*N)       = -params%dt
       tdiff(2*N+1:2*N+L+1) = -xi*params%dtau
       Utime(1:N)           = U 
       Utime(N+1:2*N)       = U
       Utime(2*N+1:2*N+L+1) = Ui


       call calc_sigma_4th(G0,Sigma4,params)

       ! !Sigma^(4a):
       ! print*,"get 4a..."
       ! call kb_contour_gf2kb_matrix(G0,N,L,G0mat)
       ! do i=1,Ntot
       !    do j=1,Ntot
       !       G2mat(i,j) = G0mat(i,j)*G0mat(i,j)
       !       G3mat(i,j) = G0mat(i,j)*G0mat(i,j)*G0mat(i,j)
       !    enddo
       ! enddo


       ! ! print*,"Get 4a"
       ! Smat  = zero
       ! do ic=1,3!Ntot
       !    i = cont(ic)
       !    !$OMP PARALLEL DO PRIVATE(j,int_ij,ia,intb,ib)
       !    do j=1,Ntot
       !       if (j.ne.N .and. j.ne.2*N .and. j.ne.L) then
       !          int_ij=zero
       !          do ia=1,Ntot
       !             intb=zero
       !             do ib=1,Ntot
       !                intb = intb + Utime(i)*Utime(j)*Utime(ia)*Utime(ib)* &
       !                     G2mat(i,ia)*G2mat(ia,ib)*G2mat(ib,j)*tdiff(ib)
       !             enddo
       !             int_ij = int_ij + intb*tdiff(ia)
       !          enddo
       !          Smat(i,j) = G0mat(i,j)*int_ij
       !       endif
       !    enddo
       !    !$OMP END PARALLEL DO
       ! enddo
       ! do ic=1,3
       !    j=cont(ic)
       !    !$OMP PARALLEL DO PRIVATE(i,int_ij,ia,intb,ib)
       !    do i=1,Ntot
       !       int_ij=zero
       !       do ia=1,Ntot
       !          intb=zero
       !          do ib=1,Ntot
       !             intb = intb + Utime(i)*Utime(j)*Utime(ia)*Utime(ib)* &
       !                  G2mat(i,ia)*G2mat(ia,ib)*G2mat(ib,j)*tdiff(ib)
       !          enddo
       !          int_ij = int_ij + intb*tdiff(ia)
       !       enddo
       !       Smat(i,j) = G0mat(i,j)*int_ij
       !    enddo
       !    !$OMP END PARALLEL DO
       ! enddo

       ! print*,"extract Sigma^a:"
       ! call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(1))

       ! !Sigma^(4b):
       ! print*,"Get 4b"
       ! Smat  = zero
       ! do ic=1,3!i=1,Ntot
       !    i=cont(ic)
       !    !$OMP PARALLEL DO PRIVATE(j,int_ij,ia,intb,ib)
       !    do j=1,Ntot
       !       if (j.ne.N .and. j.ne.2*N .and. j.ne.L) then
       !          int_ij=zero
       !          do ia=1,Ntot
       !             intb=zero
       !             do ib=1,Ntot
       !                intb = intb + Utime(i)*Utime(j)*Utime(ia)*Utime(ib)* &
       !                     G0mat(i,ia)*G0mat(ib,j)*G3mat(ia,ib)*tdiff(ib)
       !             enddo
       !             int_ij = int_ij + intb*tdiff(ia)
       !          enddo
       !          Smat(i,j) = G2mat(i,j)*int_ij
       !       endif
       !    enddo
       !    !$OMP END PARALLEL DO
       ! enddo
       ! do ic=1,3
       !    j=cont(ic)
       !    !$OMP PARALLEL DO PRIVATE(i,int_ij,ia,intb,ib)
       !    do i=1,Ntot
       !       int_ij=zero
       !       do ia=1,Ntot
       !          intb=zero
       !          do ib=1,Ntot
       !             intb = intb + Utime(i)*Utime(j)*Utime(ia)*Utime(ib)* &
       !                  G0mat(i,ia)*G0mat(ib,j)*G3mat(ia,ib)*tdiff(ib)
       !          enddo
       !          int_ij = int_ij + intb*tdiff(ia)
       !       enddo

       !       Smat(i,j) = G2mat(i,j)*int_ij
       !    enddo
       !    !$OMP END PARALLEL DO
       ! enddo

       ! print*,"extract Sigma^b:"
       ! call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(2))


       ! !Sigma^(4d):
       ! print*,"Get 4d"
       ! Smat  = zero
       ! do ic=1,3!i=1,Ntot
       !    i=cont(ic)
       !    !$OMP PARALLEL DO PRIVATE(j,int_ij,ia,intb,ib)
       !    do j=1,Ntot
       !       if (j.ne.N .and. j.ne.2*N .and. j.ne.L) then
       !          int_ij=zero
       !          do ia=1,Ntot
       !             intb=zero
       !             do ib=1,Ntot
       !                intb = intb + Utime(i)*Utime(j)*Utime(ia)*Utime(ib)* &
       !                     G0mat(i,ia)*G0mat(ia,j)*G0mat(i,ib)*G0mat(ib,j)*G2mat(ia,ib)*tdiff(ib)
       !             enddo
       !             int_ij = int_ij + intb*tdiff(ia)
       !          enddo
       !          Smat(i,j) = G0mat(i,j)*int_ij
       !       endif
       !    enddo
       !    !$OMP END PARALLEL DO
       ! enddo
       ! do ic=1,3
       !    j=cont(ic)
       !    !$OMP PARALLEL DO PRIVATE(i,int_ij,ia,intb,ib)
       !    do i=1,Ntot
       !       int_ij=zero
       !       do ia=1,Ntot
       !          intb=zero
       !          do ib=1,Ntot
       !             intb = intb + Utime(i)*Utime(j)*Utime(ia)*Utime(ib)* &
       !                  G0mat(i,ia)*G0mat(ia,j)*G0mat(i,ib)*G0mat(ib,j)*G2mat(ia,ib)*tdiff(ib)
       !          enddo
       !          int_ij = int_ij + intb*tdiff(ia)
       !       enddo
       !       Smat(i,j) = G0mat(i,j)*int_ij
       !    enddo
       !    !$OMP END PARALLEL DO
       ! enddo
       ! print*,"extract Sigma^d:"
       ! call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(4))
       !============================================

       print*,"Sum up Sigma2+Sigma4:"
       do j=1,N
          Sigma%less(N,j) = Sigma%less(N,j) + &
               3.d0*(Sigma4(1)%less(N,j) + Sigma4(2)%less(N,j) + Sigma4(3)%less(N,j) - Sigma4(4)%less(N,j))
          Sigma%ret(N,j) = Sigma%ret(N,j) + &
               3.d0*(Sigma4(1)%ret(N,j) + Sigma4(2)%ret(N,j) + Sigma4(3)%ret(N,j) - Sigma4(4)%ret(N,j))
       enddo
       do i=1,N-1
          Sigma%less(i,N) = Sigma%less(i,N) + &
               3.d0*(Sigma4(1)%less(i,N) + Sigma4(2)%less(i,N) + Sigma4(3)%less(i,N) - Sigma4(4)%less(i,N))          
       enddo
       !!ACTHUNG!! U_i might not be correct in this expression!!
       do j=0,L
          Sigma%lmix(N,j) = Sigma%lmix(N,j) + &
               3.d0*Ui*(Sigma4(1)%lmix(N,j) + Sigma4(2)%lmix(N,j) + Sigma4(3)%lmix(N,j) - Sigma4(4)%lmix(N,j))
       enddo


       print*,"Deallocating:"
       do i=1,4
          call deallocate_kb_contour_gf(sigma4(i))
       enddo
       deallocate(G0mat,G2mat,G3mat)
       deallocate(Chimat,Smat,Sfoo)
       deallocate(tdiff)
    endif


  end subroutine neq_solve_ipt

  include "sigma_4th_v1.f90"









  subroutine test_(G,Self,params)
    type(kb_contour_gf)                   :: G,Self,GxS
    type(kb_contour_params)               :: params
    integer                               :: N,L,i,j,Ntot
    complex(8),dimension(:,:),allocatable :: Gmat,Smat,GxSmat
    real(8),dimension(:),allocatable      :: time
    N=params%Ntime
    L=params%Ntau
    Ntot=2*N+L+1
    allocate(time(Ntot))
    allocate(Gmat(Ntot,Ntot),Smat(Ntot,Ntot),GxSmat(Ntot,Ntot))
    time(1:N)=params%t(1:N)
    time(N+1:2*N)=time(N) + params%t(1:N)
    time(2*N+1:2*N+L+1)=time(2*N)+params%tau(0:L)/beta
    call kb_contour_gf2kb_matrix(G,N,L,Gmat)
    call splot3d("Gmat_test.plot",time,time,Gmat)
    call splot3d("Gmat_13.plot",params%t(1:N),params%tau(0:L),Gmat(1:N,2*N+1:2*N+1+L))
    call splot3d("Gmat_12.plot",params%t(1:N),params%t(1:N),Gmat(1:N,N+1:2*N))
    G=zero
    call kb_matrix2kb_contour_gf(Gmat,N,L,G)
    call plot_kb_contour_gf("G_test",G,params)


    !
    !get convolution in ordinary way:
    call allocate_kb_contour_gf(GxS,params)
    do i=1,params%Ntime
       params%Nt=i
       call convolute_kb_contour_gf(G,Self,GxS,params)
    enddo
    call plot_kb_contour_gf("GxS",GxS,params)
    !
    !map convolution to a matrix:
    call kb_contour_gf2kb_matrix(GxS,N,L,Smat)
    call splot3d("GxSmat_test1.plot",time,time,Smat)
    !
    !Evaluate convolution using Keldysh Matrices
    call kb_contour_gf2kb_matrix(Self,N,L,Smat)
    call convolute_kb_matrix_gf(Gmat,Smat,N,L,params,GxSmat)
    call splot3d("GxSmat_test2.plot",time,time,GxSmat)
    call splot3d("GxS_13.plot",params%t(1:N),params%tau(0:L),GxSmat(1:N,2*N+1:2*N+1+L))
    call splot3d("GxS_31.plot",params%tau(0:L),params%t(1:N),GxSmat(2*N+1:2*N+1+L,1:N))


    GxS=zero
    call kb_matrix2kb_contour_gf(GxSmat,N,L,GxS)
    call plot_kb_contour_gf("GxS_test",GxS,params)
  end subroutine test_




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



  !+-------------------------------------------------------------------+
  !PURPOSE: measure some observables and print them
  !+-------------------------------------------------------------------+
  subroutine measure_observables(g,self,params)
    type(kb_contour_gf)     :: g
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    integer                 :: unit,itime
    real(8)                 :: dens,docc,ekin,epot,etot
    itime = params%Nt
    unit = free_unit()
    open(unit,file="observables.plot",position="append")
    dens = measure_dens(g,self,params)
    docc = measure_docc(g,self,params)
    ekin = measure_ekin(g,self,params)
    epot = measure_epot(g,self,params)
    etot = ekin + epot
    write(unit,"(6F20.12)")params%t(itime),dens,docc,ekin,epot,etot
    close(unit)
  end subroutine measure_observables



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the density at a given istant of time
  ! n(t)=-xi*G^<(t,t)
  !+-------------------------------------------------------------------+
  function measure_dens(g,self,params) result(dens)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: dens
    integer                             :: N
    N = params%Nt
    dens = dimag(G%less(N,N))
  end function measure_dens


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the double occupancy at a given istant of time
  ! d(t)=n_up(t)*n_do(t)-1/U0*[Self^M*G^M]
  !      n_up(t)*n_do(t)-i/U*[Self^R*G^< + Self^<*G^A + Self^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_docc(g,self,params) result(docc)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: docc
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: SxG
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    nt   = dimag(G%less(N,N))
    allocate(SxG(0:max(N,L)))
    docc = nt**2
    if(N==1)then
       if(ui/=0.d0)then
          do k=0,L
             SxG(k)=Self%mats(L-i)*G%mats(i)
          end do
          docc=docc-1.d0/Ui*params%dtau*kb_trapz(SxG(0:),0,L)
       endif
    else
       if(u/=0.d0)then
          do k=0,L
             SxG(k)=Self%lmix(N,k)*conjg(G%lmix(N,L-k))
          end do
          docc=docc + 1.d0/U*params%dtau*dimag( (-xi)*kb_trapz(SxG(0:),0,L) )
          do k=1,N
             SxG(k)=Self%ret(N,k)*G%less(k,N)
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
          do k=1,N
             SxG(k)=Self%less(N,k)*conjg(G%ret(N,k))
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
       endif
    endif
    deallocate(SxG)
  end function measure_docc


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_ekin(g,self,params) result(ekin)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: ekin
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: Ker
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    allocate(Ker(0:max(N,L)))
    if(N==1)then
       do k=0,L
          Ker(k)=G%mats(L-k)*G%mats(k)
       end do
       ekin = -2.d0*params%dtau*kb_trapz(Ker(0:),0,L)
    else
       do k=0,L
          Ker(k)=G%lmix(N,k)*conjg(G%lmix(N,L-k))
       end do
       ekin=2.d0*params%dtau*dimag( (-xi)*kb_trapz(Ker(0:),0,L) )
       do k=1,N
          Ker(k)=G%ret(N,k)*G%less(k,N)
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
       do k=1,N
          Ker(k)=G%less(N,k)*conjg(G%ret(N,k))
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
    endif
    deallocate(Ker)
  end function measure_ekin



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! U(t)= U*docc(t) - n(t) + 1/4
  !+-------------------------------------------------------------------+
  function measure_epot(g,self,params) result(epot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: epot,docc,nt
    integer                             :: i,k,j,N,L
    N = params%Nt
    L = params%Ntau
    !
    if(N==1)then
       nt   = measure_dens(g,self,params)
       docc = measure_docc(g,self,params)
       epot = Ui*(docc - nt + 0.25d0)
    else
       nt   = measure_dens(g,self,params)
       docc = measure_docc(g,self,params)
       epot = U*(docc - nt + 0.25d0)
    endif
  end function measure_epot



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_etot(g,self,params) result(etot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: etot,ekin,epot
    ekin = measure_ekin(g,self,params)
    epot = measure_epot(g,self,params)
    etot = ekin + epot
  end function measure_etot



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

