MODULE NEQ_CONTOUR_GF
  USE NEQ_CONTOUR
  !SCIFOR
  USE CONSTANTS, only: one,xi,zero,pi
  USE ARRAYS, only: linspace,arange
  USE IOTOOLS, only: reg,free_unit,splot,splot3d,read_data,store_data
  USE FUNCTIONS
  implicit none
  private


  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS:
  !====================================================
  type,public                          :: kb_contour_gf
     complex(8),dimension(:,:),pointer :: less
     complex(8),dimension(:,:),pointer :: ret
     complex(8),dimension(:,:),pointer :: lmix
     real(8),dimension(:),pointer      :: mats
     complex(8),dimension(:),pointer   :: iw
     logical                           :: status=.false.
     integer                           :: N=0
     integer                           :: L=0
     integer                           :: LF=0
  end type kb_contour_gf


  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS DERIVATIVE
  !====================================================
  type,public                          :: kb_contour_dgf
     complex(8),dimension(:),pointer   :: less,gtr
     complex(8),dimension(:),pointer   :: ret
     complex(8),dimension(:),pointer   :: lmix
     logical                           :: status=.false.
     integer                           :: N=0
     integer                           :: L=0
  end type kb_contour_dgf


  interface operator(*)
     module procedure &
          kb_contour_gf_scalarL_d,kb_contour_gf_scalarL_c,&
          kb_contour_dgf_scalarL_d,kb_contour_dgf_scalarL_c,&
          kb_contour_gf_scalarR_d,kb_contour_gf_scalarR_c,&
          kb_contour_dgf_scalarR_d,kb_contour_dgf_scalarR_c
  end interface operator(*)


  interface assignment(=)
     module procedure &
          kb_contour_gf_equality_,kb_contour_dgf_equality_,&
          kb_contour_gf_equality__,kb_contour_dgf_equality__
  end interface assignment(=)


  interface check_dimension_kb_contour
     module procedure &
          check_dimension_kb_contour_gf, &
          check_dimension_kb_contour_gf_, &
          check_dimension_kb_contour_dgf, &
          check_dimension_kb_contour_dgf_
  end interface check_dimension_kb_contour

  interface kb_trapz
     module procedure kb_trapz_d, kb_trapz_c
  end interface kb_trapz

  interface kb_half_trapz
     module procedure kb_half_trapz_d, kb_half_trapz_c
  end interface kb_half_trapz


  public :: allocate_kb_contour_gf
  public :: allocate_kb_contour_dgf
  public :: deallocate_kb_contour_gf
  public :: deallocate_kb_contour_dgf
  !
  public :: add_kb_contour_gf
  public :: convolute_kb_contour_gf
  public :: vie_kb_contour_gf
  public :: vide_kb_contour_gf
  !
  public :: gtr_kb_contour_gf
  public :: rmix_kb_contour_gf
  public :: save_kb_contour_gf
  public :: inquire_kb_contour_gf
  public :: read_kb_contour_gf
  public :: plot_kb_contour_gf
  public :: check_dimension_kb_contour
  !
  public :: kb_trapz,kb_half_trapz
  ! 
  public :: operator(*)
  public :: assignment(=)




contains





  !======= ALLOCATE ======= 
  subroutine allocate_kb_contour_gf(G,params)
    type(kb_contour_gf)     :: G
    type(kb_contour_params) :: params
    integer                 :: i,j,N,L,Lf
    nullify(G%less,G%ret,G%lmix,G%mats)
    N = params%Ntime            !<== allocate at maximum time
    L = params%Ntau
    Lf= params%Niw
    G%N = N
    G%L = L
    G%Lf= Lf
    Allocate(G%less(N,N))  ; G%less=zero
    allocate(G%ret(N,N))   ; G%ret=zero
    allocate(G%lmix(N,L))  ; G%lmix=zero
    allocate(G%mats(L))    ; G%mats=0d0
    allocate(G%iw(Lf))     ; G%iw=zero
    G%status=.true.
  end subroutine allocate_kb_contour_gf
  !
  subroutine allocate_kb_contour_dgf(dG,params,wgtr)
    type(kb_contour_dgf)    :: dG
    type(kb_contour_params) :: params
    integer                 :: i,j,N,L
    logical,optional        :: wgtr
    nullify(dG%less,dG%ret,dG%lmix)
    N=params%Ntime           !<== allocate at maximum time
    L=params%Ntau
    dG%N=N
    dG%L=L
    allocate(dG%less(N))  ; dG%less=zero
    allocate(dG%ret(N))   ; dG%ret=zero
    allocate(dG%lmix(L))  ; dG%lmix=zero
    if(present(wgtr).AND.wgtr)then
       allocate(dG%gtr(N))
       dG%gtr=zero
    endif
    dG%status=.true.
  end subroutine allocate_kb_contour_dgf




  !======= DEALLOCATE ======= 
  subroutine deallocate_kb_contour_gf(G)
    type(kb_contour_gf) :: G
    if(.not.G%status)stop "neq_contour_gf/deallocate_kb_contour_gf: G not allocated"
    deallocate(G%less,G%ret,G%lmix,G%mats,G%iw)
    G%N=0
    G%L=0
    G%Lf=0
    G%status=.false.
  end subroutine deallocate_kb_contour_gf
  !
  subroutine deallocate_kb_contour_dgf(dG)
    type(kb_contour_dgf) :: dG
    if(.not.dG%status)stop "contour_gf/deallocate_kb_contour_dgf: dG not allocated"
    deallocate(dG%less,dG%ret,dG%lmix)
    if(associated(dG%gtr))deallocate(dG%gtr)
    dG%N=0
    dG%L=0
    dG%status=.false.
  end subroutine deallocate_kb_contour_dgf








  !======= GET OTHER COMPONENTS ======= 
  function gtr_kb_contour_gf(G,N) result(Ggtr)
    type(kb_contour_gf)       :: G
    integer                   :: N
    complex(8),dimension(N,N) :: Ggtr
    integer                   :: i,j
    Ggtr = zero
    forall(i=1:N,j=1:N,i>=j)Ggtr(i,j) = G%less(i,j) + G%ret(i,j)
    forall(i=1:N,j=1:N,i<j)Ggtr(i,j)=-conjg(Ggtr(j,i))
  end function gtr_kb_contour_gf

  function rmix_kb_contour_gf(G,N,L) result(Grmix)
    type(kb_contour_gf)       :: G
    integer                   :: N,L
    complex(8),dimension(L,N) :: Grmix
    integer                   :: itau,j
    Grmix = zero
    forall(itau=1:L,j=1:N)Grmix(itau,j) = conjg(G%lmix(j,L-itau+1))
  end function rmix_kb_contour_gf










  !======= CHECK DIMENSION ======= 
  function check_dimension_kb_contour_gf(G,params) result(bool)
    type(kb_contour_gf)     :: G
    type(kb_contour_params) :: params
    logical                 :: bool
    integer                 :: N,L,Lf
    bool=.false.
    N=params%Ntime              !<== check size at maximum time
    L=params%Ntau
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         stop "neq_contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
    if( size(G%lmix)/=N*L )stop "neq_contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
    if( size(G%mats)/=L )stop "neq_contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
    bool=.true.
  end function check_dimension_kb_contour_gf
  !
  function check_dimension_kb_contour_gf_(G,N,L) result(bool)
    type(kb_contour_gf)     :: G
    logical                 :: bool
    integer                 :: N,L
    bool=.false.
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         stop "neq_contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
    if( size(G%lmix)/=N*L )stop "neq_contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
    if( size(G%mats)/=L )stop "neq_contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
    bool=.true.
  end function check_dimension_kb_contour_gf_
  !
  function check_dimension_kb_contour_dgf(dG,params) result(bool)
    type(kb_contour_dgf)    :: dG
    type(kb_contour_params) :: params
    logical                 :: bool
    integer                 :: N,L
    bool=.false.
    N=params%Ntime              !<== check size at maximum time
    L=params%Ntau
    if( (size(dG%less)/=N) .OR. (size(dG%ret)/=N) )&
         stop "neq_contour_gf/check_dimension_kb_contour_dgf: wrong dimensions less/ret"
    if( size(dG%lmix)/=L )&
         stop "neq_contour_gf/check_dimension_kb_contour_dgf: wrong dimensions lmix"
    bool=.true.
  end function check_dimension_kb_contour_dgf
  !
  function check_dimension_kb_contour_dgf_(dG,N,L) result(bool)
    type(kb_contour_dgf)    :: dG
    logical                 :: bool
    integer                 :: N,L
    bool=.false.
    if( (size(dG%less)/=N) .OR. (size(dG%ret)/=N) )&
         stop "neq_contour_gf/check_dimension_kb_contour_dgf: wrong dimensions less/ret"
    if( size(dG%lmix)/=L )&
         stop "neq_contour_gf/check_dimension_kb_contour_dgf: wrong dimensions lmix"
    bool=.true.
  end function check_dimension_kb_contour_dgf_



















  !======= SAVE ======= 
  subroutine save_kb_contour_gf(G,file)
    type(kb_contour_gf) :: G
    character(len=*)    :: file
    integer             :: unit
    if(.not.G%status)stop "neq_contour_gf/save_kb_contour_gf: G not allocated"
    unit = free_unit()
    open(unit,file=reg(file)//"_dimension.data")
    write(unit,"(A1,3A4,3X)")"#","N","L","Lf"
    write(unit,"(1X,3I4,3X)")G%N,G%L,G%Lf
    close(unit)
    call store_data(reg(file)//"_less.data",G%less(:,:))
    call store_data(reg(file)//"_ret.data", G%ret(:,:))
    call store_data(reg(file)//"_lmix.data",G%lmix(:,:))
    call store_data(reg(file)//"_mats.data",G%mats(:))
    call store_data(reg(file)//"_iw.data",G%iw(:))
  end subroutine save_kb_contour_gf
  !
  subroutine save_kb_contour_dgf(dG,file)
    type(kb_contour_dgf)  :: dG
    character(len=*)      :: file
    integer :: unit
    if(.not.dG%status)stop "neq_contour_gf/save_kb_contour_dgf: dG not allocated"
    unit = free_unit()
    open(unit,file=reg(file)//"_dimension.data")
    write(unit,"(A1,2A4,2X)")"#","N","L"
    write(unit,"(1X,2I4,2X)")dG%N,dG%L
    close(unit)
    call store_data(reg(file)//"_less.data",dG%less(:))
    call store_data(reg(file)//"_ret.data", dG%ret(:))
    call store_data(reg(file)//"_lmix.data",dG%lmix(:))
  end subroutine save_kb_contour_dgf









  !======= READ ======= 
  subroutine read_kb_contour_gf(G,file)
    type(kb_contour_gf)  :: G
    character(len=*)     :: file
    logical              :: check
    check = inquire_kb_contour_gf(file)
    if(.not.G%status.OR..not.check)stop "neq_contour_gf/read_kb_contour_gf: G not allocated"
    call read_data(trim(file)//"_less.data",G%less(:,:))
    call read_data(trim(file)//"_ret.data",G%ret(:,:))
    call read_data(trim(file)//"_lmix.data",G%lmix(:,:))
    call read_data(trim(file)//"_mats.data",G%mats(:))
    call read_data(trim(file)//"_iw.data",G%iw(:))
  end subroutine read_kb_contour_gf
  !
  subroutine read_kb_contour_dgf(dG,file)
    type(kb_contour_dgf) :: dG
    character(len=*)     :: file
    logical              :: check
    check = inquire_kb_contour_dgf(file)
    if(.not.dG%status.OR..not.check)stop "neq_contour_gf/read_kb_contour_dgf: dG not allocated"
    call read_data(trim(file)//"_less.data",dG%less(:))
    call read_data(trim(file)//"_ret.data",dG%ret(:))
    call read_data(trim(file)//"_lmix.data",dG%lmix(:))
  end subroutine read_kb_contour_dgf
  !
  function inquire_kb_contour_gf(file) result(check)
    integer                :: i
    logical                :: check,bool(5)
    character(len=*)       :: file
    character,dimension(5) :: ctype=([character(len=5) :: 'less','ret','lmix','mats','iw'])
    check=.true.
    do i=1,5
       inquire(file=reg(file)//"_"//reg(ctype(i))//".data",exist=bool(i))
       if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.gz",exist=bool(i))
       check=check.AND.bool(i)
    enddo
  end function inquire_kb_contour_gf
  !
  function inquire_kb_contour_dgf(file) result(check)
    integer                :: i
    logical                :: check,bool(3)
    character(len=*)       :: file
    character,dimension(3) :: ctype=([character(len=5) :: 'less','ret','lmix'])
    check=.true.
    do i=1,3
       inquire(file=reg(file)//"_"//reg(ctype(i))//".data",exist=bool(i))
       if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.gz",exist=bool(i))
       check=check.AND.bool(i)
    enddo
  end function inquire_kb_contour_dgf










  !======= PLOT ======= 
  subroutine plot_kb_contour_gf(file,G,params)
    character(len=*)        :: file
    type(kb_contour_gf)     :: G
    type(kb_contour_params) :: params
    integer                 :: Nt
    if(.not.G%status)stop "neq_contour_gf/plot_kb_contour_gf: G is not allocated" 
    Nt=params%Itime
    call splot3d(reg(file)//"_less_t_t.plot",params%t(:Nt),params%t(:Nt),G%less(:Nt,:Nt))
    call splot3d(reg(file)//"_ret_t_t.plot",params%t(:Nt),params%t(:Nt),G%ret(:Nt,:Nt))
    call splot3d(reg(file)//"_lmix_t_tau.plot",params%t(:Nt),params%tau(:),G%lmix(:Nt,:))
    call splot(reg(file)//"_mats_tau.plot",params%tau(:),G%mats(:))
    call splot(reg(file)//"_mats_iw.plot",params%wm(:),G%iw(:))
  end subroutine plot_kb_contour_gf
  !
  subroutine plot_kb_contour_dgf(file,dG,params)
    character(len=*)        :: file
    type(kb_contour_dgf)    :: dG
    type(kb_contour_params) :: params
    integer                 :: Nt
    if(.not.dG%status)stop "neq_contour_gf/plot_kb_contour_gf: G is not allocated" 
    Nt=params%Itime
    call splot(reg(file)//"_less_t.plot",params%t(:Nt),dG%less(:Nt))
    call splot(reg(file)//"_ret_t.plot",params%t(:Nt),dG%ret(:Nt))
    call splot(reg(file)//"_lmix_tau.plot",params%tau(:),dG%lmix(:))
  end subroutine plot_kb_contour_dgf











  !======= CONVOLUTION ======= 
  !C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max
  !t_max_index==N
  subroutine convolute_kb_contour_gf(A,B,C,params,dcoeff,ccoeff)
    type(kb_contour_gf)                   :: A,B,C,Cdum
    type(kb_contour_params)               :: params
    real(8),optional                      :: dcoeff
    complex(8),optional                   :: ccoeff
    integer                               :: N,L
    real(8)                               :: dt,dtau
    complex(8),dimension(:),allocatable   :: AxB    
    integer                               :: i,j,k,itau,jtau
    logical                               :: checkA,checkB,checkC
    complex(8),dimension(:,:),allocatable :: Amat,Bmat,Cmat
    integer                               :: Ntot
    if(  (.not.A%status).OR.&
         (.not.B%status).OR.&
         (.not.C%status))stop "neq_contour_gf/convolute_kb_contour_gf: A,B,C not allocated"
    N   = params%Itime      !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    dt  = params%dt
    dtau= params%dtau
    !
    checkA=check_dimension_kb_contour(A,params%Ntime,L) 
    checkB=check_dimension_kb_contour(B,params%Ntime,L)
    checkC=check_dimension_kb_contour(C,params%Ntime,L)
    !
    allocate(AxB(max(L,N)))
    !
    !Ret. component
    !C^R(t,t')=\int_{t'}^t ds A^R(t,s)*B^R(s,t')
    !-------------------------------------------------------------------
    C%ret(N,1:N)=zero
    do j=1,N
       AxB  = zero
       do k=j,N
          AxB(k) = A%ret(N,k)*B%ret(k,j)
       enddo
       C%ret(n,j) = C%ret(n,j) + dt*kb_trapz(AxB,j,N)
    enddo
    !
    !Lmix. component
    !C^\lmix(t,tau')=\int_0^{beta} ds A^\lmix(t,s)*B^M(s,tau')
    !               +\int_0^{t} ds A^R(t,s)*B^\lmix(s,tau')
    !               = I1 + I2
    !-------------------------------------------------------------------
    C%lmix(N,L)=zero
    do jtau=1,L
       !I1:
       AxB = zero
       do k=1,jtau
          AxB(k)=A%lmix(N,k)*B%mats(k + L-jtau)
       end do
       C%lmix(N,jtau)=C%lmix(N,jtau) - dtau*kb_trapz(AxB,1,jtau)
       do k=jtau,L
          AxB(k)=A%lmix(N,k)*B%mats(k - jtau+1)
       end do
       C%lmix(n,jtau)=C%lmix(n,jtau) + dtau*kb_trapz(AxB,jtau,L)
       !
       !I2:
       AxB = zero
       do k=1,N
          AxB(k) = A%ret(N,k)*B%lmix(k,jtau)
       enddo
       C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB,1,N)
    enddo
    !
    !Less component
    !C^<(t,t')=-i\int_0^{beta} ds A^\lmix(t,s)*B^\rmix(s,t')
    !           +\int_0^{t'} ds A^<(t,s)*B^A(s,t')
    !           +\int_0^{t} ds A^R(t,s)*B^<(s,t')
    !         = I1 + I2 + I3
    ! (t,t')=>(N,j) <==> Vertical side, no tip (j=1,N-1)
    !-------------------------------------------------------------------
    do j=1,N-1
       C%less(N,j)=zero
       !I1:
       do k=1,L
          AxB(k)=A%lmix(N,k)*conjg(B%lmix(j,L-k+1))
       end do
       C%less(N,j)=C%less(N,j) - xi*dtau*kb_trapz(AxB,1,L)
       !I2:
       do k=1,j
          AxB(k)=A%less(N,k)*conjg(B%ret(j,k))
       end do
       C%less(N,j)=C%less(N,j) + dt*kb_trapz(AxB,1,j)
       !I3:
       do k=1,N
          AxB(k)=A%ret(N,k)*B%less(k,j)
       end do
       C%less(N,j)=C%less(N,j) + dt*kb_trapz(AxB,1,N)
    end do
    ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
    do i=1,N
       C%less(i,N)=zero
       do k=1,L
          AxB(k)=A%lmix(i,k)*conjg(B%lmix(n,L-k+1))
       end do
       C%less(i,N)=C%less(i,N) - xi*dtau*kb_trapz(AxB,1,L)
       !
       do k=1,N
          AxB(k)=A%less(i,k)*conjg(B%ret(N,k))
       end do
       C%less(i,N)=C%less(i,N) + dt*kb_trapz(AxB,1,N)
       !
       do k=1,i
          AxB(k)=A%ret(i,k)*B%less(k,N)
       end do
       C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB,1,i)
    end do
    deallocate(AxB)
    !    
    if(present(dcoeff))then
       C%lmix(N,L)     = dcoeff*C%lmix(N,L)
       C%less(N,1:N-1) = dcoeff*C%less(N,1:N-1)
       C%less(1:N,N)   = dcoeff*C%less(1:N,N)
       C%ret(N,1:N)    = dcoeff*C%ret(N,1:N)
    endif
    if(present(ccoeff))then
       C%lmix(N,L)     = ccoeff*C%lmix(N,L)
       C%less(N,1:N-1) = ccoeff*C%less(N,1:N-1)
       C%less(1:N,N)   = ccoeff*C%less(1:N,N)
       C%ret(N,1:N)    = ccoeff*C%ret(N,1:N)
    endif
    !
  end subroutine convolute_kb_contour_gf





  !======= ADD ======= 
  !C(t,t')=A(t,t') + B(t,t'), with t=t_max && t'=0,t_max
  !t_max_index==N
  !IF YOU WANT THIS TO BE OVERLOADED IN + OPERATOR
  !C SHOULD BE ALLOCATED INSIDE THE FUNCTION ADD_KB...
  !BECAUSE POINTERS IS THE RESULT OF FUNCTION ITSELF (C).
  !ANYWAY THIS REQUIRE ONE TO CHECK ALL THE SIZES AND ALLOCATION
  !I RATHER PREFER TO USE THIS ROUTINE NOW...
  subroutine add_kb_contour_gf(A,B,C,params)
    type(kb_contour_gf)     :: A,B,C
    type(kb_contour_params) :: params
    integer                 :: N,L
    logical                 :: checkA,checkB,checkC
    if(  (.not.A%status).OR.&
         (.not.B%status).OR.&
         (.not.C%status))stop "neq_contour_gf/add_kb_contour_gf: A,B,C not allocated"
    N   = params%Ntime
    L   = params%Ntau
    !
    checkA=check_dimension_kb_contour(A,N,L) 
    checkB=check_dimension_kb_contour(B,N,L)
    checkC=check_dimension_kb_contour(C,N,L)
    C%less(:,:) = A%less(:,:) + B%less(:,:)
    C%ret(:,:)  = A%ret(:,:)  + B%ret(:,:)
    C%lmix(:,:) = A%lmix(:,:) + B%lmix(:,:)
    C%mats(:)   = A%mats(:)   + B%mats(:)
  end subroutine  add_kb_contour_gf

  subroutine add_kb_contour_dgf(A,B,C,params)
    type(kb_contour_dgf)    :: A,B,C
    type(kb_contour_params) :: params
    integer                 :: N,L
    logical                 :: checkA,checkB,checkC
    if(  (.not.A%status).OR.&
         (.not.B%status).OR.&
         (.not.C%status))stop "neq_contour_gf/add_kb_contour_gf: A,B,C not allocated"
    N   = params%Itime                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    checkA=check_dimension_kb_contour(A,N,L)
    checkB=check_dimension_kb_contour(B,N,L)
    checkC=check_dimension_kb_contour(C,N,L)
    C%less(:) = A%less(:) + B%less(:)
    C%ret(:)  = A%ret(:)  + B%ret(:)  
    C%lmix(:) = A%lmix(:) + B%lmix(:)
  end subroutine add_kb_contour_dgf











  !======= VOLTERRA INTEGRAL EQUATION ======= 
  !----------------------------------------------------------------------------
  !  This subroutine solves a Volterra integral equation of the second kind,
  !              G(t,t')+(K*G)(t,t')=Q(t,t')
  !  for t=n*dt or t'=n*dt, using 2^nd *implicit* Runge-Kutta method.
  !----------------------------------------------------------------------------
  subroutine vie_kb_contour_gf(Q,K,G,params)
    type(kb_contour_gf),intent(in)      :: Q
    type(kb_contour_gf),intent(in)      :: K
    type(kb_contour_gf),intent(inout)   :: G
    type(kb_contour_params),intent(in)  :: params
    integer                             :: N,L
    real(8)                             :: dt,dtau
    integer                             :: i,j,s,itau,jtau
    complex(8),dimension(:),allocatable :: KxG
    !
    N   = params%Itime                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    dt  = params%dt
    dtau= params%dtau
    !
    allocate(KxG(max(N,L)))
    !
    ! Ret component
    ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
    !-------------------------------------------------------------------
    G%ret(N,N)=Q%ret(N,N)
    do j=1,N-1
       G%ret(N,j)=Q%ret(N,j)
       do s=j,N-1
          KxG(s)=K%ret(N,s)*G%ret(s,j)
       eNd do
       G%ret(N,j)=G%ret(N,j) + dt*kb_half_trapz(KxG,j,N-1)
       G%ret(N,j)=G%ret(N,j)/(1d0 - 0.5d0*dt*K%ret(N,N))
    end do
    !
    ! Lmix component
    ! G^\lmix(t,tau') - \int_0^t K^R(t,s)*G^\lmix(s,tau')ds
    !    = Q^\lmix(t,tau') + \int_0^\beta K^\lmix(t,s)*G^M(s,tau')ds
    !-------------------------------------------------------------------
    do jtau=1,L
       G%lmix(N,jtau)=Q%lmix(N,jtau)
       do s=1,jtau
          KxG(s)=K%lmix(N,s)*G%mats(s + L-jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau) - dtau*kb_trapz(KxG,1,jtau)
       do s=jtau,L
          KxG(s)=K%lmix(N,s)*G%mats(s - jtau+1)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau)+dtau*kb_trapz(KxG,jtau,L)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau) + dt*kb_half_trapz(KxG,1,N-1)
       G%lmix(N,jtau)=G%lmix(N,jtau)/(1d0 - 0.5d0*dt*K%ret(N,N))
    end do
    !
    ! Less component
    ! G^<(t,t') - \int_0^t K^{R}(t,s)*G^{<}(s,t')ds
    !    = Q^<(t,t') - i\int_0^\beta K^\lmix(t,s)*G^\rmix(s,t')ds
    !      + \int_0^{t'} K^<(t,s)*G^A(s,t')ds
    !-------------------------------------------------------------------
    ! G^<(t_{N},t_{j})
    !-------------------------------------------------------------------
    do j=1,N-1
       G%less(N,j)=Q%less(N,j)
       do s=1,L
          KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s+1))
       enddo
       G%less(N,j)=G%less(N,j) - xi*dtau*kb_trapz(KxG,1,L)
       !
       do s=1,j
          KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
       enddo
       G%less(N,j)=G%less(N,j) + dt*kb_trapz(KxG,1,j)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%less(s,j)
       enddo
       G%less(N,j)=G%less(N,j) + dt*kb_half_trapz(KxG,1,N-1)
       G%less(N,j)=G%less(N,j)/(1d0 - 0.5d0*dt*K%ret(N,N))
    end do
    !-------------------------------------------------------------------
    ! G^<(t_{i},t_{N}) <= Hermite conjugate
    !-------------------------------------------------------------------
    do i=1,N-1
       G%less(i,N) = -conjg(G%less(N,i))
    end do
    !-------------------------------------------------------------------
    ! G^{<}(t_{n},t_{n}) <= Diagonal
    !-------------------------------------------------------------------
    G%less(N,N)=Q%less(N,N)
    do s=1,L
       KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s+1))
    end do
    G%less(N,N)=G%less(N,N) - xi*dtau*kb_trapz(KxG,1,L)
    !
    do s=1,N
       KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
    end do
    G%less(N,N)=G%less(N,N)+dt*kb_trapz(KxG,1,N)
    !
    do s=1,N-1
       KxG(s)=K%ret(N,s)*G%less(s,N)
    end do
    G%less(N,N)=G%less(N,N) + dt*kb_half_trapz(KxG,1,N-1)
    G%less(N,N)=G%less(N,N)/(1d0 - 0.5d0*dt*K%ret(N,N))
    !
    deallocate(KxG)
    !
  end subroutine vie_kb_contour_gf








  !======= VOLTERRA INTEGRO-DIFFERENTIAL EQUATION ======= 
  !----------------------------------------------------------------------------
  !  This subroutine solves a Volterra integro-differential equation of 
  !  the second kind,
  !              [i*d/dt-h(t)]G(t,t') = delta(t,t') + (K*G)(t,t'),
  !  for t=n*dt or t'=n*dt, using 2^nd implicit Runge-Kutta method.
  !----------------------------------------------------------------------------
  subroutine vide_kb_contour_gf(H,K,G,dG,dG_new,params)
    complex(8),dimension(:)             :: H
    type(kb_contour_gf),intent(in)      :: K
    type(kb_contour_gf),intent(inout)   :: G
    type(kb_contour_dgf),intent(inout)  :: dG
    type(kb_contour_dgf)                :: dG_new
    type(kb_contour_params),intent(in)  :: params
    integer                             :: N,L
    real(8)                             :: dt,dtau
    integer                             :: i,j,itau,jtau,s
    complex(8),dimension(:),allocatable :: KxG
    complex(8)                          :: dG_less,A
    !
    N   = params%Itime                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    dt  = params%dt
    dtau= params%dtau
    !
    allocate(KxG(max(N,L)))
    !
    !TYPE: X=RET,LESS,LMIX (R,<,\LMIX)
    !I D/DT G^X(T,:)  = H(T)G^X(T,:) + Q^X(T,:) + INT_{0,T'}^{T}K^R(T,S)G^X(S,:)DS
    !
    !Ret component
    ! d/dt G^R(t,:) = -i*h(t)*G^R(t,:) -i*delta(t,:) - i\int_{:}^t K^R(t,s)*G^R(s,:)ds
    !-------------------------------------------------------------------------------------
    G%ret(N,N)    = -xi
    dG_new%ret(N) = -xi*H(N)*G%ret(N,N)
    do j=1,N-1
       G%ret(N,j)=G%ret(N-1,j) + 0.5d0*dt*dG%ret(j)
       do s=j,N-1
          KxG(s)=K%ret(N,s)*G%ret(s,j)
       enddo
       dG_new%ret(j) = -xi*dt*kb_half_trapz(KxG,j,N-1)
       !
       G%ret(N,j)=G%ret(N,j) + 0.5d0*dt*dG_new%ret(j)
       G%ret(N,j)=G%ret(N,j)/(1d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
       !
       dG_new%ret(j)=dG_new%ret(j) - xi*H(N)*G%ret(N,j) - 0.5d0*xi*dt*K%ret(N,N)*G%ret(N,j)
    enddo
    !
    !Lmix component
    !d/dt G^\lmix(t,:) = -i*H(t)*G^\lmix(t,:) -i\int_0^\beta K^\lmix(t,s)*G^M(s,:)ds -i\int_0^t K^R(t,s)G^\lmix(s,:)ds
    !-------------------------------------------------------------------------------------
    do jtau=1,L
       G%lmix(N,jtau)=G%lmix(N-1,jtau) + 0.5d0*dt*dG%lmix(jtau)
       do s=1,jtau
          KxG(s)=K%lmix(N,s)*G%mats(s + L-jtau)
       end do
       dG_new%lmix(jtau) = xi*dtau*kb_trapz(KxG,1,jtau)
       do s=jtau,L
          KxG(s)=K%lmix(N,s)*G%mats(s-jtau+1)
       end do
       dG_new%lmix(jtau) = dG_new%lmix(jtau) - xi*dtau*kb_trapz(KxG,jtau,L)!<= add -iQ(t)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
       end do
       dG_new%lmix(jtau) = dG_new%lmix(jtau) - xi*dt*kb_half_trapz(KxG,1,N-1)
       !
       G%lmix(N,jtau)=G%lmix(N,jtau) + 0.5d0*dt*dG_new%lmix(jtau)
       G%lmix(N,jtau)=G%lmix(N,jtau)/(1d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
       !
       dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*H(N)*G%lmix(N,jtau) - 0.5d0*xi*dt*K%ret(N,N)*G%lmix(N,jtau)
    end do
    !
    !Less component
    !d/dt G^<(t,:) = -i*H(t)*G^<(t,:) -i*[ (-i)*\int_0^\beta K^\lmix(t,s)*G^\rmix(s,:)ds + \int_0^{:}K^<(t,s)*G^A(s,:)ds ]
    !                                 -i*\int_0^t K^R(t,s)*G^<(s,:)ds
    !
    ! G^<(t_{N},t_{j}), d/dt G^<(t_{N},t_{j}) <== lower-right triangle
    !-------------------------------------------------------------------------------------
    do j=1,N-1
       G%less(N,j)=G%less(N-1,j) + 0.5d0*dt*dG%less(j)
       do s=1,L
          KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s+1))
       end do
       dG_new%less(j) = -xi*(-xi)*dtau*kb_trapz(KxG,1,L)
       do s=1,j
          KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
       end do
       dG_new%less(j)=dG_new%less(j) - xi*dt*kb_trapz(KxG,1,j)!<= -iQ(t)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%less(s,j)
       end do
       dG_new%less(j) = dG_new%less(j) - xi*dt*kb_half_trapz(KxG,1,N-1)
       !
       G%less(N,j)=G%less(N,j) + 0.5d0*dt*dG_new%less(j)
       G%less(N,j)=G%less(N,j)/(1d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
       !
       dG_new%less(j)=dG_new%less(j) - xi*H(N)*G%less(N,j) - 0.5d0*xi*dt*K%ret(N,N)*G%less(N,j)
    end do
    !
    ! G^<(t_{i},t_{N}), d/dt G^<(t_{i},t_{N}) <== upper left triangle
    ! Hermitian conjugate G
    !-------------------------------------------------------------------------------------
    do i=1,N-1
       G%less(i,N)=-conjg(G%less(N,i))
    end do
    !
    ! d/dt G^<(t_{N-1},t_{N})
    !-------------------------------------------------------------------------------------
    dG_less=-xi*H(N-1)*G%less(N-1,N)
    do s=1,L
       KxG(s)=K%lmix(N-1,s)*conjg(G%lmix(N,L-s+1))
    end do
    dG_less=dG_less - xi*(-xi)*dtau*kb_trapz(KxG,1,L)
    do s=1,N
       KxG(s)=K%less(N-1,s)*conjg(G%ret(N,s))
    end do
    dG_less=dG_less - xi*dt*kb_trapz(KxG,1,N)
    do s=1,N-1
       KxG(s)=K%ret(N-1,s)*G%less(s,N)
    end do
    dG_less=dG_less - xi*dt*kb_trapz(KxG,1,N-1)
    !
    !G^<(N,N), d/dt G^<(N,N)
    !-------------------------------------------------------------------------------------
    G%less(N,N)=G%less(N-1,N) + 0.5d0*dt*dG_less
    do s=1,L
       KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s+1))
    end do
    dG_new%less(N) = -xi*(-xi)*dtau*kb_trapz(KxG,1,L)
    !
    do s=1,N
       KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
    end do
    dG_new%less(N) = dG_new%less(N) - xi*dt*kb_trapz(KxG,1,N)
    !
    do s=1,N-1
       KxG(s)=K%ret(N,s)*G%less(s,N)
    end do
    dG_new%less(N)=dG_new%less(N)-xi*dt*kb_half_trapz(KxG,1,N-1)
    !
    G%less(N,N)=G%less(N,N) + 0.5d0*dt*dG_new%less(N)
    G%less(N,N)=G%less(N,N)/(1d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
    !
    dG_new%less(N)=dG_new%less(N) - xi*H(N)*G%less(N,N) - 0.5d0*xi*dt*K%ret(N,N)*G%less(N,N)
    !
    !d/dt G <== d/dt G_new
    ! dG=dG_new
    ! call deallocate_kb_contour_dgf(dG_new)
    deallocate(KxG)
    !
  end subroutine vide_kb_contour_gf
















  !----------------------------------------------------------------------------
  !  This function calculates the sum
  !    \sum_{k=i}^{j} w_{k}^{i,j} f_{k}.
  !    w_{k}^{i,j} = 1/2  for k=i,j
  !                = 1    for i<k<j
  !                      OR
  !                = 1    for i<k<=j (half-edge)
  !----------------------------------------------------------------------------
  function kb_trapz_d(f,ia,ib) result(sum)
    real(8),dimension(:),intent(in) :: f
    integer,intent(in)              :: ia, ib
    integer                         :: k
    real(8)                         :: sum
    sum=0.d0
    if(ia==ib)then
       return
    else
       sum=sum+0.5d0*f(ia)
       do k=ia+1,ib-1
          sum=sum+f(k)
       end do
       sum=sum+0.5d0*f(ib)
    end if
  end function kb_trapz_d

  function kb_trapz_c(f,ia,ib,origin) result(sum)
    complex(8),dimension(:),intent(in) :: f
    integer,intent(in)                 :: ia, ib
    integer                            :: k
    complex(8)                         :: sum
    logical,optional,intent(in)        :: origin
    logical                            :: w0
    real(8)                            :: w

    w0=.true.
    if(present(origin)) w0=origin

    w = 0.5d0
    if(.not.w0) w = 1d0

    sum=zero
    if(ia==ib)then
       return
    else
       sum=sum+w*f(ia)
       do k=ia+1,ib-1
          sum=sum+f(k)
       end do
       sum=sum+w*f(ib)
    end if
  end function kb_trapz_c

  function kb_half_trapz_d(f,ia,ib) result(sum)
    real(8),dimension(:),intent(in) :: f
    integer,intent(in)              :: ia, ib
    integer                         :: k
    real(8)                         :: sum
    sum=0.5d0*f(ia)
    do k=ia+1,ib
       sum=sum+f(k)
    end do
  end function kb_half_trapz_d
  !
  function kb_half_trapz_c(f,ia,ib) result(sum)
    complex(8),dimension(:),intent(in) :: f
    integer,intent(in)                 :: ia, ib
    integer                            :: k
    complex(8)                         :: sum
    sum=0.5d0*f(ia)
    do k=ia+1,ib
       sum=sum+f(k)
    end do
  end function kb_half_trapz_c











  !======= OPERATIONS + & * ======= 
  subroutine kb_contour_gf_equality_(G1,C)
    type(kb_contour_gf),intent(inout) :: G1
    complex(8),intent(in)             :: C
    G1%less(:,:)  = C
    G1%ret(:,:)   = C
    G1%lmix(:,:)  = C
    G1%mats(:)    = C
    G1%iw(:)      = C
  end subroutine kb_contour_gf_equality_
  !
  subroutine kb_contour_gf_equality__(G1,G2)
    type(kb_contour_gf),intent(inout) :: G1
    type(kb_contour_gf),intent(in)    :: G2
    G1%less(:,:)  = G2%less(:,:)
    G1%ret(:,:)   = G2%ret(:,:)
    G1%lmix(:,:)  = G2%lmix(:,:)
    G1%mats(:)    = G2%mats(:)
    G1%iw(:)      = G2%iw(:)
  end subroutine kb_contour_gf_equality__
  !
  subroutine kb_contour_dgf_equality_(dG1,C)
    type(kb_contour_dgf),intent(inout) :: dG1
    complex(8),intent(in)             :: C
    dG1%less(:)  = C
    dG1%ret(:)   = C
    dG1%lmix(:)  = C
  end subroutine kb_contour_dgf_equality_
  !
  subroutine kb_contour_dgf_equality__(dG1,dG2)
    type(kb_contour_dgf),intent(inout) :: dG1
    type(kb_contour_dgf),intent(in)    :: dG2
    dG1%less(:)  = dG2%less(:)
    dG1%ret(:)   = dG2%ret(:)
    dG1%lmix(:)  = dG2%lmix(:)
  end subroutine kb_contour_dgf_equality__





  function kb_contour_gf_scalarL_d(C,G) result(F)
    real(8),intent(in)             :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = C*G%less(:,:)
    F%ret(:,:)  = C*G%ret(:,:)
    F%lmix(:,:) = C*G%lmix(:,:)
    F%mats(:)   = C*G%mats(:)
    F%iw(:)     = C*G%iw(:)
  end function kb_contour_gf_scalarL_d
  !
  function kb_contour_dgf_scalarL_d(C,dG) result(dF)
    real(8),intent(in)             :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = C*dG%less(:)
    dF%ret(:)  = C*dG%ret(:)
    dF%lmix(:) = C*dG%lmix(:)
  end function kb_contour_dgf_scalarL_d





  function kb_contour_gf_scalarL_c(C,G) result(F)
    complex(8),intent(in)          :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = C*G%less(:,:)
    F%ret(:,:)  = C*G%ret(:,:)
    F%lmix(:,:) = C*G%lmix(:,:)
    F%mats(:)   = C*G%mats(:)
    F%iw(:)     = C*G%iw(:)
  end function kb_contour_gf_scalarL_c
  !
  function kb_contour_dgf_scalarL_c(C,dG) result(dF)
    complex(8),intent(in)           :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = C*dG%less(:)
    dF%ret(:)  = C*dG%ret(:)
    dF%lmix(:) = C*dG%lmix(:)
  end function kb_contour_dgf_scalarL_c





  function kb_contour_gf_scalarR_d(G,C) result(F)
    real(8),intent(in)             :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = G%less(:,:)*C
    F%ret(:,:)  = G%ret(:,:)*C
    F%lmix(:,:) = G%lmix(:,:)*C
    F%mats(:)   = G%mats(:)*C
    F%iw(:)     = G%iw(:)*C
  end function kb_contour_gf_scalarR_d
  !
  function kb_contour_dgf_scalarR_d(dG,C) result(dF)
    real(8),intent(in)             :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = dG%less(:)*C
    dF%ret(:)  = dG%ret(:)*C
    dF%lmix(:) = dG%lmix(:)*C
  end function kb_contour_dgf_scalarR_d




  function kb_contour_gf_scalarR_c(G,C) result(F)
    complex(8),intent(in)             :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = G%less(:,:)*C
    F%ret(:,:)  = G%ret(:,:)*C
    F%lmix(:,:) = G%lmix(:,:)*C
    F%mats(:)   = G%mats(:)*C
    F%iw(:)     = G%iw(:)*C
  end function kb_contour_gf_scalarR_c
  !
  function kb_contour_dgf_scalarR_c(dG,C) result(dF)
    complex(8),intent(in)             :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = dG%less(:)*C
    dF%ret(:)  = dG%ret(:)*C
    dF%lmix(:) = dG%lmix(:)*C
  end function kb_contour_dgf_scalarR_c

END MODULE NEQ_CONTOUR_GF

