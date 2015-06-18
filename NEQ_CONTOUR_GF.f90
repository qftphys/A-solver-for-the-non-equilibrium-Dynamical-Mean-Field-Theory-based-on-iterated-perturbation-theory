MODULE NEQ_CONTOUR_GF
  USE NEQ_CONTOUR
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS
  !USE SF_SPECIAL
  implicit none
  private

  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS:
  !====================================================
  type,public :: kb_contour_gf
     complex(8),dimension(:,:),allocatable :: less
     complex(8),dimension(:,:),allocatable :: ret
     complex(8),dimension(:,:),allocatable :: lmix
     real(8),dimension(:),allocatable      :: mats
     complex(8),dimension(:),allocatable   :: iw
     logical                               :: status=.false.
     logical                               :: anomalous=.false.
     integer                               :: N=0
     integer                               :: L=0
     integer                               :: LF=0
  end type kb_contour_gf
  !
  type,public :: kb_contour_sigma
     type(kb_contour_gf)                   :: self
     complex(8),dimension(:),allocatable   :: hfb
  end type kb_contour_sigma
  !
  !
  !
  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS DERIVATIVE
  !====================================================
  type,public :: kb_contour_dgf
     complex(8),dimension(:),allocatable :: less,gtr
     complex(8),dimension(:),allocatable :: ret
     complex(8),dimension(:),allocatable :: lmix
     logical                             :: status=.false.
     logical                             :: anomalous=.false.
     integer                             :: N=0
     integer                             :: L=0
  end type kb_contour_dgf
  !
  !
  !
  !
  !ALLOCATION ROUTINES:
  interface allocate_kb_contour_gf
     module procedure allocate_kb_contour_gf_main
     module procedure allocate_kb_contour_gf_nambu_redux
     module procedure allocate_kb_contour_gf_nambu
  end interface allocate_kb_contour_gf
  !
  interface allocate_kb_contour_dgf
     module procedure allocate_kb_contour_dgf_main
     module procedure allocate_kb_contour_dgf_nambu_redux
  end interface allocate_kb_contour_dgf
  !
  public :: allocate_kb_contour_gf
  public :: allocate_kb_contour_dgf
  !
  !
  !DEALLOCATION ROUTINES:
  public :: deallocate_kb_contour_gf
  public :: deallocate_kb_contour_dgf
  !
  !
  !CHECK DIMENSIONS:
  interface check_dimension_kb_contour
     module procedure check_dimension_kb_contour_gf
     module procedure check_dimension_kb_contour_gf_
     module procedure check_dimension_kb_contour_dgf
     module procedure check_dimension_kb_contour_dgf_
  end interface check_dimension_kb_contour
  public :: check_dimension_kb_contour
  !
  !
  !ADD (TOTAL DOMAIN) ROUTINES:
  interface add_kb_contour_gf
     module procedure add_kb_contour_gf_simple
     module procedure add_kb_contour_gf_recursive
     module procedure add_kb_contour_gf_delta_d
     module procedure add_kb_contour_gf_delta_c
  end interface add_kb_contour_gf
  public :: add_kb_contour_gf
  public :: add_kb_contour_dgf
  !
  !
  !SUM (PERIMETER) ROUTINES:
  interface sum_kb_contour_gf
     module procedure sum_kb_contour_gf_simple
     module procedure sum_kb_contour_gf_recursive
     module procedure sum_kb_contour_gf_delta_d
     module procedure sum_kb_contour_gf_delta_c
  end interface sum_kb_contour_gf
  public :: sum_kb_contour_gf
  !
  !
  !DELETE (RESET PERIMETER) ROUTINES:
  public :: del_kb_contour_gf
  !
  !CONVOLUTION:
  interface convolute_kb_contour_gf
     module procedure convolute_kb_contour_gf_simple
     module procedure convolute_kb_contour_gf_recursive
     module procedure convolute_kb_contour_gf_delta_left
     module procedure convolute_kb_contour_gf_delta_right
  end interface convolute_kb_contour_gf
  public :: convolute_kb_contour_gf
  !
  !
  !OTHER ROUTINES && PLOT:
  public :: extrapolate_kb_contour_gf
  public :: save_kb_contour_gf
  public :: inquire_kb_contour_gf
  public :: read_kb_contour_gf
  public :: plot_kb_contour_gf
  !
  !
  !INTEGRATION ROUTINES:
  interface kb_trapz
     module procedure kb_trapz_d
     module procedure kb_trapz_c
  end interface kb_trapz
  !
  interface kb_half_trapz
     module procedure kb_half_trapz_d
     module procedure kb_half_trapz_c
  end interface kb_half_trapz
  !
  public :: kb_trapz,kb_half_trapz
  !
  !
  !VIE/VIDE SOLVER:
  interface vie_kb_contour_gf
     module procedure vie_kb_contour_gf_q
     module procedure vie_kb_contour_gf_delta
  end interface vie_kb_contour_gf
  public :: vie_kb_contour_gf
  public :: vide_kb_contour_gf
  !
  ! 
  !SYMMETRIES & COMPONENTS:
  public :: get_gtr
  public :: get_adv
  public :: get_rmix
  public :: get_bar
  !
  !
  !OVERLOAD OPERATORS
  interface operator(*)
     module procedure kb_contour_gf_scalarL_d
     module procedure kb_contour_gf_scalarL_c
     module procedure kb_contour_dgf_scalarL_d
     module procedure kb_contour_dgf_scalarL_c
     module procedure kb_contour_gf_scalarR_d
     module procedure kb_contour_gf_scalarR_c
     module procedure kb_contour_dgf_scalarR_d
     module procedure kb_contour_dgf_scalarR_c
  end interface operator(*)
  !
  interface assignment(=)
     module procedure kb_contour_gf_equality_
     module procedure kb_contour_dgf_equality_
     module procedure kb_contour_gf_equality__
     module procedure kb_contour_dgf_equality__
  end interface assignment(=)
  !
  !
  public :: operator(*)
  public :: assignment(=)



contains



  !======= ALLOCATE ======= 
  subroutine allocate_kb_contour_gf_main(G,params)
    type(kb_contour_gf)     :: G
    type(kb_contour_params) :: params
    integer                 :: i,j,N,L,Lf
    if(allocated(G%less))deallocate(G%less)
    if(allocated(G%ret)) deallocate(G%ret)
    if(allocated(G%lmix))deallocate(G%lmix)
    if(allocated(G%mats))deallocate(G%mats)
    N = params%Ntime            !<== allocate at maximum time
    L = params%Ntau
    Lf= params%Niw
    G%N = N
    G%L = L
    G%Lf= Lf
    allocate(G%less(N,N))  ; G%less=zero
    allocate(G%ret(N,N))   ; G%ret=zero
    allocate(G%lmix(N,0:L)); G%lmix=zero
    allocate(G%mats(0:L))  ; G%mats=0d0
    allocate(G%iw(Lf))     ; G%iw=zero
    G%status=.true.
  end subroutine allocate_kb_contour_gf_main
  !
  subroutine allocate_kb_contour_gf_nambu_redux(G,params)
    type(kb_contour_gf),dimension(2) :: G
    type(kb_contour_params)          :: params
    call allocate_kb_contour_gf_main(G(1),params)
    call allocate_kb_contour_gf_main(G(2),params)
    G(1)%anomalous=.false.
    G(2)%anomalous=.true.
  end subroutine allocate_kb_contour_gf_nambu_redux
  !
  subroutine allocate_kb_contour_gf_nambu(G,params)
    type(kb_contour_gf),dimension(2,2) :: G
    type(kb_contour_params)            :: params
    integer                            :: i,j
    do i=1,2
       do j=1,2
          call allocate_kb_contour_gf_main(G(i,j),params)
       enddo
    enddo
    G(1,1)%anomalous=.false.
    G(1,2)%anomalous=.true.
    G(2,1)%anomalous=.true.
    G(2,2)%anomalous=.false.
  end subroutine allocate_kb_contour_gf_nambu
  !
  !
  subroutine allocate_kb_contour_dgf_main(dG,params,wgtr)
    type(kb_contour_dgf)    :: dG
    type(kb_contour_params) :: params
    integer                 :: i,j,N,L
    logical,optional        :: wgtr
    if(allocated(dG%less))deallocate(dG%less)
    if(allocated(dG%ret)) deallocate(dG%ret)
    if(allocated(dG%lmix))deallocate(dG%lmix)
    N=params%Ntime           !<== allocate at maximum time
    L=params%Ntau
    dG%N=N
    dG%L=L
    allocate(dG%less(N))  ; dG%less=zero
    allocate(dG%ret(N))   ; dG%ret=zero
    allocate(dG%lmix(0:L)); dG%lmix=zero
    if(present(wgtr).AND.wgtr)then
       allocate(dG%gtr(N))
       dG%gtr=zero
    endif
    dG%status=.true.
  end subroutine allocate_kb_contour_dgf_main
  !
  subroutine allocate_kb_contour_dgf_nambu_redux(dG,params)
    type(kb_contour_dgf),dimension(2)     :: dG
    type(kb_contour_params) :: params
    call allocate_kb_contour_dgf_main(dG(1),params)
    call allocate_kb_contour_dgf_main(dG(2),params)
    dG(1)%anomalous=.false.
    dG(2)%anomalous=.true.
  end subroutine allocate_kb_contour_dgf_nambu_redux









  !======= DEALLOCATE ======= 
  subroutine deallocate_kb_contour_gf(G)
    type(kb_contour_gf) :: G
    if(.not.G%status)stop "contour_gf/deallocate_kb_contour_gf: G not allocated"
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
    if(allocated(dG%gtr))deallocate(dG%gtr)
    dG%N=0
    dG%L=0
    dG%status=.false.
  end subroutine deallocate_kb_contour_dgf






  !======= CHECK DIMENSION ======= 
  function check_dimension_kb_contour_gf(G,params) result(bool)
    type(kb_contour_gf)     :: G
    type(kb_contour_params) :: params
    logical                 :: bool
    integer                 :: N,L
    bool=.false.
    N=params%Ntime              !<== check size at maximum time
    L=params%Ntau
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
    if( size(G%lmix)/=N*(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
    if( size(G%mats)/=(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
    bool=.true.
  end function check_dimension_kb_contour_gf
  !
  function check_dimension_kb_contour_gf_(G,N,L) result(bool)
    type(kb_contour_gf)     :: G
    logical                 :: bool
    integer                 :: N,L
    bool=.false.
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         stop "ERROR contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
    if( size(G%lmix)/=N*(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
    if( size(G%mats)/=(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
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
         stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions less/ret"
    if( size(dG%lmix)/=(L+1) )&
         stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions lmix"
    bool=.true.
  end function check_dimension_kb_contour_dgf
  !
  function check_dimension_kb_contour_dgf_(dG,N,L) result(bool)
    type(kb_contour_dgf)    :: dG
    logical                 :: bool
    integer                 :: N,L
    bool=.false.
    if( (size(dG%less)/=N) .OR. (size(dG%ret)/=N) )&
         stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions less/ret"
    if( size(dG%lmix)/=(L+1) )&
         stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions lmix"
    bool=.true.
  end function check_dimension_kb_contour_dgf_
























  !======= SAVE ======= 
  subroutine save_kb_contour_gf(G,file)
    type(kb_contour_gf) :: G
    character(len=*)    :: file
    integer             :: unit
    if(.not.G%status)stop "contour_gf/save_kb_contour_gf: G not allocated"
    unit = free_unit()
    open(unit,file=reg(file)//"_dimension.data.neqipt")
    write(unit,"(A1,3A4,3X)")"#","N","L","Lf"
    write(unit,"(1X,3I4,3X)")G%N,G%L,G%Lf
    close(unit)
    call store_data(reg(file)//"_less.data.neqipt",G%less(:,:))
    call store_data(reg(file)//"_ret.data.neqipt", G%ret(:,:))
    call store_data(reg(file)//"_lmix.data.neqipt",G%lmix(:,0:))
    call store_data(reg(file)//"_mats.data.neqipt",G%mats(0:))
    call store_data(reg(file)//"_iw.data.neqipt",G%iw(:))
  end subroutine save_kb_contour_gf
  !
  subroutine save_kb_contour_dgf(dG,file)
    type(kb_contour_dgf)  :: dG
    character(len=*)      :: file
    integer :: unit
    if(.not.dG%status)stop "contour_gf/save_kb_contour_dgf: dG not allocated"
    unit = free_unit()
    open(unit,file=reg(file)//"_dimension.data.neqipt")
    write(unit,"(A1,2A4,2X)")"#","N","L"
    write(unit,"(1X,2I4,2X)")dG%N,dG%L
    close(unit)
    call store_data(reg(file)//"_less.data.neqipt",dG%less(:))
    call store_data(reg(file)//"_ret.data.neqipt", dG%ret(:))
    call store_data(reg(file)//"_lmix.data.neqipt",dG%lmix(0:))
  end subroutine save_kb_contour_dgf









  !======= READ ======= 
  subroutine read_kb_contour_gf(G,file)
    type(kb_contour_gf)  :: G
    character(len=*)     :: file
    logical              :: check
    check = inquire_kb_contour_gf(file)
    if(.not.G%status.OR..not.check)stop "contour_gf/read_kb_contour_gf: G not allocated"
    call read_data(trim(file)//"_less.data.neqipt",G%less(:,:))
    call read_data(trim(file)//"_ret.data.neqipt",G%ret(:,:))
    call read_data(trim(file)//"_lmix.data.neqipt",G%lmix(:,0:))
    call read_data(trim(file)//"_mats.data.neqipt",G%mats(0:))
    call read_data(trim(file)//"_iw.data.neqipt",G%iw(:))
  end subroutine read_kb_contour_gf
  !
  subroutine read_kb_contour_dgf(dG,file)
    type(kb_contour_dgf) :: dG
    character(len=*)     :: file
    logical              :: check
    check = inquire_kb_contour_dgf(file)
    if(.not.dG%status.OR..not.check)stop "contour_gf/read_kb_contour_dgf: dG not allocated"
    call read_data(trim(file)//"_less.data.neqipt",dG%less(:))
    call read_data(trim(file)//"_ret.data.neqipt",dG%ret(:))
    call read_data(trim(file)//"_lmix.data.neqipt",dG%lmix(0:))
  end subroutine read_kb_contour_dgf
  !
  function inquire_kb_contour_gf(file) result(check)
    integer          :: i
    logical          :: check,bool(5)
    character(len=*) :: file
    character(len=16),dimension(5)  :: ctype=([ character(len=5) :: 'less','ret','lmix','mats','iw'])
    check=.true.
    do i=1,5
       inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt",exist=bool(i))
       if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt.gz",exist=bool(i))
       check=check.AND.bool(i)
    enddo
  end function inquire_kb_contour_gf
  !
  function inquire_kb_contour_dgf(file) result(check)
    integer          :: i
    logical          :: check,bool(3)
    character(len=*) :: file
    character(len=16),dimension(3)  :: ctype=([character(len=5) :: 'less','ret','lmix'])
    check=.true.
    do i=1,3
       inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt",exist=bool(i))
       if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt.gz",exist=bool(i))
       check=check.AND.bool(i)
    enddo
  end function inquire_kb_contour_dgf










  !======= PLOT ======= 
  subroutine plot_kb_contour_gf(file,G,params)
    character(len=*)        :: file
    type(kb_contour_gf)     :: G
    type(kb_contour_params) :: params
    integer                 :: Nt,i
    if(.not.G%status)stop "contour_gf/plot_kb_contour_gf: G is not allocated" 
    Nt=params%Ntime
    call splot3d(reg(file)//"_less_t_t.data.neqipt",params%t(:Nt),params%t(:Nt),G%less(:Nt,:Nt))
    call splot3d(reg(file)//"_ret_t_t.data.neqipt",params%t(:Nt),params%t(:Nt),G%ret(:Nt,:Nt))
    call splot3d(reg(file)//"_lmix_t_tau.data.neqipt",params%t(:Nt),params%tau(0:),G%lmix(:Nt,0:))
    call splot(reg(file)//"_mats_tau.data.neqipt",params%tau(0:),G%mats(0:))
    call splot(reg(file)//"_tau_tau.data.neqipt",(/(i*params%beta/params%Niw,i=0,params%Niw)/),G%mats(0:))
    call splot(reg(file)//"_mats_iw.data.neqipt",params%wm(:),G%iw(:))
  end subroutine plot_kb_contour_gf
  !
  subroutine plot_kb_contour_dgf(file,dG,params)
    character(len=*)        :: file
    type(kb_contour_dgf)    :: dG
    type(kb_contour_params) :: params
    integer                 :: Nt
    if(.not.dG%status)stop "contour_gf/plot_kb_contour_gf: G is not allocated" 
    Nt=params%Ntime
    call splot(reg(file)//"_less_t.data.neqipt",params%t(:Nt),dG%less(:Nt))
    call splot(reg(file)//"_ret_t.data.neqipt",params%t(:Nt),dG%ret(:Nt))
    call splot(reg(file)//"_lmix_tau.data.neqipt",params%tau(0:),dG%lmix(0:))
  end subroutine plot_kb_contour_dgf






  !======= ADD ======= 
  !C(t,t')=A(t,t') + B(t,t'), with t=t_max && t'=0,t_max
  !t_max_index==N
  subroutine add_kb_contour_gf_simple(A,B,C,params)
    type(kb_contour_gf)     :: A,B,C
    type(kb_contour_params) :: params
    integer                 :: N,L
    logical                 :: checkA,checkB,checkC
    if(  (.not.A%status).OR.&
         (.not.B%status).OR.&
         (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
    N   = params%Ntime   !<== work with the TOTAL size of the contour
    L   = params%Ntau
    !
    checkA=check_dimension_kb_contour(A,N,L) 
    checkB=check_dimension_kb_contour(B,N,L)
    checkC=check_dimension_kb_contour(C,N,L)
    !
    C%less(:,:) = A%less(:,:) + B%less(:,:)
    C%ret(:,:)  = A%ret(:,:)  + B%ret(:,:)
    C%lmix(:,0:)= A%lmix(:,0:)+ B%lmix(:,0:)
    C%mats(0:)  = A%mats(0:)  + B%mats(0:)
    C%iw(:)     = A%iw(:)      + B%iw(:)
  end subroutine  add_kb_contour_gf_simple

  subroutine add_kb_contour_gf_recursive(A,C,params)
    type(kb_contour_gf)     :: A(:)
    type(kb_contour_gf)     :: C
    type(kb_contour_params) :: params
    integer                 :: i,Na,N,L
    logical                 :: checkA,checkC
    !
    Na=size(A)
    do i=1,Na
       if(  (.not.A(i)%status) )stop "contour_gf/add_kb_contour_gf: A(i) not allocated"
    enddo
    if(  (.not.C%status) )stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
    N   = params%Ntime    !<== work with the TOTAL size of the contour
    L   = params%Ntau
    !
    do i=1,Na
       checkA=check_dimension_kb_contour(A(i),N,L) 
    enddo
    checkC=check_dimension_kb_contour(C,N,L)
    !
    C=zero
    do i=1,Na
       C%less(:,:) = C%less(:,:) + A(i)%less(:,:)
       C%ret(:,:)  = C%ret(:,:)  + A(i)%ret(:,:)
       C%lmix(:,0:)= C%lmix(:,0:)+ A(i)%lmix(:,0:)
       C%mats(0:)  = C%mats(0:)  + A(i)%mats(0:)
       C%iw(:)     = C%iw(:)     + A(i)%iw(:)
    enddo
  end subroutine  add_kb_contour_gf_recursive

  subroutine add_kb_contour_dgf(A,B,C,params)
    type(kb_contour_dgf)    :: A,B,C
    type(kb_contour_params) :: params
    integer                 :: N,L
    logical                 :: checkA,checkB,checkC
    if(  (.not.A%status).OR.&
         (.not.B%status).OR.&
         (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    checkA=check_dimension_kb_contour(A,N,L)
    checkB=check_dimension_kb_contour(B,N,L)
    checkC=check_dimension_kb_contour(C,N,L)
    C%less(:) = A%less(:) + B%less(:)
    C%ret(:)  = A%ret(:)  + B%ret(:)  
    C%lmix(0:)= A%lmix(0:)+ B%lmix(0:)
  end subroutine add_kb_contour_dgf

  subroutine add_kb_contour_gf_delta_d(A,B,C,params)
    type(kb_contour_gf)     :: A,C
    real(8),dimension(:)    :: B
    type(kb_contour_params) :: params
    integer                 :: N,L,i
    logical                 :: checkA,checkC
    if(  (.not.A%status).OR.&
         (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
    N   = params%Ntime
    L   = params%Ntau
    !
    checkA=check_dimension_kb_contour(A,N,L) 
    checkC=check_dimension_kb_contour(C,N,L)
    if(size(B)<N)stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
    C%less(:,:) = A%less(:,:)
    C%ret(:,:)  = A%ret(:,:)
    C%lmix(:,0:)= A%lmix(:,0:)
    C%mats(0:)  = A%mats(0:)
    !C%iw(:)     = A%iw(:) + B(1)
    do i=1,N
       C%ret(i,i) = C%ret(i,i) + B(i)
    enddo
    C%mats(0) = C%mats(0) + B(1)
  end subroutine  add_kb_contour_gf_delta_d

  subroutine add_kb_contour_gf_delta_c(A,B,C,params)
    type(kb_contour_gf)     :: A,C
    complex(8),dimension(:) :: B
    type(kb_contour_params) :: params
    integer                 :: N,L,i
    logical                 :: checkA,checkC
    if(  (.not.A%status).OR.&
         (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
    N   = params%Ntime
    L   = params%Ntau
    !
    checkA=check_dimension_kb_contour(A,N,L) 
    checkC=check_dimension_kb_contour(C,N,L)
    if(size(B)<N)stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
    C%less(:,:) = A%less(:,:)
    C%ret(:,:)  = A%ret(:,:)
    C%lmix(:,0:)= A%lmix(:,0:)
    C%mats(0:)  = A%mats(0:)
    !C%iw(:)     = A%iw(:) + B(1)
    do i=1,N
       C%ret(i,i) = C%ret(i,i) + B(i)
    enddo
    C%mats(0) = C%mats(0) + B(1)
  end subroutine  add_kb_contour_gf_delta_c











  !======= SUM ======= 
  ! performs the sum a*A + b*B and stores it in C along the perimeter
  ! can be called multiple times to add up. 
  subroutine sum_kb_contour_gf_simple(A,ak,B,bk,C,params)
    type(kb_contour_gf)               :: A,B
    type(kb_contour_gf),intent(inout) :: C
    real(8)                           :: ak,bk
    type(kb_contour_params)           :: params
    integer                           :: N,L
    if(  (.not.A%status).OR.&
         (.not.B%status).OR.&
         (.not.C%status))stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
    !
    N   = params%Nt   !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    if(N==1)then
       C%mats(0:) = ak*A%mats(0:) + bk*B%mats(0:)
       C%iw(:)    = ak*A%iw(:)    + bk*B%iw(:)
    endif
    !
    C%ret(N,1:N)   = ak*A%ret(N,1:N)   + bk*B%ret(N,1:N)
    C%less(N,1:N)  = ak*A%less(N,1:N)  + bk*B%less(N,1:N)
    C%lmix(N,0:)   = ak*A%lmix(N,0:)   + bk*B%lmix(N,0:)
    !
    if(.not.C%anomalous)then
       C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
    else
       C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
    endif
    !
  end subroutine sum_kb_contour_gf_simple

  subroutine sum_kb_contour_gf_recursive(A,ak,C,params,iaddup)
    type(kb_contour_gf)               :: A(:)
    real(8)                           :: ak(size(A))
    type(kb_contour_gf),intent(inout) :: C
    type(kb_contour_params)           :: params
    integer                           :: N,L,Na,i
    logical,optional                  :: iaddup
    logical                           :: iaddup_
    !
    iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
    !
    Na=size(A)
    !
    do i=1,Na
       if((.not.A(i)%status) )stop "contour_gf/sum_kb_contour_gf: G or Gk not allocated"
    enddo
    if( (.not.C%status) )stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
    !
    N   = params%Nt   !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    if(.not.iaddup_)call del_kb_contour_gf(C,params)
    !
    if(N==1)then
       do i=1,Na
          C%mats(0:) = C%mats(0:) + ak(i)*A(i)%mats(0:)
          C%iw(:)    = C%iw(:)    + ak(i)*A(i)%iw(:)
       enddo
    endif
    !
    do i=1,Na
       C%ret(N,1:N)   = C%ret(N,1:N)  + ak(i)*A(i)%ret(N,1:N)
       C%less(N,1:N)  = C%less(N,1:N) + ak(i)*A(i)%less(N,1:N)
       C%lmix(N,0:)   = C%lmix(N,0:)  + ak(i)*A(i)%lmix(N,0:)
    enddo
    !
    if(.not.C%anomalous)then
       C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
    else
       C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
    endif
    !
  end subroutine sum_kb_contour_gf_recursive

  subroutine sum_kb_contour_gf_delta_d(A,ak,B,bk,C,params)
    type(kb_contour_gf)               :: A
    real(8),dimension(:)              :: B
    type(kb_contour_gf),intent(inout) :: C
    real(8)                           :: ak,bk
    type(kb_contour_params)           :: params
    integer                           :: N,L
    if(  (.not.A%status).OR.&
         (.not.C%status))stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
    !
    N   = params%Nt   !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    call del_kb_contour_gf(C,params)
    !
    if(N==1)then 
       C%mats(0)  = ak*A%mats(0)  + B(1)
       C%mats(1:) = ak*A%mats(1:)
    endif
    !
    C%ret(N,1:N)   = ak*A%ret(N,1:N)
    C%less(N,1:N)  = ak*A%less(N,1:N)
    C%lmix(N,0:)   = ak*A%lmix(N,0:)
    !
    if(.not.C%anomalous)then
       C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
    else
       C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
    endif
    !
    C%ret(N,N) = C%ret(N,N) + B(N)
    !
  end subroutine sum_kb_contour_gf_delta_d

  subroutine sum_kb_contour_gf_delta_c(A,ak,B,bk,C,params)
    type(kb_contour_gf)               :: A
    complex(8),dimension(:)           :: B
    type(kb_contour_gf),intent(inout) :: C
    real(8)                           :: ak,bk
    type(kb_contour_params)           :: params
    integer                           :: N,L

    if(  (.not.A%status).OR.&
         (.not.C%status))stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
    !
    N   = params%Nt   !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    call del_kb_contour_gf(C,params)
    !
    if(N==1)then 
       C%mats(0)  = ak*A%mats(0)  + B(1)
       C%mats(1:) = ak*A%mats(1:)
    endif
    !
    C%ret(N,1:N)   = ak*A%ret(N,1:N) 
    C%less(N,1:N)  = ak*A%less(N,1:N)
    C%lmix(N,0:)   = ak*A%lmix(N,0:)
    !
    if(.not.C%anomalous)then
       C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
    else
       C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
    endif
    !
    C%ret(N,N) = C%ret(N,N) + B(N)
    ! C%mats(0) = C%mats(0) + B(1)
    !
  end subroutine sum_kb_contour_gf_delta_c







  !======= DEL ======= 
  ! reset the function along the actual perimeter to zero:
  subroutine del_kb_contour_gf(G,params)
    type(kb_contour_gf)     :: G
    type(kb_contour_params) :: params
    integer                 :: N,L
    logical                 :: checkG
    if(  (.not.G%status)) stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
    !
    N   = params%Nt   !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    if(N==1)then
       G%iw = zero
       G%mats = 0d0
    endif
    !
    G%ret(N,1:N)   = zero
    G%less(N,1:N)  = zero
    G%lmix(N,0:)   = zero
  end subroutine del_kb_contour_gf

  subroutine del_kb_contour_dgf(dG,params)
    type(kb_contour_dgf)    :: dG
    type(kb_contour_params) :: params
    integer                 :: N,L
    logical                 :: checkA,checkB,checkC
    if(  (.not.dG%status))stop "contour_gf/del_kb_contour_dgf: dG not allocated"
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    dG%less(1:N) = zero
    dG%ret(1:N)  = zero
    dG%lmix(0:)  = zero
  end subroutine del_kb_contour_dgf


























  !======= EXTRAPOLATION ======= 
  !extrapolate a function from a given time
  !to the next one:
  subroutine extrapolate_kb_contour_gf(g,params)
    type(kb_contour_gf)     :: g
    type(kb_contour_params) :: params
    integer                 :: i,j,k,N,L
    if(.not.g%status)     stop "extrapolate_kb_contour_gf: g is not allocated"
    if(.not.params%status)stop "extrapolate_kb_contour_gf: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    !
    select case(N)
    case(1)
       return
    case(2)
       !GUESS G AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
       do j=1,N
          g%ret(N,j) =g%ret(1,1)
          g%less(N,j)=g%less(1,1)
       end do
       do i=1,N-1
          g%less(i,N)=g%less(1,1)
       end do
       do j=0,L
          g%lmix(N,j)=g%lmix(1,j)
       end do
    case default
       !EXTEND G FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
       !USING QUADRATIC EXTRAPOLATION
       do k=1,N-1
          g%less(N,k)=2.d0*g%less(N-1,k)-g%less(N-2,k)
          g%less(k,N)=2.d0*g%less(k,N-1)-g%less(k,N-2)
       end do
       g%less(N,N)=2.d0*g%less(N-1,N-1)-g%less(N-2,N-2)
       !
       do k=0,L
          g%lmix(N,k)=2.d0*g%lmix(N-1,k)-g%lmix(N-2,k)
       end do
       !
       g%ret(N,N)=-xi
       do k=1,N-2
          g%ret(N,k)=2.d0*g%ret(N-1,k)-g%ret(N-2,k)
       end do
       g%ret(N,N-1)=0.5d0*(g%ret(N,N)+g%ret(N,N-2))
    end select
  end subroutine extrapolate_kb_contour_gf











  !======= CONVOLUTION ======= 
  !C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max
  !t_max_index==N
  subroutine convolute_kb_contour_gf_simple(A,B,C,params,dcoeff,ccoeff)
    type(kb_contour_gf)                 :: A,B,C
    type(kb_contour_params)             :: params
    real(8),optional                    :: dcoeff
    complex(8),optional                 :: ccoeff
    integer                             :: N,L
    real(8)                             :: dt,dtau
    complex(8),dimension(:),allocatable :: AxB    
    integer                             :: i,j,k,itau,jtau
    logical                             :: checkA,checkB,checkC
    if(  (.not.A%status).OR.&
         (.not.B%status).OR.&
         (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: A,B,C not allocated"
    N   = params%Nt      !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    dt  = params%dt
    dtau= params%dtau
    !
    checkA=check_dimension_kb_contour(A,params%Ntime,L) 
    checkB=check_dimension_kb_contour(B,params%Ntime,L)
    checkC=check_dimension_kb_contour(C,params%Ntime,L)
    !
    allocate(AxB(0:max(L,N)))
    !
    if(N==1) then
       ! Matsubara frequencies
       C%iw=A%iw*B%iw
       ! Matsubara times
       AxB(0:) = zero
       do jtau=0,L
          do k=0,jtau
             AxB(k)=A%mats(k)*B%mats(L+k-jtau)
          end do
          C%mats(jtau)=C%mats(jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
          do k=jtau,L
             AxB(k)=A%mats(k)*B%mats(k-jtau)
          end do
          C%mats(jtau)=C%mats(jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
       enddo
    endif
    !
    !Ret. component
    !C^R(t,t')=\int_{t'}^t ds A^R(t,s)*B^R(s,t')
    C%ret(N,1:N)=zero
    do j=1,N       !for all t' in 0:t_max
       AxB(0:)  = zero !AxB should be set to zero before integration
       do k=j,N    !store the convolution between t'{=j} and t{=N}
          AxB(k) = A%ret(N,k)*B%ret(k,j)
       enddo
       C%ret(n,j) = C%ret(n,j) + dt*kb_trapz(AxB(0:),j,N)
    enddo
    !
    !Lmix. component
    !C^\lmix(t,tau')=\int_0^{beta} ds A^\lmix(t,s)*B^M(s,tau')
    !                  +\int_0^{t} ds A^R(t,s)*B^\lmix(s,tau')
    !               = I1 + I2
    C%lmix(N,0:L)=zero
    do jtau=0,L
       !I1:
       AxB(0:) = zero
       !break the integral I1 in two parts to take care of the 
       !sign of (tau-tau').
       do k=0,jtau
          AxB(k)=A%lmix(N,k)*B%mats(L+k-jtau)
       end do
       C%lmix(N,jtau)=C%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
       do k=jtau,L
          AxB(k)=A%lmix(N,k)*B%mats(k-jtau)
       end do
       C%lmix(n,jtau)=C%lmix(n,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
       !
       !I2:
       AxB(0:) = zero
       do k=1,N
          AxB(k) = A%ret(N,k)*B%lmix(k,jtau)
       enddo
       C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
    enddo
    !
    !Less component
    !C^<(t,t')=-i\int_0^{beta} ds A^\lmix(t,s)*B^\rmix(s,t')
    !             +\int_0^{t'} ds A^<(t,s)*B^A(s,t')
    !             +\int_0^{t} ds A^R(t,s)*B^<(s,t')
    ! (t,t')=>(N,j) <==> Vertical side, with tip (j=1,N) !!no tip (j=1,N-1)
    do j=1,N!-1
       C%less(N,j)=zero
       do k=0,L
          !         AxB(k)=A%lmix(N,k)*conjg(B%lmix(j,L-k))
          AxB(k)=A%lmix(N,k)*get_rmix(B,k,j,L)
       end do
       C%less(N,j)=C%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
       !
       do k=1,j
          !         AxB(k)=A%less(N,k)*conjg(B%ret(j,k))
          AxB(k)=A%less(N,k)*get_adv(B,k,j)
       end do
       C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
       !
       do k=1,N
          AxB(k)=A%ret(N,k)*B%less(k,j)
       end do
       C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
    end do
    !
    ! (t,t`)=>(i,N) <==> Horizontal side, w/ no tip (i=1,N-1)
    do i=1,N-1
       C%less(i,N)=zero
    enddo
    ! do i=1,N
    !    C%less(i,N)=zero
    !    do k=0,L
    !       !         AxB(k)=A%lmix(i,k)*conjg(B%lmix(n,L-k))
    !       AxB(k)=A%lmix(i,k)*get_rmix(B,k,N,L)
    !    end do
    !    C%less(i,N)=C%less(i,N)-xi*dtau*kb_trapz(AxB(0:),0,L)
    !    !
    !    do k=1,N
    !       !         AxB(k)=A%less(i,k)*conjg(B%ret(N,k))
    !       AxB(k)=A%less(i,k)*get_adv(B,k,N)
    !    end do
    !    C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,N)
    !    !
    !    do k=1,i
    !       AxB(k)=A%ret(i,k)*B%less(k,N)
    !    end do
    !    C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,i)
    ! end do
    !    
    if(present(dcoeff))then
       C%lmix(N,0:L) = dcoeff*C%lmix(N,0:L)
       C%less(N,1:N-1)=dcoeff*C%less(N,1:N-1)
       C%less(1:N,N)=dcoeff*C%less(1:N,N)
       C%ret(N,1:N)=dcoeff*C%ret(N,1:N)
    endif
    if(present(ccoeff))then
       C%lmix(N,0:L) = ccoeff*C%lmix(N,0:L)
       C%less(N,1:N-1)=ccoeff*C%less(N,1:N-1)
       C%less(1:N,N)=ccoeff*C%less(1:N,N)
       C%ret(N,1:N)=ccoeff*C%ret(N,1:N)
    endif
    deallocate(AxB)
  end subroutine convolute_kb_contour_gf_simple


  subroutine convolute_kb_contour_gf_recursive(A,C,params,dcoeff,ccoeff)
    type(kb_contour_gf)                 :: A(:)
    type(kb_contour_gf)                 :: C
    type(kb_contour_gf)                 :: Knew,Kold
    type(kb_contour_params)             :: params
    real(8),optional                    :: dcoeff
    complex(8),optional                 :: ccoeff
    integer                             :: N,L
    integer                             :: i,j,k,itau,jtau,Na
    logical                             :: checkA,checkC
    !
    Na = size(A)
    !
    if(  (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: C not allocated"
    do i=1,Na
       if(  (.not.A(i)%status))stop "contour_gf/convolute_kb_contour_gf: A(i) not allocated"
    enddo
    !
    call allocate_kb_contour_gf(Knew,params)
    call allocate_kb_contour_gf(Kold,params)
    !
    N   = params%Nt      !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    checkC=check_dimension_kb_contour(C,params%Ntime,L)
    do i=1,Na
       checkA=check_dimension_kb_contour(A(i),params%Ntime,L) 
    enddo
    !
    !Perform the recursive convolution:
    Kold = A(1)
    do i=2,Na-1
       call convolute_kb_contour_gf(Kold,A(i),Knew,params)
    enddo
    call convolute_kb_contour_gf(Knew,A(Na),C,params)
    !
    if(present(dcoeff))then
       C%lmix(N,0:L) = dcoeff*C%lmix(N,0:L)
       C%less(N,1:N-1)=dcoeff*C%less(N,1:N-1)
       C%less(1:N,N)=dcoeff*C%less(1:N,N)
       C%ret(N,1:N)=dcoeff*C%ret(N,1:N)
    endif
    if(present(ccoeff))then
       C%lmix(N,0:L) = ccoeff*C%lmix(N,0:L)
       C%less(N,1:N-1)=ccoeff*C%less(N,1:N-1)
       C%less(1:N,N)=ccoeff*C%less(1:N,N)
       C%ret(N,1:N)=ccoeff*C%ret(N,1:N)
    endif
    call deallocate_kb_contour_gf(Knew)
    call deallocate_kb_contour_gf(Kold)
    !    
  end subroutine convolute_kb_contour_gf_recursive


  !======= CONVOLUTION WITH A FUNCTION MULTIPLIED BY THE DELTA FUNCTION======= 
  !C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max
  !t_max_index==N
  ! The case A(t,t')=A(t)delta(t,t')
  subroutine convolute_kb_contour_gf_delta_left(A,B,C,params)
    type(kb_contour_gf)                 :: B,C
    complex(8),dimension(:)             :: A
    type(kb_contour_params)             :: params
    integer                             :: N,L
    integer                             :: i,j,k,itau,jtau
    logical                             :: checkB,checkC
    if(  (.not.B%status).OR.&
         (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: A,B,C not allocated"
    N   = params%Nt      !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    checkB=check_dimension_kb_contour(B,params%Ntime,L) 
    checkC=check_dimension_kb_contour(C,params%Ntime,L)
    !
    if (N.eq.1) then
       ! Matsubara
       C%mats = A(1)*B%mats
    endif
    !
    do j=1,N
       !Ret. component
       C%ret(n,j) = A(n)*B%ret(n,j)
       !Less. component
       C%less(n,j) = A(n)*B%less(n,j)
    enddo
    ! do i=1,N-1
    !    !Less. component
    !    C%less(i,n)=A%less(i,N)*B(N)
    ! enddo
    !
    !Lmix. component
    do jtau=0,L
       C%lmix(n,jtau) = A(n)*B%lmix(n,jtau)
    enddo
  end subroutine convolute_kb_contour_gf_delta_left
  ! The case B(t,t')=B(t)delta(t,t')
  subroutine convolute_kb_contour_gf_delta_right(A,B,C,params)
    type(kb_contour_gf)                 :: A,C
    complex(8),dimension(:)             :: B
    type(kb_contour_params)             :: params
    integer                             :: N,L
    integer                             :: i,j,k,itau,jtau
    logical                             :: checkA,checkC
    !
    if(  (.not.A%status).OR.&
         (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: A,B,C not allocated"
    N   = params%Nt      !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    checkA=check_dimension_kb_contour(A,params%Ntime,L) 
    checkC=check_dimension_kb_contour(C,params%Ntime,L)
    !
    if (N.eq.1) then
       ! Matsubara
       C%mats(0:) = A%mats(0:)*B(1)
    endif
    !
    do j=1,N
       !Ret. component
       C%ret(n,j) = A%ret(n,j)*B(j)
       !Less. component
       C%less(n,j) = A%less(n,j)*B(j)
    enddo
    ! do i=1,N-1
    !    !Less. component
    !    C%less(i,n)=A%less(i,N)*B(N)
    ! enddo
    !
    !Lmix. component
    do jtau=0,L
       C%lmix(N,jtau) = A%lmix(N,jtau)*B(1)
    enddo
  end subroutine convolute_kb_contour_gf_delta_right













  !======= VOLTERRA INTEGRAL EQUATION ======= 
  !----------------------------------------------------------------------------
  !  This subroutine solves a Volterra integral equation of the second kind,
  !              G(t,t') = Q(t,t') + (K*G)(t,t')
  !  for t=n*dt or t'=n*dt, using 2^nd *implicit* Runge-Kutta method.
  !----------------------------------------------------------------------------
  subroutine vie_kb_contour_gf_q(Q,K,G,params)
    type(kb_contour_gf),intent(in)      :: Q
    type(kb_contour_gf),intent(in)      :: K
    type(kb_contour_gf),intent(inout)   :: G
    type(kb_contour_params),intent(in)  :: params
    integer                             :: N,L
    real(8)                             :: dt,dtau
    integer                             :: i,j,s,itau,jtau
    complex(8),dimension(:),allocatable :: KxG
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    dt  = params%dt
    dtau= params%dtau
    !
    allocate(KxG(0:max(N,L)))
    !
    ! Ret component
    ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
    G%ret(N,N)=Q%ret(N,N)
    do j=1,N-1
       G%ret(N,j)=Q%ret(N,j)
       do s=j,N-1
          KxG(s)=K%ret(N,s)*G%ret(s,j)
       eNd do
       G%ret(N,j)=G%ret(N,j) + dt*kb_half_trapz(KxG(0:),j,N-1)
       G%ret(N,j)=G%ret(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
    end do
    !
    ! Lmix component
    ! G^\lmix(t,tau') - \int_0^t K^R(t,s)*G^\lmix(s,tau')ds
    !    = Q^\lmix(t,tau') + \int_0^\beta K^\lmix(t,s)*G^M(s,tau')ds
    do jtau=0,L
       G%lmix(N,jtau)=Q%lmix(N,jtau)
       do s=0,jtau
          KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau)-dtau*kb_trapz(KxG(0:),0,jtau)
       do s=jtau,L
          KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau)+dtau*kb_trapz(KxG(0:),jtau,L)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau) + dt*kb_half_trapz(KxG(0:),1,N-1)
       G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0-0.5d0*dt*K%ret(N,N))
    end do
    !
    ! Less component
    ! G^<(t,t') - \int_0^t K^{R}(t,s)*G^{<}(s,t')ds
    !    = Q^<(t,t') - i\int_0^\beta K^\lmix(t,s)*G^\rmix(s,t')ds
    !      + \int_0^{t'} K^<(t,s)*G^A(s,t')ds
    ! G^<(t_{N},t_{j})
    do j=1, N-1
       G%less(N,j)=Q%less(N,j)
       do s=0,L
          !         KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s))
          KxG(s)=K%lmix(N,s)*get_rmix(G,s,j,L)
       enddo
       G%less(N,j)=G%less(N,j)-xi*dtau*kb_trapz(KxG(0:),0,L)
       !
       do s=1,j
          !         KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
          KxG(s)=K%less(N,s)*get_adv(G,s,j)
       enddo
       G%less(N,j)=G%less(N,j)+dt*kb_trapz(KxG(0:),1,j)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%less(s,j)
       enddo
       G%less(N,j)=G%less(N,j) + dt*kb_half_trapz(KxG(0:),1,N-1)
       G%less(N,j)=G%less(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
    end do
    !
    ! G^<(t_{i},t_{N}) <= Hermite conjugate
    if (.not.G%anomalous) then
       do i=1,N-1
          G%less(i,N) = -conjg(G%less(N,i))
       end do
    else
       do i=1,N-1
          G%less(i,N) = G%less(N,i)+G%ret(N,i)
       end do
    endif
    !
    ! G^{<}(t_{n},t_{n}) <= Diagonal
    G%less(N,N)=Q%less(N,N)
    do s=0,L
       !      KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s))
       KxG(s)=K%lmix(N,s)*get_rmix(G,s,N,L)
    end do
    G%less(N,N)=G%less(N,N)-xi*dtau*kb_trapz(KxG(0:),0,L)
    !
    do s=1,N
       !      KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
       KxG(s)=K%less(N,s)*get_adv(G,s,N)
    end do
    G%less(N,N)=G%less(N,N)+dt*kb_trapz(KxG(0:),1,N)
    !
    do s=1,N-1
       KxG(s)=K%ret(N,s)*G%less(s,N)
    end do
    G%less(N,N)=G%less(N,N) + dt*kb_half_trapz(KxG(0:),1,N-1)
    G%less(N,N)=G%less(N,N)/(1.d0-0.5d0*dt*K%ret(N,N))
    !
    deallocate(KxG)
    !
  end subroutine vie_kb_contour_gf_q

  subroutine vie_kb_contour_gf_delta(K,G,params)
    type(kb_contour_gf),intent(in)      :: K
    type(kb_contour_gf),intent(inout)   :: G
    type(kb_contour_params),intent(in)  :: params
    integer                             :: N,L
    real(8)                             :: dt,dtau
    integer                             :: i,j,s,itau,jtau
    complex(8),dimension(:),allocatable :: KxG
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    dt  = params%dt
    dtau= params%dtau
    !
    allocate(KxG(0:max(N,L)))
    !
    ! Ret component
    ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
    G%ret(N,N)=one!Q%ret(N,N)
    do j=1,N-1
       G%ret(N,j)=zero!Q%ret(N,j)
       do s=j,N-1
          KxG(s)=K%ret(N,s)*G%ret(s,j)
       eNd do
       G%ret(N,j)=G%ret(N,j) + dt*kb_half_trapz(KxG(0:),j,N-1)
       G%ret(N,j)=G%ret(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
    end do
    !
    ! Lmix component
    ! G^\lmix(t,tau') - \int_0^t K^R(t,s)*G^\lmix(s,tau')ds
    !    = Q^\lmix(t,tau') + \int_0^\beta K^\lmix(t,s)*G^M(s,tau')ds
    do jtau=0,L
       G%lmix(N,jtau)=zero!Q%lmix(N,jtau)
       do s=0,jtau
          KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau)-dtau*kb_trapz(KxG(0:),0,jtau)
       do s=jtau,L
          KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau)+dtau*kb_trapz(KxG(0:),jtau,L)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
       end do
       G%lmix(N,jtau)=G%lmix(N,jtau) + dt*kb_half_trapz(KxG(0:),1,N-1)
       G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0-0.5d0*dt*K%ret(N,N))
    end do
    !
    ! Less component
    ! G^<(t,t') - \int_0^t K^{R}(t,s)*G^{<}(s,t')ds
    !    = Q^<(t,t') - i\int_0^\beta K^\lmix(t,s)*G^\rmix(s,t')ds
    !      + \int_0^{t'} K^<(t,s)*G^A(s,t')ds
    ! G^<(t_{N},t_{j})
    do j=1, N-1
       G%less(N,j)=zero!Q%less(N,j)
       do s=0,L
          !         KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s))
          KxG(s)=K%lmix(N,s)*get_rmix(G,s,j,L)
       enddo
       G%less(N,j)=G%less(N,j)-xi*dtau*kb_trapz(KxG(0:),0,L)
       !
       do s=1,j
          !         KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
          KxG(s)=K%less(N,s)*get_adv(G,s,j)
       enddo
       G%less(N,j)=G%less(N,j)+dt*kb_trapz(KxG(0:),1,j)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%less(s,j)
       enddo
       G%less(N,j)=G%less(N,j) + dt*kb_half_trapz(KxG(0:),1,N-1)
       G%less(N,j)=G%less(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
    end do
    !
    ! G^<(t_{i},t_{N}) <= Hermite conjugate
    if (.not.G%anomalous) then
       do i=1,N-1
          G%less(i,N) = -conjg(G%less(N,i))
       end do
    else
       do i=1,N-1
          G%less(i,N) = G%less(N,i)+G%ret(N,i)
       end do
    endif
    !
    ! G^{<}(t_{n},t_{n}) <= Diagonal
    G%less(N,N)=zero!Q%less(N,N)
    do s=0,L
       !      KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s))
       KxG(s)=K%lmix(N,s)*get_rmix(G,s,N,L)
    end do
    G%less(N,N)=G%less(N,N)-xi*dtau*kb_trapz(KxG(0:),0,L)
    !
    do s=1,N
       !      KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
       KxG(s)=K%less(N,s)*get_adv(G,s,N)
    end do
    G%less(N,N)=G%less(N,N)+dt*kb_trapz(KxG(0:),1,N)
    !
    do s=1,N-1
       KxG(s)=K%ret(N,s)*G%less(s,N)
    end do
    G%less(N,N)=G%less(N,N) + dt*kb_half_trapz(KxG(0:),1,N-1)
    G%less(N,N)=G%less(N,N)/(1.d0-0.5d0*dt*K%ret(N,N))
    !
    deallocate(KxG)
    !
  end subroutine vie_kb_contour_gf_delta








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
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    dt  = params%dt
    dtau= params%dtau
    !call allocate_kb_contour_dgf(dG_new,params)
    !
    allocate(KxG(0:max(N,L)))
    !
    !TYPE: X=RET,LESS,LMIX (R,<,\LMIX)
    !I D/DT G^X(T,:)  = H(T)G^X(T,:) + Q^X(T,:) + INT_{0,T'}^{T}K^R(T,S)G^X(S,:)DS
    !
    !Ret component
    ! d/dt G^R(t,:) = -i*h(t)*G^R(t,:) -i*delta(t,:) - i\int_{:}^t K^R(t,s)*G^R(s,:)ds
    G%ret(N,N)=-xi
    dG_new%ret(N)=-xi*H(N)*G%ret(N,N)
    do j=1,N-1
       G%ret(N,j)=G%ret(N-1,j) + 0.5d0*dt*dG%ret(j)
       do s=j,N-1
          KxG(s)=K%ret(N,s)*G%ret(s,j)
       enddo
       dG_new%ret(j)=-xi*dt*kb_half_trapz(KxG(0:),j,N-1)
       !
       G%ret(N,j)=G%ret(N,j) + 0.5d0*dt*dG_new%ret(j)
       G%ret(N,j)=G%ret(N,j)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
       !
       dG_new%ret(j)=dG_new%ret(j) - xi*H(N)*G%ret(N,j) - 0.5d0*xi*dt*K%ret(N,N)*G%ret(N,j)
    enddo
    !
    !Lmix component
    !d/dt G^\lmix(t,:) = -i*H(t)*G^\lmix(t,:) -i\int_0^\beta K^\lmix(t,s)*G^M(s,:)ds -i\int_0^t K^R(t,s)G^\lmix(s,:)ds
    do jtau=0,L
       G%lmix(N,jtau)=G%lmix(N-1,jtau)+0.5d0*dt*dG%lmix(jtau)
       do s=0,jtau
          KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
       end do
       dG_new%lmix(jtau)=xi*dtau*kb_trapz(KxG(0:),0,jtau)
       do s=jtau,L
          KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
       end do
       dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*dtau*kb_trapz(KxG(0:),jtau,L)!<= add -iQ(t)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
       end do
       dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
       !
       G%lmix(N,jtau)=G%lmix(N,jtau) + 0.5d0*dt*dG_new%lmix(jtau)
       G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
       !
       dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*H(N)*G%lmix(N,jtau)-0.5d0*xi*dt*K%ret(N,N)*G%lmix(N,jtau)
    end do
    !
    !Less component
    !d/dt G^<(t,:) = -i*H(t)*G^<(t,:) -i*[ (-i)*\int_0^\beta K^\lmix(t,s)*G^\rmix(s,:)ds + \int_0^{:}K^<(t,s)*G^A(s,:)ds ]
    !                                 -i*\int_0^t K^R(t,s)*G^<(s,:)ds
    !
    ! G^<(t_{N},t_{j}), d/dt G^<(t_{N},t_{j}) <== lower-right triangle
    do j=1,N-1
       G%less(N,j)=G%less(N-1,j) + 0.5d0*dt*dG%less(j)
       do s=0,L
          !         KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s))
          KxG(s)=K%lmix(N,s)*get_rmix(G,s,j,L)
       end do
       dG_new%less(j)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
       do s=1,j
          !         KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
          KxG(s)=K%less(N,s)*get_adv(G,s,j)
       end do
       dG_new%less(j)=dG_new%less(j)-xi*dt*kb_trapz(KxG(0:),1,j)!<= -iQ(t)
       !
       do s=1,N-1
          KxG(s)=K%ret(N,s)*G%less(s,j)
       end do
       dG_new%less(j)=dG_new%less(j)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
       !
       G%less(N,j)=G%less(N,j) + 0.5d0*dt*dG_new%less(j)
       G%less(N,j)=G%less(N,j)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
       !
       dG_new%less(j)=dG_new%less(j)-xi*H(N)*G%less(N,j)-0.5d0*xi*dt*K%ret(N,N)*G%less(N,j)
    end do
    !
    ! G^<(t_{i},t_{N}), d/dt G^<(t_{i},t_{N}) <== upper left triangle
    ! Hermitian conjugate G
    if (.not.G%anomalous) then
       do i=1,N-1
          G%less(i,N)=-conjg(G%less(N,i))
       end do
    else
       do i=1,N-1
          G%less(i,N)=G%less(N,i)+G%ret(N,i)
       end do
    endif
    !
    ! d/dt G^<(t_{N-1},t_{N})
    dG_less=-xi*H(N-1)*G%less(N-1,N)
    do s=0,L
       !      KxG(s)=K%lmix(N-1,s)*conjg(G%lmix(N,L-s))
       KxG(s)=K%lmix(N-1,s)*get_rmix(G,s,N,L)
    end do
    dG_less=dG_less-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
    do s=1,N
       !      KxG(s)=K%less(N-1,s)*conjg(G%ret(N,s))
       KxG(s)=K%less(N-1,s)*get_adv(G,s,N)
    end do
    dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N)
    do s=1,N-1
       KxG(s)=K%ret(N-1,s)*G%less(s,N)
    end do
    dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N-1)
    !
    !G^<(N,N), d/dt G^<(N,N)
    G%less(N,N)=G%less(N-1,N)+0.5d0*dt*dG_less
    do s=0,L
       !      KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s))
       KxG(s)=K%lmix(N,s)*get_rmix(G,s,N,L)
    end do
    dG_new%less(N)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
    !
    do s=1,N
       !      KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
       KxG(s)=K%less(N,s)*get_adv(G,s,N)
    end do
    dG_new%less(N)=dG_new%less(N)-xi*dt*kb_trapz(KxG(0:),1,N)
    !
    do s=1,N-1
       KxG(s)=K%ret(N,s)*G%less(s,N)
    end do
    dG_new%less(N)=dG_new%less(N)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
    !
    G%less(N,N)=G%less(N,N)+0.5d0*dt*dG_new%less(N)
    G%less(N,N)=G%less(N,N)/(1.d0+0.5d0*xi*dt*H(N)+0.25d0*xi*dt**2*K%ret(N,N))
    !
    dG_new%less(N)=dG_new%less(N)-xi*H(N)*G%less(N,N)-0.5d0*xi*dt*K%ret(N,N)*G%less(N,N)
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
    real(8),dimension(0:),intent(in) :: f
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
    complex(8),dimension(0:),intent(in) :: f
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
    real(8),dimension(0:),intent(in) :: f
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
    complex(8),dimension(0:),intent(in) :: f
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
    G1%lmix(:,0:) = C
    G1%mats(0:)   = C
    G1%iw(:)      = C
  end subroutine kb_contour_gf_equality_
  !
  subroutine kb_contour_gf_equality__(G1,G2)
    type(kb_contour_gf),intent(inout) :: G1
    type(kb_contour_gf),intent(in)    :: G2
    G1%less(:,:)  = G2%less(:,:)
    G1%ret(:,:)   = G2%ret(:,:)
    G1%lmix(:,0:) = G2%lmix(:,0:)
    G1%mats(0:)   = G2%mats(0:)
    G1%iw(:)      = G2%iw(:)
  end subroutine kb_contour_gf_equality__
  !
  subroutine kb_contour_dgf_equality_(dG1,C)
    type(kb_contour_dgf),intent(inout) :: dG1
    complex(8),intent(in)             :: C
    dG1%less(:)  = C
    dG1%ret(:)   = C
    dG1%lmix(0:) = C
  end subroutine kb_contour_dgf_equality_
  !
  subroutine kb_contour_dgf_equality__(dG1,dG2)
    type(kb_contour_dgf),intent(inout) :: dG1
    type(kb_contour_dgf),intent(in)    :: dG2
    dG1%less(:)  = dG2%less(:)
    dG1%ret(:)   = dG2%ret(:)
    dG1%lmix(0:) = dG2%lmix(0:)
  end subroutine kb_contour_dgf_equality__



  function kb_contour_gf_scalarL_d(C,G) result(F)
    real(8),intent(in)             :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = C*G%less(:,:)
    F%ret(:,:)  = C*G%ret(:,:)
    F%lmix(:,0:)= C*G%lmix(:,0:)
    F%mats(0:)  = C*G%mats(0:)
    F%iw(:)     = C*G%iw(:)
  end function kb_contour_gf_scalarL_d
  !
  function kb_contour_dgf_scalarL_d(C,dG) result(dF)
    real(8),intent(in)             :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = C*dG%less(:)
    dF%ret(:)  = C*dG%ret(:)
    dF%lmix(0:)= C*dG%lmix(0:)
  end function kb_contour_dgf_scalarL_d



  function kb_contour_gf_scalarL_c(C,G) result(F)
    complex(8),intent(in)          :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = C*G%less(:,:)
    F%ret(:,:)  = C*G%ret(:,:)
    F%lmix(:,0:)= C*G%lmix(:,0:)
    F%mats(0:)  = C*G%mats(0:)
    F%iw(:)     = C*G%iw(:)
  end function kb_contour_gf_scalarL_c
  !
  function kb_contour_dgf_scalarL_c(C,dG) result(dF)
    complex(8),intent(in)           :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = C*dG%less(:)
    dF%ret(:)  = C*dG%ret(:)
    dF%lmix(0:)= C*dG%lmix(0:)
  end function kb_contour_dgf_scalarL_c





  function kb_contour_gf_scalarR_d(G,C) result(F)
    real(8),intent(in)             :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = G%less(:,:)*C
    F%ret(:,:)  = G%ret(:,:)*C
    F%lmix(:,0:)= G%lmix(:,0:)*C
    F%mats(0:)  = G%mats(0:)*C
    F%iw(:)     = G%iw(:)*C
  end function kb_contour_gf_scalarR_d
  !
  function kb_contour_dgf_scalarR_d(dG,C) result(dF)
    real(8),intent(in)             :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = dG%less(:)*C
    dF%ret(:)  = dG%ret(:)*C
    dF%lmix(0:)= dG%lmix(0:)*C
  end function kb_contour_dgf_scalarR_d




  function kb_contour_gf_scalarR_c(G,C) result(F)
    complex(8),intent(in)             :: C
    type(kb_contour_gf),intent(in) :: G
    type(kb_contour_gf)            :: F
    F%less(:,:) = G%less(:,:)*C
    F%ret(:,:)  = G%ret(:,:)*C
    F%lmix(:,0:)= G%lmix(:,0:)*C
    F%mats(0:)  = G%mats(0:)*C
    F%iw(:)     = G%iw(:)*C
  end function kb_contour_gf_scalarR_c
  !
  function kb_contour_dgf_scalarR_c(dG,C) result(dF)
    complex(8),intent(in)             :: C
    type(kb_contour_dgf),intent(in) :: dG
    type(kb_contour_dgf)            :: dF
    dF%less(:) = dG%less(:)*C
    dF%ret(:)  = dG%ret(:)*C
    dF%lmix(0:)= dG%lmix(0:)*C
  end function kb_contour_dgf_scalarR_c

  function get_adv(G,i,j) result (adv)
    implicit none
    type(kb_contour_gf),intent(in) :: G
    integer :: i,j
    complex(8) :: adv
    if (.not.G%anomalous) then
       adv=conjg(G%ret(j,i))
       return
    else
       adv=G%ret(j,i)
       return
    endif
  end function get_adv

  function get_rmix(G,i,j,L) result (rmix)
    implicit none
    type(kb_contour_gf),intent(in) :: G
    integer :: i,j,L
    complex(8) :: rmix
    if (.not.G%anomalous) then
       rmix=conjg(G%lmix(j,L-i))
       return
    else
       rmix=G%lmix(j,i)
       return
    endif
  end function get_rmix

  function get_gtr(G,i,j) result (gtr)
    implicit none
    type(kb_contour_gf),intent(in) :: G
    integer :: i,j
    complex(8) :: gtr
    if (.not.G%anomalous) then
       if (i.ge.j) then
          gtr=G%less(i,j)+G%ret(i,j)
          return
       else
          gtr=G%less(i,j)-conjg(G%ret(j,i))
          return
       endif
    else
       gtr=G%less(j,i)
       return
    endif
  end function get_gtr


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  comment
  !+-----------------------------------------------------------------------------+!
  subroutine get_bar(A,B,params)
    type(kb_contour_gf)                 :: A,B
    type(kb_contour_params)             :: params
    integer                             :: N,L,Lf
    real(8)                             :: dt,dtau
    integer                             :: i,j,k,itau,jtau
    logical                             :: checkA,checkB,checkC

    if(  (.not.A%status).OR.&
         (.not.B%status))stop "contour_gf/get_bar_kb_contour_gf: A,B,C not allocated"
    N   = params%Nt      !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    Lf  = params%Niw
    dt  = params%dt
    dtau= params%dtau
    !
    checkA=check_dimension_kb_contour(A,params%Ntime,L) 
    checkB=check_dimension_kb_contour(B,params%Ntime,L)
    !
    if(N==1)then
       ! Matsubara imaginary time component
       if (.not.B%anomalous) then
          A%mats(0:L) = B%mats(L:0:-1)
       else
          A%mats(0:L) = B%mats(0:L)
       endif
       !
       ! Matsubara frequencies component
       if (.not.B%anomalous) then
          A%iw = -conjg(B%iw)
       else
          A%iw = B%iw
       endif
    endif
    !
    !
    !Ret. component
    if (.not.B%anomalous) then
       do j=1,N
          A%ret(n,j) = -conjg(B%ret(n,j))
       enddo
    else
       do j=1,N
          A%ret(n,j) =  conjg(B%ret(n,j))
       enddo
    endif
    !
    !Lmix. component
    if (.not.B%anomalous) then
       do jtau=0,L
          A%lmix(N,jtau) = -conjg(B%lmix(N,L-jtau))
       enddo
    else
       do jtau=0,L
          A%lmix(N,jtau) =  B%lmix(N,L-jtau)
       enddo
    endif
    !
    !Less component
    if (.not.B%anomalous) then
       do j=1,N-1
          A%less(N,j)=-B%less(j,N)+conjg(B%ret(N,j))
       end do
    else
       do j=1,N-1
          A%less(N,j)=-conjg(B%less(j,N))
       end do
    endif
    !
    ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
    if (.not.B%anomalous) then
       do j=1,N-1
          A%less(j,N)=-B%less(N,j)-B%ret(N,j)
       end do
    else
       do j=1,N-1
          B%less(N,j)=-conjg(B%less(j,N))
       end do
    endif
  end subroutine get_bar
END MODULE NEQ_CONTOUR_GF

