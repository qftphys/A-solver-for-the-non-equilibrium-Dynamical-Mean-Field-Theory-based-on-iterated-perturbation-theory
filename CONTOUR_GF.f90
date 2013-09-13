MODULE CONTOUR_GF
  USE COMMON_VARS
  USE TOOLS
  USE IOTOOLS
  implicit none
  private


  !##################################################################
  ! KELDYSH CONTOUR GREEN'S FUNCTIONS
  !##################################################################
  type :: keldysh_contour_gf
     complex(8),dimension(:,:),pointer  :: less,ret
     logical                            :: status=.false.
     integer                            :: N=0
  end type keldysh_contour_gf
  public :: keldysh_contour_gf
  public :: allocate_keldysh_contour_gf
  public :: deallocate_keldysh_contour_gf
  public :: write_keldysh_contour_gf
  public :: read_keldysh_contour_gf
  public :: plot_keldysh_contour_gf
  public :: inquire_keldysh_contour_gf




  !##################################################################
  ! KELDYSH-BAYM-MATSUBARS CONTOUR GREEN'S FUNCTIONS:
  !##################################################################
  type :: kbm_contour_gf
     complex(8),dimension(:,:),pointer  :: less,ret
     complex(8),dimension(:,:),pointer  :: lmix
     real(8),dimension(:,:),pointer     :: mats
     logical                            :: status=.false.
     integer                            :: N=0,L=0
  end type kbm_contour_gf
  public :: kbm_contour_gf
  public :: allocate_kbm_contour_gf
  public :: deallocate_kbm_contour_gf
  public :: write_kbm_contour_gf
  public :: inquire_kbm_contour_gf
  public :: read_kbm_contour_gf
  public :: plot_kbm_contour_gf




  interface operator(*)
     module procedure &
          keldysh_contour_gf_scalarL_d,keldysh_contour_gf_scalarL_c,&
          keldysh_contour_gf_scalarR_d,keldysh_contour_gf_scalarR_c,&
          kbm_contour_gf_scalarL_d,kbm_contour_gf_scalarL_c,&
          kbm_contour_gf_scalarR_d,kbm_contour_gf_scalarR_c
  end interface operator(*)
  interface assignment(=)
     module procedure &
          keldysh_contour_gf_equality,&
          keldysh_contour_gf_equality_,&
          kbm_contour_gf_equality_
  end interface assignment(=)
  public :: assignment(=)
  public :: operator(*)


contains



  !################################################################################
  !########### KELDYSH CONTOUR GREEN'S FUNCTION (REAL-TIME ONLY) ##################
  !################################################################################
  subroutine keldysh_contour_gf_equality(G1,G2)
    type(keldysh_contour_gf),intent(inout) :: G1
    type(keldysh_contour_gf),intent(in)    :: G2
    G1%less = G2%less
    G1%ret = G2%ret
  end subroutine keldysh_contour_gf_equality

  subroutine keldysh_contour_gf_equality_(G1,C)
    type(keldysh_contour_gf),intent(inout) :: G1
    complex(8),intent(in) :: C
    G1%less = C
    G1%ret  = C
  end subroutine keldysh_contour_gf_equality_

  subroutine allocate_keldysh_contour_gf(G,N)
    type(keldysh_contour_gf) :: G
    integer                  :: i,j,N
    nullify(G%less,G%ret)
    G%N=N
    allocate(G%less(N,N),G%ret(N,N))
    G%less=zero
    G%ret =zero
    G%status=.true.
  end subroutine allocate_keldysh_contour_gf

  subroutine deallocate_keldysh_contour_gf(G)
    type(keldysh_contour_gf) :: G
    deallocate(G%less,G%ret)
    G%N=0
    G%status=.false.
  end subroutine deallocate_keldysh_contour_gf

  subroutine write_keldysh_contour_gf(G,file)
    type(keldysh_contour_gf) :: G
    character(len=*)         :: file
    integer                  :: N
    N=G%N
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         call error("ERROR contour_gf/write_keldysh_contour_gf: wrong dimensions")
    call store_data(trim(file)//"_less.data",G%less(:,:))
    call store_data(trim(file)//"_ret.data",G%ret(:,:))
  end subroutine write_keldysh_contour_gf

  subroutine read_keldysh_contour_gf(G,file)
    type(keldysh_contour_gf) :: G
    character(len=*)     :: file
    integer              :: N,L
    N=G%N
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         call error("ERROR contour_gf/write_keldysh_contour_gf: wrong dimensions")
    call read_data(trim(file)//"_less.data",G%less(:,:))
    call read_data(trim(file)//"_ret.data",G%ret(:,:))
  end subroutine read_keldysh_contour_gf

  function inquire_keldysh_contour_gf(file) result(check)
    logical          :: check,bool1,bool2
    character(len=*) :: file
    inquire(file=trim(file)//"_less.data",exist=bool1)
    if(.not.bool1)inquire(file=trim(file)//"_less.data.gz",exist=bool1)
    !if(.not.bool1)call warning("Can not read "//trim(file)//"_less.data")
    inquire(file=trim(file)//"_ret.data",exist=bool2)
    if(.not.bool2)inquire(file=trim(file)//"_ret.data.gz",exist=bool2)
    !if(.not.bool2)call warning("Can not read "//trim(file)//"_ret.data")
    check=bool1.AND.bool2
  end function inquire_keldysh_contour_gf

  subroutine plot_keldysh_contour_gf(G,t,file)
    type(keldysh_contour_gf)  :: G
    character(len=*)      :: file
    real(8),dimension(:) :: t
    integer               :: N
    N=G%N
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) .OR. size(t)/=N )&
         call error("ERROR contour_gf/plot_keldysh_contour_gf: wrong dimensions")
    call splot3d(trim(file)//"_less_t_t",t(:),t(:),G%less(:,:))
    call splot3d(trim(file)//"_ret_t_t",t(:),t(:),G%ret(:,:))
  end subroutine plot_keldysh_contour_gf

  function keldysh_contour_gf_scalarL_d(C,G) result(F)
    real(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(:,:)= C*G%less(:,:)
    F%ret(:,:) = C*G%ret(:,:)
  end function keldysh_contour_gf_scalarL_d

  function keldysh_contour_gf_scalarL_c(C,G) result(F)
    complex(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(:,:)=C*G%less(:,:)
    F%ret(:,:)=C*G%ret(:,:)
  end function keldysh_contour_gf_scalarL_c


  function keldysh_contour_gf_scalarR_d(G,C) result(F)
    real(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(:,:)=G%less(:,:)*C
    F%ret(:,:)=G%ret(:,:)*C
  end function keldysh_contour_gf_scalarR_d


  function keldysh_contour_gf_scalarR_c(G,C) result(F)
    complex(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less(:,:)=G%less(:,:)*C
    F%ret(:,:)=G%ret(:,:)*C
  end function keldysh_contour_gf_scalarR_c


  !******************************************************************
  !******************************************************************
  !******************************************************************






  !################################################################################
  !############ KADANOFF-BAYM-MATSUBARA CONTOUR GREEN'S FUNCTION ##################
  !################################################################################
  subroutine allocate_kbm_contour_gf(G,N,L)
    type(kbm_contour_gf)  :: G
    integer                 :: i,j,N,L
    nullify(G%less,G%ret,G%lmix,G%mats)
    G%N=N
    G%L=L
    allocate(G%less(N,N)) ; G%less=zero
    allocate(G%ret(N,N))  ; G%ret=zero
    allocate(G%lmix(N,L)) ; G%lmix=zero
    allocate(G%mats(L,L)) ; G%mats=zero
    G%status=.true.
  end subroutine allocate_kbm_contour_gf

  subroutine deallocate_kbm_contour_gf(G)
    type(kbm_contour_gf) :: G
    deallocate(G%less,G%ret,G%lmix,G%mats)
    G%N=0
    G%L=0
    G%status=.false.
  end subroutine deallocate_kbm_contour_gf

  subroutine write_kbm_contour_gf(G,file)
    type(kbm_contour_gf) :: G
    character(len=*)     :: file
    integer              :: N,L
    N=G%N ; L=G%L
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 1")
    if( size(G%lmix)/=N*L)&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 2")
    if( size(G%mats)/=L*L)&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 3")
    call store_data(trim(file)//"_less.data",G%less(:,:))
    call store_data(trim(file)//"_ret.data", G%ret(:,:))
    call store_data(trim(file)//"_lmix.data",G%lmix(:,:))
    call store_data(trim(file)//"_mats.data",G%mats(:,:))
  end subroutine write_kbm_contour_gf

  function inquire_kbm_contour_gf(file) result(check)
    integer          :: i
    logical          :: check,bool(5)
    character(len=*) :: file
    character(len=16),dimension(4)  :: ctype=(['less','ret','lmix','mats'])
    check=.true.
    do i=1,4
       inquire(file=trim(file)//"_"//trim(ctype(i))//".data",exist=bool(i))
       if(.not.bool(i))inquire(file=trim(file)//"_"//trim(ctype(i))//".data.gz",exist=bool(i))
       !if(.not.bool(i))call warning("Can not read "//trim(file)//"_"//trim(ctype(i))//".data")
       check=check.AND.bool(i)
    enddo
  end function inquire_kbm_contour_gf

  subroutine read_kbm_contour_gf(G,file)
    type(kbm_contour_gf) :: G
    character(len=*)     :: file
    integer              :: N,L
    N=G%N ; L=G%L
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 1")
    if( size(G%lmix)/=N*L)&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 2")
    if( size(G%mats)/=L*L)&
         call error("ERROR contour_gf/write_kbm_contour_gf: wrong dimensions 3")
    call read_data(trim(file)//"_less.data",G%less(:,:))
    call read_data(trim(file)//"_ret.data",G%ret(:,:))
    call read_data(trim(file)//"_lmix.data",G%lmix(:,:))
    call read_data(trim(file)//"_mats.data",G%mats(:,:))
  end subroutine read_kbm_contour_gf

  subroutine plot_kbm_contour_gf(G,t,tau,file)
    type(kbm_contour_gf)  :: G
    character(len=*)      :: file
    real(8),dimension(:)  :: t,tau
    integer               :: N,L
    N=G%N ; L=G%L
    if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 1")
    if( size(G%lmix)/=N*L)&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 2")
    if( size(G%mats)/=L*L)&
         call error("ERROR contour_gf/plot_kbm_contour_gf: wrong dimensions 3")
    call splot3d(trim(file)//"_less_t_t",t(:),t(:),G%less(:,:))
    call splot3d(trim(file)//"_ret_t_t",t(:),t(:),G%ret(:,:))
    call splot3d(trim(file)//"_lmix_t_tau",t(:),tau(:),G%lmix(:,:))
    call splot3d(trim(file)//"_mats_tau_tau",tau(:),tau(:),G%mats(:,:))
  end subroutine plot_kbm_contour_gf

  subroutine kbm_contour_gf_equality_(G1,C)
    type(kbm_contour_gf),intent(inout) :: G1
    complex(8),intent(in) :: C
    G1%less(:,:) = C
    G1%ret(:,:) = C
    G1%lmix(:,:) = C
    G1%mats(:,:) = C
  end subroutine kbm_contour_gf_equality_

  function kbm_contour_gf_scalarL_d(C,G) result(F)
    real(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(:,:)= C*G%less(:,:)
    F%ret(:,:) = C*G%ret(:,:)
    F%lmix(:,:)= C*G%lmix(:,:)
    F%mats(:,:)= C*G%mats(:,:)
  end function kbm_contour_gf_scalarL_d

  function kbm_contour_gf_scalarL_c(C,G) result(F)
    complex(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(:,:)=C*G%less(:,:)
    F%ret(:,:)=C*G%ret(:,:)
    F%lmix(:,:)=C*G%lmix(:,:)
    F%mats(:,:)=C*G%mats(:,:)
  end function kbm_contour_gf_scalarL_c

  function kbm_contour_gf_scalarR_d(G,C) result(F)
    real(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(:,:)=G%less(:,:)*C
    F%ret(:,:)=G%ret(:,:)*C
    F%lmix(:,:)=G%lmix(:,:)*C
    F%mats(:,:)=G%mats(:,:)*C
  end function kbm_contour_gf_scalarR_d

  function kbm_contour_gf_scalarR_c(G,C) result(F)
    complex(8),intent(in) :: C
    type(kbm_contour_gf),intent(in) :: G
    type(kbm_contour_gf) :: F
    F%less(:,:)=G%less(:,:)*C
    F%ret(:,:)=G%ret(:,:)*C
    F%lmix(:,:)=G%lmix(:,:)*C
    F%mats(:,:)=G%mats(:,:)*C  
  end function kbm_contour_gf_scalarR_c


END MODULE CONTOUR_GF
