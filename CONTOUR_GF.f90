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
     !one-time index wrapping the 2D array (lower-right triangle)
     complex(8),dimension(:),allocatable  :: less,gtr 
     logical                              :: status=.false.
     integer                              :: N=0
  end type keldysh_contour_gf

  public :: keldysh_contour_gf
  public :: allocate_keldysh_contour_gf
  public :: deallocate_keldysh_contour_gf
  public :: write_keldysh_contour_gf
  public :: read_keldysh_contour_gf
  public :: plot_keldysh_contour_gf
  public :: plot_keldysh_contour_tri_gf
  public :: inquire_keldysh_contour_gf

  public :: pack_index
  public :: pack_index_tri
  public :: unpack_index_i,unpack_index_j

  public :: pack_less_tri , pack_less
  public :: pack_gtr_tri , pack_gtr

  interface operator(*)
     module procedure &
          keldysh_contour_gf_scalarL_d,&
          keldysh_contour_gf_scalarL_c,&
          keldysh_contour_gf_scalarR_d,&
          keldysh_contour_gf_scalarR_c
  end interface operator(*)

  interface assignment(=)
     module procedure &
          keldysh_contour_gf_equality,&
          keldysh_contour_gf_equality_
  end interface assignment(=)

  public :: assignment(=)
  public :: operator(*)


contains

  !##################################################################
  ! WRAPPING INDEX MAPS:
  !##################################################################
  !map 2d(i,j) coordinates to 1Dlong array k
  function pack_index(i,j,n) result(k)
    integer,intent(in) :: i,j,n
    integer            :: k
    k = (i-1)*n + j
  end function pack_index
  !--------------------------------------------

  !--------------------------------------------
  !map the 2d(i,j) of lower-right triangle of t,t' matrices to single 1D index k
  function pack_index_tri(i,j) result(k)
    integer,intent(in) :: i,j
    integer            :: k
    k = j + i*(i-1)/2
  end function  pack_index_tri
  !--------------------------------------------

  !--------------------------------------------
  !unpack the wrapping index: first coordinate (t)
  function unpack_index_i(k,N) result(i)
    integer :: k,i,N
    i=int((k-1)/N+1.d-5)+1
  end function unpack_index_i
  !--------------------------------------------

  !--------------------------------------------
  !unpack the wrapping index: second coordinate (t')
  function unpack_index_j(k,N) result(j)
    integer k,N,j
    j=mod(k-1,N)+1
  end function unpack_index_j
  !--------------------------------------------

  !--------------------------------------------
  function pack_less_tri(i,j,G) result(gc)
    integer                    :: i,j,k,Nt
    type(keldysh_contour_gf)   :: G
    complex(8)                 :: gc
    if(i>=j)then
       k = pack_index_tri(i,j)
       gc = G%less(k)
    else
       k = pack_index_tri(j,i)
       gc=-conjg(G%less(k))
    endif
  end function pack_less_tri
  !--------------------------------------------

  !--------------------------------------------
  function pack_less(i,j,Nt,G) result(gc)
    integer                    :: i,j,k,Nt
    type(keldysh_contour_gf)   :: G
    complex(8)                 :: gc
    if(i>=j)then
       k = pack_index(i,j,Nt)
       gc = G%less(k)
    else
       k = pack_index(j,i,Nt)
       gc=-conjg(G%less(k))
    endif
  end function pack_less
  !--------------------------------------------

  !--------------------------------------------
  function pack_gtr_tri(i,j,G) result(gc)
    integer                  :: i,j,k
    type(keldysh_contour_gf) :: G
    complex(8)               :: gc
    if(i>=j)then
       k = pack_index_tri(i,j)
       gc = G%gtr(k)
    else
       k = pack_index_tri(j,i)
       gc=-conjg(G%gtr(k))
    endif
  end function pack_gtr_tri
  !--------------------------------------------

  !--------------------------------------------
  function pack_gtr(i,j,Nt,G) result(gc)
    integer                    :: i,j,k,Nt
    type(keldysh_contour_gf)   :: G
    complex(8)                 :: gc
    if(i>=j)then
       k = pack_index(i,j,Nt)
       gc = G%gtr(k)
    else
       k = pack_index(j,i,Nt)
       gc=-conjg(G%gtr(k))
    endif
  end function pack_gtr

  !******************************************************************
  !******************************************************************
  !******************************************************************

  !##################################################################
  !####### KELDYSH CONTOUR GREEN'S FUNCTION (REAL-TIME ONLY) ########
  !##################################################################

  subroutine allocate_keldysh_contour_gf(G,N)
    type(keldysh_contour_gf) :: G
    integer                  :: i,j,N
    if(allocated(G%less))deallocate(G%less)
    if(allocated(G%gtr))deallocate(G%gtr)
    G%N=N
    allocate(G%less(N),G%gtr(N))
    G%less=zero
    G%gtr =zero
    G%status=.true.
  end subroutine allocate_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine deallocate_keldysh_contour_gf(G)
    type(keldysh_contour_gf) :: G
    deallocate(G%less,G%gtr)
    G%N=0
    G%status=.false.
  end subroutine deallocate_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine write_keldysh_contour_gf(G,file)
    type(keldysh_contour_gf) :: G
    character(len=*)         :: file
    integer                  :: N
    N=G%N
    if( (size(G%less)/=N) .OR. (size(G%gtr)/=N) )&
         call error("ERROR contour_gf/write_keldysh_contour_gf: wrong dimensions")
    call store_data(trim(file)//"_less.data",G%less)
    call store_data(trim(file)//"_gtr.data",G%gtr)
  end subroutine write_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine read_keldysh_contour_gf(G,file)
    type(keldysh_contour_gf) :: G
    character(len=*)         :: file
    integer                  :: N,L
    logical                  :: check
    N=G%N
    check = inquire_keldysh_contour_gf(file)
    if(.not.check)call error("READ_KELDYSH_CONTOUR_GF: can not read data from "//reg(file))
    call read_data(trim(file)//"_less.data",G%less)
    call read_data(trim(file)//"_gtr.data",G%gtr)
  end subroutine read_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function inquire_keldysh_contour_gf(file) result(check)
    logical          :: check,bool1,bool2
    character(len=*) :: file
    inquire(file=trim(file)//"_less.data",exist=bool1)
    if(.not.bool1)inquire(file=trim(file)//"_less.neq.gz",exist=bool1)
    !if(.not.bool1)call msg("Can not read "//reg(file)//"_less.neq")
    inquire(file=trim(file)//"_gtr.data",exist=bool2)
    if(.not.bool2)inquire(file=trim(file)//"_gtr.neq.gz",exist=bool2)
    !if(.not.bool2)call msg("Can not read "//reg(file)//"_gtr.neq")
    check=bool1.AND.bool2
  end function inquire_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine plot_keldysh_contour_gf(t,G,file)
    type(keldysh_contour_gf)              :: G
    complex(8),dimension(:,:),allocatable :: Gless,Ggtr
    character(len=*)                      :: file
    real(8),dimension(:)                  :: t
    integer                               :: N,i,j,k
    N=int(sqrt(dble(G%N)))
    if( size(t)/=N )&
         call error("ERROR plot_keldysh_contour_gf: wrong dimensions"//reg(file))
    allocate(Gless(N,N),Ggtr(N,N))
    do i=1,N
       do j=1,N
          Gless(i,j) = pack_less(i,j,N,G)
          Ggtr(i,j)  = pack_gtr(i,j,N,G)
       enddo
    enddo
    call splot3d(trim(file)//"_less_t_t.neq",t,t,Gless)
    call splot3d(trim(file)//"_gtr_t_t.neq",t,t,Ggtr)
    deallocate(Gless,Ggtr)
  end subroutine plot_keldysh_contour_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine plot_keldysh_contour_tri_gf(t,G,N,file)    
    type(keldysh_contour_gf)  :: G
    real(8),dimension(N)      :: t
    integer                   :: N
    complex(8),dimension(N,N) :: Gless,Ggtr
    character(len=*)          :: file
    integer                   :: i,j,k

    if( G%N /= N*(N+1)/2 )&
         call error("ERROR plot_keldysh_contour_tri_gf: wrong dimensions"//reg(file))
    do i=1,N
       do j=1,N
          Gless(i,j) = pack_less_tri(i,j,G)
          Ggtr(i,j)  = pack_gtr_tri(i,j,G)
       enddo
    enddo
    call splot3d(trim(file)//"_less_t_t.neq",t,t,Gless)
    call splot3d(trim(file)//"_gtr_t_t.neq",t,t,Ggtr)
  end subroutine plot_keldysh_contour_tri_gf

  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine keldysh_contour_gf_equality(G1,G2)
    type(keldysh_contour_gf),intent(inout) :: G1
    type(keldysh_contour_gf),intent(in)    :: G2
    G1%less = G2%less
    G1%gtr = G2%gtr
  end subroutine keldysh_contour_gf_equality
  !--------------------------------------------

  !--------------------------------------------
  subroutine keldysh_contour_gf_equality_(G1,C)
    type(keldysh_contour_gf),intent(inout) :: G1
    complex(8),intent(in) :: C
    G1%less = C
    G1%gtr = C
  end subroutine keldysh_contour_gf_equality_

  !******************************************************************
  !******************************************************************
  !******************************************************************


  function keldysh_contour_gf_scalarL_d(C,G) result(F)
    real(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less= C*G%less
    F%gtr = C*G%gtr
  end function keldysh_contour_gf_scalarL_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function keldysh_contour_gf_scalarL_c(C,G) result(F)
    complex(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less=C*G%less
    F%gtr=C*G%gtr
  end function keldysh_contour_gf_scalarL_c

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function keldysh_contour_gf_scalarR_d(G,C) result(F)
    real(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less=G%less*C
    F%gtr=G%gtr*C
  end function keldysh_contour_gf_scalarR_d

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function keldysh_contour_gf_scalarR_c(G,C) result(F)
    complex(8),intent(in) :: C
    type(keldysh_contour_gf),intent(in) :: G
    type(keldysh_contour_gf) :: F
    F%less=G%less*C
    F%gtr=G%gtr*C
  end function keldysh_contour_gf_scalarR_c


  !******************************************************************
  !******************************************************************
  !******************************************************************



END MODULE CONTOUR_GF
