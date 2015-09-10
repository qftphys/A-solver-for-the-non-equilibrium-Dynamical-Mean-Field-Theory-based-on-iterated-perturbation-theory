MODULE NEQ_CONTOUR_GF
  USE NEQ_CONTOUR
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS
  implicit none
  private


  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS:
  type,public :: kb_contour_gf
     complex(8),dimension(:,:),allocatable :: less
     complex(8),dimension(:,:),allocatable :: ret
     complex(8),dimension(:,:),allocatable :: lmix
     real(8),dimension(:),allocatable      :: mats
     complex(8),dimension(:),allocatable   :: iw
     logical                               :: status=.false.
  end type kb_contour_gf
  !
  !
  !
  ! KADANOFF-BAYM CONTOUR SIGMA FUNCTION:
  type,public :: kb_contour_sigma
     type(kb_contour_gf)                   :: reg
     complex(8),dimension(:),allocatable   :: hf
     logical                               :: status=.false.
  end type kb_contour_sigma
  !
  !
  !  
  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS DERIVATIVE
  type,public :: kb_contour_dgf
     complex(8),dimension(:),allocatable   :: less,gtr
     complex(8),dimension(:),allocatable   :: ret
     complex(8),dimension(:),allocatable   :: lmix
     logical                               :: status=.false.
     logical                               :: anomalous=.false.
  end type kb_contour_dgf
  !
  !
  !
  ! ALLOCATION/DEALLOCATION ROUTINES:
  interface allocate_kb_contour_gf
     module procedure allocate_kb_contour_gf_main
  end interface allocate_kb_contour_gf
  interface deallocate_kb_contour_gf
     module procedure deallocate_kb_contour_gf_main
  end interface deallocate_kb_contour_gf
  public :: allocate_kb_contour_gf
  public :: deallocate_kb_contour_gf
  !
  interface allocate_kb_contour_sigma
     module procedure allocate_kb_contour_sigma_main
  end interface allocate_kb_contour_sigma
  interface deallocate_kb_contour_sigma
     module procedure deallocate_kb_contour_sigma_main
  end interface deallocate_kb_contour_sigma
  public :: allocate_kb_contour_sigma
  public :: deallocate_kb_contour_sigma
  !
  interface allocate_kb_contour_dgf
     module procedure allocate_kb_contour_dgf_main
  end interface allocate_kb_contour_dgf
  interface deallocate_kb_contour_dgf
     module procedure deallocate_kb_contour_dgf_main
  end interface deallocate_kb_contour_dgf
  public :: allocate_kb_contour_dgf
  public :: deallocate_kb_contour_dgf
  !
  !
  !VIE/VIDE SOLVER:
  interface vie_kb_contour_gf
     module procedure vie_kb_contour_gf_Q
     module procedure vie_kb_contour_gf_Sigma
     module procedure vie_kb_contour_gf_Delta
  end interface vie_kb_contour_gf
  public :: vie_kb_contour_gf
  !
  interface vide_kb_contour_gf
     module procedure vide_kb_contour_gf_K
     module procedure vide_kb_contour_gf_Sigma
  end interface vide_kb_contour_gf
  public :: vide_kb_contour_gf
  !
  !
  ! ADD (TOTAL DOMAIN) ROUTINES:
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
  ! SUM (PERIMETER) ROUTINES:
  interface sum_kb_contour_gf
     module procedure sum_kb_contour_gf_simple
     module procedure sum_kb_contour_gf_delta_d
     module procedure sum_kb_contour_gf_delta_c
     module procedure sum_kb_contour_gf_recursive
  end interface sum_kb_contour_gf
  public :: sum_kb_contour_gf
  !
  !
  ! DELETE (RESET PERIMETER) ROUTINES:
  public :: del_kb_contour_gf
  !
  !
  ! CONVOLUTION:
  interface convolute_kb_contour_gf
     module procedure convolute_kb_contour_gf_gf
     module procedure convolute_kb_contour_sigma_gf
     module procedure convolute_kb_contour_gf_sigma
     module procedure convolute_kb_contour_gf_gf_recursive
     module procedure convolute_kb_contour_gf_sigma_recursive
     ! module procedure convolute_kb_contour_gf_delta_recursive
     module procedure convolute_kb_contour_gf_delta_left
     module procedure convolute_kb_contour_gf_delta_right
  end interface convolute_kb_contour_gf
  public :: convolute_kb_contour_gf
  !
  !
  ! OTHER ROUTINES && PLOT:  
  public :: extrapolate_kb_contour_gf
  public :: save_kb_contour_gf
  public :: inquire_kb_contour_gf
  public :: read_kb_contour_gf
  public :: plot_kb_contour_gf
  public :: extrapolate_kb_contour_sigma
  public :: save_kb_contour_sigma
  public :: inquire_kb_contour_sigma
  public :: read_kb_contour_sigma
  public :: plot_kb_contour_sigma
  interface check_dimension_kb_contour
     module procedure check_dimension_kb_contour_gf
     module procedure check_dimension_kb_contour_gf_
     module procedure check_dimension_kb_contour_sigma
     module procedure check_dimension_kb_contour_sigma_
     module procedure check_dimension_kb_contour_dgf
     module procedure check_dimension_kb_contour_dgf_
  end interface check_dimension_kb_contour
  public :: check_dimension_kb_contour
  !
  ! 
  !SYMMETRIES & COMPONENTS:
  public :: get_gtr
  public :: get_adv
  public :: get_rmix
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
  interface assignment(=)
     module procedure kb_contour_gf_equality_
     module procedure kb_contour_dgf_equality_
     module procedure kb_contour_gf_equality__
     module procedure kb_contour_dgf_equality__
  end interface assignment(=)
  public :: operator(*)
  public :: assignment(=)
  !
  !
  !INTEGRATION ROUTINES:
  interface kb_trapz
     module procedure kb_trapz_d
     module procedure kb_trapz_c
  end interface kb_trapz
  public :: kb_trapz
  !
  interface kb_half_trapz
     module procedure kb_half_trapz_d
     module procedure kb_half_trapz_c
  end interface kb_half_trapz
  public :: kb_half_trapz



contains


  !======= ALLOCATE/DEALLOCATE ======= 
  include "neq_contour_gf_allocate.f90"



  !======= ADD =======
  !C(t,t')=A(t,t') + B(t,t'), with t,t'=0,t_max
  include "neq_contour_gf_add.f90"



  !======= SUM =======
  ! performs the sum C=a*A + b*B along the perimeter t=N*dt & t'=0,N*dt
  ! when called multiple times it sums up to a given array. 
  include "neq_contour_gf_sum.f90"



  !======= CONVOLUTION ======= 
  include "neq_contour_gf_convolute.f90"



  !======= VOLTERRA INTEGRAL EQUATION =======
  include "neq_contour_gf_vie.f90"



  !======= VOLTERRA INTEGRO-DIFFERENTIAL EQUATION =======
  include "neq_contour_gf_vide.f90"



  !======= MISCELLANEOUS ======= 
  ! CHECK DIMENSION
  ! SAVE
  ! READ
  ! PLOT
  ! EXTRAPOLATION
  include "neq_contour_gf_misc.f90"



  !======= OPERATIONS =======
  ! .*. times a scalar
  ! .=. equality among GF
  ! +get_Adv
  ! +get_Gtr
  ! +get_Rmix
  ! +get_Bar
  include "neq_contour_gf_operations.f90"



  !======= INTEGRATION =======
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
  !
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
  !
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






END MODULE NEQ_CONTOUR_GF

