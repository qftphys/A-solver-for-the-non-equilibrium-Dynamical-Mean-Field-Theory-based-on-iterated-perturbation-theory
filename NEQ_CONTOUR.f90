MODULE NEQ_CONTOUR
  USE CONSTANTS, only: one,xi,zero,pi
  USE ARRAYS, only: linspace,arange
  USE IOTOOLS, only: reg,free_unit
  implicit none
  private

  ! KADANOFF-BAYM-KELDYSH-SCHWINGER CONTOUR PARAMETERS
  !====================================================
  type,public                          :: kb_contour_params
     integer                           :: Ntime=0        !Largest dimension of the contour
     integer                           :: Ntau=0         !Largest dimension of the contour
     integer                           :: Niw=0          !Number of Matsubara Frequencies
     real(8)                           :: dt=0.d0        !real-time step
     real(8)                           :: beta=0.d0      !length of the im-time interval
     integer                           :: itime=0        !Actual time index
     real(8)                           :: time=0.d0      !actual time at it=Itime
     real(8)                           :: dtau=0.d0      !im-time step
     real(8)                           :: tmax=0.d0      !length of the real-time interval
     real(8),dimension(:),allocatable  :: t              !real-time array
     real(8),dimension(:),allocatable  :: tau            !im-time array
     real(8),dimension(:),allocatable  :: wm             !matsubara freq. array
     logical                           :: status=.false. !allocation status
  end type kb_contour_params

  public :: allocate_kb_contour_params
  public :: deallocate_kb_contour_params
  public :: setup_kb_contour_params
  public :: write_kb_contour_params
  public :: read_kb_contour_params
  public :: print_kb_contour_params
  !

  interface assignment(=)
     module procedure kb_contour_params_equality
  end interface assignment(=)
  public :: assignment(=)


contains



  !======= ALLOCATE ======= 
  subroutine allocate_kb_contour_params(params,Ntime,Ntau,Niw)
    type(kb_contour_params) :: params
    integer,intent(in)      :: Ntime,Ntau,Niw
    if(params%status)stop "neq_contour/allocate_kb_contour_params: Contour is already allocated"
    if(allocated(params%t))deallocate(params%t)
    if(allocated(params%tau))deallocate(params%tau)
    if(allocated(params%wm))deallocate(params%wm)
    params%Ntime= Ntime
    params%Ntau = Ntau
    params%Niw  = Niw
    allocate(params%t(Ntime))
    allocate(params%tau(Ntau))
    allocate(params%wm(Niw))
    params%status=.true.
  end subroutine allocate_kb_contour_params



  !======= DEALLOCATE ======= 
  subroutine deallocate_kb_contour_params(params)
    type(kb_contour_params) :: params
    if(.not.params%status)stop "neq_contour/deallocate_kb_contour_params: Contour not allocated"
    if(allocated(params%t))deallocate(params%t)
    if(allocated(params%tau))deallocate(params%tau)
    params%Ntime = 0
    params%Ntau  = 0
    params%Itime = 0
    params%dt    = 0.d0
    params%dtau  = 0.d0
    params%status=.false.
  end subroutine deallocate_kb_contour_params



  !======= BUILD ======= 
  subroutine setup_kb_contour_params(params,dt,beta)
    type(kb_contour_params) :: params
    real(8),intent(in)      :: dt,beta
    integer                 :: Ntime,Ntau,Niw
    real(8) :: tmax_,dtau_,dt_
    if(.not.params%status)stop "neq_contour/set_kb_contour_params: Contour not allocated"
    Ntime       = params%Ntime
    Ntau        = params%Ntau
    Niw         = params%Niw
    params%Itime= 1         !<== set the actual time_step to the minimum
    params%dt   = dt
    tmax_       = dt*(Ntime-1)
    params%tmax = tmax_
    params%beta = beta
    params%dtau = beta/(Ntau-1)
    params%t    = linspace(0d0,tmax_,Ntime,mesh=dt_)
    params%tau  = linspace(0d0,beta,Ntau,mesh=dtau_)
    params%wm   = pi/beta*(2*arange(1,Niw)-1)
    if(dt_/=params%dt)stop "neq_contour/set_kb_contour_params: dt != dt_ "
    if(dtau_/=params%dtau)stop "neq_contour/set_kb_contour_params: dtau != dtau_ "
    call print_kb_contour_params(params)
  end subroutine setup_kb_contour_params




  !======= WRITE ======= 
  subroutine write_kb_contour_params(params,file)
    type(kb_contour_params)  :: params
    character(len=*)     :: file
    integer              :: unit
    if(.not.params%status)stop "neq_contour/write_kb_contour_params: Contour not allocated"
    unit=free_unit()
    open(unit,file=reg(file))
    write(unit,*)params%Ntime
    write(unit,*)params%Ntau
    write(unit,*)params%Niw
    write(unit,*)params%dt
    write(unit,*)params%beta
    write(unit,*)params%itime
    write(unit,*)params%time
    write(unit,*)params%dtau
    write(unit,*)params%tmax
    write(unit,*)params%t
    write(unit,*)params%tau
    write(unit,*)params%wm
    close(unit)
  end subroutine write_kb_contour_params


  !======= PRINT ======= 
  subroutine print_kb_contour_params(params)
    type(kb_contour_params)  :: params
    integer :: i,Ntime,Ntau
    if(.not.params%status)stop "neq_contour/print_kb_contour_params: Contour not allocated"
    Ntime=params%Ntime
    Ntau =params%Ntau
    write(*,"(A10,A1,I15)")"Ntime","=",params%Ntime
    write(*,"(A10,A1,I15)")"Ntau","=",params%Ntau
    write(*,"(A10,A1,I15)")"Niw","=",params%Niw
    write(*,"(A10,A1,I15)")"Itime","=",params%itime
    write(*,"(A10,A1,F15.8)")"time","=",params%time
    write(*,"(A10,A1,F15.8)")"dt","=",params%dt
    write(*,"(A10,A1,F15.8)")"dtau","=",params%dtau
    write(*,"(A10,A1,F15.8)")"beta","=",params%beta
    write(*,"(A10,A1,F15.8)")"tmax","=",params%tmax
    write(*,"(A10,A1,3F15.8,A,3F15.8)")"t",  "=",(params%t(i),i=1,3),"   ... ",(params%t(i),i=Ntime-2,Ntime)
    write(*,"(A10,A1,3F15.8,A,3F15.8)")"tau","=",(params%tau(i),i=1,3),"   ... ",(params%tau(i),i=Ntau-2,Ntau)
  end subroutine print_kb_contour_params




  !======= READ ======= 
  subroutine read_kb_contour_params(params,file)
    type(kb_contour_params) :: params
    character(len=*)        :: file
    integer                 :: unit
    logical                 :: check
    integer                 :: Ntime,Ntau,Niw
    real(8)                 :: dt,beta
    if(params%status)stop "neq_contour/read_kb_contour_params: Contour already allocated"
    inquire(file=reg(file),exist=check)
    if(.not.check)stop "neq_contour/read_kb_contour_params: file does not exist"
    unit=free_unit()
    open(unit,file=reg(file))
    read(unit,*)Ntime
    read(unit,*)Ntau
    read(unit,*)Niw
    call allocate_kb_contour_params(params,Ntime,Ntau,Niw)
    read(unit,*)params%dt
    read(unit,*)params%beta
    read(unit,*)params%itime
    read(unit,*)params%time
    read(unit,*)params%dtau
    read(unit,*)params%tmax
    read(unit,*)params%t
    read(unit,*)params%tau
    read(unit,*)params%wm
    close(unit)
  end subroutine read_kb_contour_params




  !======= OPERATIONS = ======= 
  subroutine kb_contour_params_equality(P1,P2)
    type(kb_contour_params),intent(inout) :: P1
    type(kb_contour_params),intent(in)    :: P2
    if(.not.P2%status)stop "neq_contour/kb_contour_params_equality: P2 not allocated"
    if(P1%status)call deallocate_kb_contour_params(P1)
    call allocate_kb_contour_params(P1,P2%Ntime,P2%Ntau,P2%Niw)
    P1%Ntime  = P2%Ntime
    P1%Ntau   = P2%Ntau
    P1%Niw    = P2%Niw
    P1%dt     = P2%dt
    P1%beta   = P2%beta
    P1%Itime  = P2%Itime
    P1%time   = P2%time
    P1%dtau   = P2%dtau
    P1%tmax   = P2%tmax
    P1%t(:)   = P2%t(:)
    P1%tau(:) = P2%tau(:)
    P1%wm(:)  = P2%wm(:)
  end subroutine kb_contour_params_equality


END MODULE NEQ_CONTOUR


! #ifdef _test
! program test_NEQ_CONTOUR
!   USE NEQ_CONTOUR
!   implicit none
!   integer                 :: Ntime,Ntau,Niw
!   real(8)                 :: dt,beta
!   type(kb_contour_params) :: cc_params,cc2
!   Ntime = 10
!   Ntau  = 10
!   Niw    = 50
!   dt    = 0.1d0
!   beta  = 10d0
!   print*,"allocating"
!   call allocate_kb_contour_params(cc_params,Ntime,Ntau,Niw)
!   print*,"setup & print"
!   call setup_kb_contour_params(cc_params,dt,beta)
!   print*,"write"
!   call write_kb_contour_params(cc_params,"cc_params_file.test")
!   print*,"deallocate"
!   call deallocate_kb_contour_params(cc_params)
!   print*,"read"
!   call read_kb_contour_params(cc_params,"cc_params_file.test")
!   call print_kb_contour_params(cc_params)
!   cc2 = cc_params
!   call print_kb_contour_params(cc2)
!   call write_kb_contour_params(cc_params,"cc2_params_file.test")
! end program test_NEQ_CONTOUR
! #endif
