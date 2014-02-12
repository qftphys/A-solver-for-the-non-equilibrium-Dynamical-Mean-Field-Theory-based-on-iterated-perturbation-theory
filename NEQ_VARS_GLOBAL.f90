!#####################################################################
!     PURPOSE  : Defines the global variables used thru all the code
!     AUTHORS  : Adriano Amaricci
!#####################################################################
MODULE NEQ_VARS_GLOBAL
  !Local:
  USE CONTOUR_GF
  !SciFor library
  USE SCIFOR_VERSION
  USE COMMON_VARS
  USE PARSE_CMD
  USE ARRAYS
  USE VECTORS
  USE IOTOOLS
  USE FUNCTIONS
  USE DERIVATE
  USE TOOLS
  implicit none

  !Version revision
  include "revision.inc"

  !Gloabl  variables
  !=========================================================
  integer                                :: Ntime         !Number of Time steps
  integer                                :: Ntau          !Imaginary time slices
  real(8)                                :: dt,dtau       !time step
  real(8)                                :: fmesh         !freq. step
  real(8)                                :: beta          !inverse temperature
  real(8)                                :: wmax          !min/max frequency
  real(8)                                :: tmax          !min/max time
  integer                                :: Lfreq         !Number of frequencies
  integer                                :: Lk            !total lattice  dimension
  integer                                :: Lkreduced     !reduced lattice dimension
  integer                                :: Nx            !lattice grid dimensions
  integer                                :: nloop         !dmft loop variables
  real(8)                                :: ts            !n.n. hopping amplitude
  real(8)                                :: Ui            !equilibrium local interaction
  real(8)                                :: Usti          !equilibrium local interaction
  real(8)                                :: U             !non-equilibrium local interaction
  real(8)                                :: Ust           !non-equilibrium local interaction
  real(8)                                :: Vbath         !coupling to the Thermostat
  real(8)                                :: Wbath         !Width of the BATH DOS
  real(8)                                :: Walpha        !exponent of the pseudo-gapped bath.
  real(8)                                :: Wgap          !gap of the gapped bath
  character(len=16)                      :: bath_type     !choose the shape of the BATH
  real(8)                                :: eps           !broadening
  real(8)                                :: dmft_error     !convergence error threshold
  character(len=16)                      :: field_type    !choose the profile of the electric field
  integer                                :: Nsuccess      !number of convergence success
  real(8)                                :: weight        !mixing weight parameter
  logical                                :: fchi
  logical                                :: plot3D,ifourth
  integer                                :: fupdate       !flag to decide WFupdate procedure


  !CONTAINER FOR THE CONTOUR PARAMETERS 
  !=========================================================  
  type(kb_contour_params)                :: cc_params


  !LATTICE (weight & dispersion) ARRAYS:
  !=========================================================  
  real(8),dimension(:),allocatable       :: wt,epsik


  !ELECTRIC FIELD VARIABLES (& NML):
  !=========================================================  
  type(vect2D)                           :: Ak,Ek         !Electric field vector potential and vector
  real(8)                                :: Efield        !Electric field strength
  real(8)                                :: Ex,Ey         !Electric field vectors as input
  real(8)                                :: t0,t1         !turn on/off time, t0 also center of the pulse
  integer                                :: Ncycles       !Number of cycles in pulsed light packet
  real(8)                                :: omega0        !parameter for the Oscilatting field and Pulsed light
  real(8)                                :: E1            !Electric field strenght for the AC+DC case (tune to resonate)


  !SUSCEPTIBILITY
  !=========================================================  
  real(8),allocatable,dimension(:,:,:,:) :: chi


  !FILES & DATA DIRECTORY:
  !=========================================================
  character(len=32)                      :: g0file
  character(len=32)                      :: plot_dir



  !NAMELISTS:
  !=========================================================
  namelist/variables/&
       dt,Ntime,Ntau,beta,Ui,U,Usti,Ust,ts,Lfreq,eps,  &    ! INPUT
       bath_type,Vbath,Wbath,Walpha,Wgap,  &    ! BATH
       Nx,Lkreduced,                       &    ! LATTICE
       Efield,Ex,Ey,field_type,t0,t1,      &    !
       Ncycles,omega0,E1,                  &    ! ELECTRIC FIELD
       nloop,dmft_error,nsuccess,weight,    &    ! DMFT
       plot3D,fchi,ifourth,fupdate,     &    ! FLAGS
       g0file,plot_dir           ! FILES & DIRS



contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input_init(inputFILE)
    character(len=*)               :: inputFILE
    character(len=256),allocatable :: help_buffer(:)
    integer                        :: i
    logical                        :: control

    call version(revision)

    !INPUT
    dt         = 0.1d0
    Ntime      = 50
    Ntau       = 25
    beta       = 10.d0
    Ui         = 2.d0
    U          = 2.d0
    ts         = 1.d0
    Lfreq      = 2048  
    eps        = 0.01d0
    !BATH
    bath_type  = 'flat'
    Vbath      = 0.d0
    Wbath      = 20.d0
    Walpha     = 1.d0
    Wgap       = 5.d0
    !
    !LATTICE
    Nx           = 25
    Lkreduced    = 300
    !FIELD:
    Efield       = 0.d0
    field_type   = 'dc'
    Ex           = 1.d0
    Ey           = 0.d0
    t0           = 0.d0
    t1           = 1.d9
    Ncycles      = 1
    omega0       = 1.d0*pi
    E1           = 0.d0
    !DMFT
    nloop        = 30
    dmft_error    = 1.d-3
    nsuccess     = 2
    weight       = 1.d0
    !FLAGS:
    plot3D       = .false.
    fchi         = .false.
    ifourth      = .false.
    fupdate      = 0
    !FILES&DIR:
    g0file   = 'g0.restart'
    plot_dir     = 'PLOT'

    inquire(file=adjustl(trim(inputFILE)),exist=control)
    if(control)then
       open(10,file=adjustl(trim(inputFILE)))
       read(10,nml=variables)
       close(10)
    else
       print*,"Can not find INPUT file"
       print*,"Dumping a default version in default."//trim(inputFILE)
       call dump_input_file("default.")
       stop 
    endif

    !GLOBAL
    call parse_cmd_variable(dt         ,"DT")
    call parse_cmd_variable(Ntau       ,"NTAU")
    call parse_cmd_variable(Nx         ,"NX")
    call parse_cmd_variable(Ntime      ,"NTIME")
    call parse_cmd_variable(beta       ,"BETA")
    call parse_cmd_variable(Ui         ,"UI")
    call parse_cmd_variable(Usti       ,"USTI")
    call parse_cmd_variable(U          ,"U")
    call parse_cmd_variable(Ust        ,"UST")
    call parse_cmd_variable(ts         ,"TS")
    call parse_cmd_variable(Lfreq      ,"LFREQ")
    call parse_cmd_variable(eps        ,"EPS")
    !BATH
    call parse_cmd_variable(Vbath      ,"VBATH")
    call parse_cmd_variable(bath_type  ,"BATH_TYPE")
    call parse_cmd_variable(wbath      ,"WBATH")
    call parse_cmd_variable(walpha     ,"WALPHA")
    call parse_cmd_variable(wgap       ,"WGAP")
    !LATTICE
    call parse_cmd_variable(Nx         ,"NX")
    call parse_cmd_variable(Lkreduced  ,"LKREDUCED")
    !EFIELD
    call parse_cmd_variable(field_type ,"FIELD_TYPE")
    call parse_cmd_variable(Efield     ,"EFIELD")
    call parse_cmd_variable(Ex         ,"EX")
    call parse_cmd_variable(Ey         ,"EY")
    call parse_cmd_variable(t0         ,"T0")
    call parse_cmd_variable(t1         ,"T1")
    call parse_cmd_variable(ncycles    ,"NCYCLES")
    call parse_cmd_variable(omega0     ,"OMEGA0")
    call parse_cmd_variable(E1         ,"E1")
    !DMFT
    call parse_cmd_variable(nloop      ,"NLOOP")
    call parse_cmd_variable(dmft_error  ,"DMFT_ERROR")
    call parse_cmd_variable(Nsuccess   ,"NSUCCESS")
    call parse_cmd_variable(weight     ,"WEIGHT")
    !FLAGS:
    call parse_cmd_variable(plot3D     ,"PLOT3D")
    call parse_cmd_variable(fchi       ,"FCHI")
    call parse_cmd_variable(ifourth    ,"IFOURTH")
    call parse_cmd_variable(fupdate    ,"FUPDATE")
    !FILES&DIR:
    call parse_cmd_variable(g0file ,"G0FILE")
    call parse_cmd_variable(plot_dir   ,"PLOT_DIR")

    !if(U==0.d0)Nloop=1

    if(mpiID==0)then
       write(*,*)"CONTROL PARAMETERS"
       write(*,nml=variables)
       write(*,*)"--------------------------------------------"
       write(*,*)""
       call dump_input_file("used.")
    endif

    if(plot3D)call create_data_dir(trim(plot_dir))

  contains
    subroutine dump_input_file(prefix)
      character(len=*) :: prefix
      if(mpiID==0)then
         open(10,file=trim(adjustl(trim(prefix)))//adjustl(trim(inputFILE)))
         write(10,nml=variables)
         close(10)
      endif
    end subroutine dump_input_file
  end subroutine read_input_init











end module NEQ_VARS_GLOBAL

