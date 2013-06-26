!#####################################################################
!     PURPOSE  : Defines the global variables used thru all the code
!     AUTHORS  : G.Mazza & A.Amaricci
!#####################################################################
!include "ReadMe"
MODULE NEQ_VARS_GLOBAL
  !Local:
  USE CONTOUR_GF
  !SciFor library
  USE SCIFOR_VERSION
  USE COMMON_VARS
  USE PARSE_CMD
  USE TIMER
  USE VECTORS
  USE SQUARE_LATTICE
  USE INTEGRATE
  USE IOTOOLS
  USE FFTGF
  USE FUNCTIONS
  USE INTERPOLATE
  USE TOOLS
  implicit none

  !Version revision
  include "revision.inc"

  !Gloabl  variables
  !=========================================================
  integer                                :: Nstep         !Number of Time steps
  integer                                :: Nfit
  integer                                :: L             !a big number
  integer                                :: Lk            !total lattice  dimension
  integer                                :: Lkreduced     !reduced lattice dimension
  integer                                :: Nx,Ny         !lattice grid dimensions
  integer                                :: iloop,nloop    !dmft loop variables
  integer                                :: eqnloop        !dmft loop of the equilibrium solution
  real(8)                                :: ts             !n.n./n.n.n. hopping amplitude
  real(8)                                :: u              !local,non-local interaction 
  real(8)                                :: Vbath          !Hopping amplitude to the BATH
  real(8)                                :: Wbath          !Width of the BATH DOS
  real(8)                                :: dt,dtau        !time step
  real(8)                                :: dtfit
  real(8)                                :: fmesh          !freq. step
  real(8)                                :: beta           !inverse temperature
  real(8)                                :: eps            !broadening
  character(len=16)                      :: int_method    !choose the integration method (rect,trapz,simps)
  character(len=16)                      :: bath_type     !choose the shape of the BATH
  character(len=16)                      :: field_type !choose the profile of the electric field
  real(8)                                :: eps_error     !convergence error threshold
  integer                                :: Nsuccess      !number of convergence success
  real(8)                                :: weight        !mixing weight parameter
  real(8)                                :: wmin,wmax     !min/max frequency
  real(8)                                :: tmin,tmax     !min/max time
  real(8)                                :: Walpha         !exponent of the pseudo-gapped bath.
  real(8)                                :: Wgap          !gap of the gapped bath
  logical                                :: plot3D,fchi
  logical                                :: solve_eq
  integer                                :: fupdate !flag to decide WFupdate procedure
  !

  !FILES TO RESTART
  !=========================================================
  character(len=32)                      :: irdSFILE,irdNkFILE


  !FREQS & TIME ARRAYS:
  !=========================================================  
  real(8),dimension(:),allocatable       :: wr,t,wm,tfit
  real(8),dimension(:),allocatable       :: tau


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


  ! !EQUILIUBRIUM (and Wigner transformed) GREEN'S FUNCTION 
  ! !=========================================================
  ! type(keldysh_equilibrium_gf)           :: gf0
  ! type(keldysh_equilibrium_gf)           :: gf
  ! type(keldysh_equilibrium_gf)           :: sf
  ! real(8),dimension(:),allocatable       :: exa


  !NON-EQUILIBRIUM FUNCTIONS:
  !=========================================================  
  !WEISS-FIELDS
  type(keldysh_contour_gf) :: G0
  !SELF-ENERGY
  type(keldysh_contour_gf) :: Sigma
  !LOCAL GF
  type(keldysh_contour_gf) :: locG
  !Bath SELF-ENERGY
  type(keldysh_contour_gf) :: S0



  !MOMENTUM-DISTRIBUTION
  !=========================================================  
  real(8),allocatable,dimension(:,:)     :: nk
  real(8),allocatable,dimension(:)       :: eq_nk


  !SUSCEPTIBILITY ARRAYS (in KADANOFF-BAYM)
  !=========================================================  
  real(8),allocatable,dimension(:,:,:,:) :: chi



  !NAMELISTS:
  !=========================================================
  namelist/variables/&
       dt,&
       beta,&
       Nstep        ,& 
       U            ,& 
       ts           ,& 
       eps          ,& 
       L            ,& 
       Lkreduced    ,& 
                                !DMFT
       nloop        ,& 
       eqnloop      ,& 
                                !BATH:
       bath_type    ,& 
       Vbath        ,& 
       Wbath        ,& 
       Walpha       ,&
       Wgap         ,&
                                !FIELD:
       Efield       ,& 
       field_type   ,& 
       Ex           ,& 
       Ey           ,& 
       t0           ,& 
       t1           ,& 
       Ncycles      ,& 
       omega0       ,& 
       E1           ,& 
                                !K-GRID
       Nx           ,& 
       Ny           ,& 
                                !CONVERGENCE:
       eps_error    ,& 
       nsuccess     ,& 
       weight       ,& 
                                !FLAGS:
       int_method   ,& 
       solve_eq     ,& 
       plot3D       ,& 
       fchi         ,& 
       fupdate      ,&
                                !FILES&DIR:
       irdSFILE      ,& 
       irdNkFILE




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

    !GLOBAL
    dt           = 0.1d0
    beta         = 10.d0
    Nstep        = 100
    U            = 4.d0
    ts           = 1.d0
    eps          = 0.01d0
    L            = 2048  
    Lkreduced    = 300
    !DMFT
    nloop        = 30
    eqnloop      = 50
    !BATH:
    bath_type    = 'flat'
    Vbath        = 0.d0
    Wbath        = 20.d0
    Walpha       = 1.d0
    Wgap         = 5.d0
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
    !K-GRID
    Nx           = 25
    Ny           = 25
    !CONVERGENCE:
    eps_error    = 1.d-4
    nsuccess     = 2
    weight       = 1.d0
    !FLAGS:
    int_method   = 'trapz'
    solve_eq     = .false. 
    plot3D       = .false.
    fchi         = .false.
    fupdate      = 0
    !FILES&DIR:
    irdSFILE      = 'restartSigma'
    irdNkFILE      = 'restartNk'

    inquire(file=adjustl(trim(inputFILE)),exist=control)
    if(control)then
       open(10,file=adjustl(trim(inputFILE)))
       read(10,nml=variables)
       close(10)
    else
       print*,"Can not find INPUT file"
       print*,"Dumping a default version in default."//trim(inputFILE)
       call dump_input_file("default.")
       call error("Can not find INPUT file, dumping a default version in default."//trim(inputFILE))
    endif

    !GLOBAL
    call parse_cmd_variable(dt           ,"DT")
    call parse_cmd_variable(beta         ,"BETA")
    call parse_cmd_variable(nstep        ,"NSTEP")
    call parse_cmd_variable(U            ,"U")
    call parse_cmd_variable(ts           ,"TS")
    call parse_cmd_variable(eps          ,"EPS")
    call parse_cmd_variable(L            ,"L")
    call parse_cmd_variable(Lkreduced    ,"LKREDUCED")
    !DMFT
    call parse_cmd_variable(nloop        ,"NLOOP")
    call parse_cmd_variable(eqnloop      ,"EQNLOOP")
    !BATH
    call parse_cmd_variable(Vbath        ,"VBATH")
    call parse_cmd_variable(bath_type    ,"BATH_TYPE")
    call parse_cmd_variable(wbath        ,"WBATH")
    call parse_cmd_variable(walpha       ,"WALPHA")
    call parse_cmd_variable(wgap         ,"WGAP")
    !EFIELD
    call parse_cmd_variable(field_type   ,"FIELD_TYPE")
    call parse_cmd_variable(Efield       ,"EFIELD")
    call parse_cmd_variable(Ex           ,"EX")
    call parse_cmd_variable(Ey           ,"EY")
    call parse_cmd_variable(t0           ,"T0")
    call parse_cmd_variable(t1           ,"T1")
    call parse_cmd_variable(ncycles      ,"NCYCLES")
    call parse_cmd_variable(omega0       ,"OMEGA0")
    call parse_cmd_variable(E1           ,"E1")
    !CONVERGENCE:
    call parse_cmd_variable(eps_error    ,"EPS_ERROR")
    call parse_cmd_variable(Nsuccess     ,"NSUCCESS")
    call parse_cmd_variable(weight       ,"WEIGHT")
    !GRID k-POINTS:
    call parse_cmd_variable(Nx           ,"NX")
    call parse_cmd_variable(Ny           ,"NY")
    !FLAGS:
    call parse_cmd_variable(int_method   ,"INT_METHOD")
    call parse_cmd_variable(solve_eq     ,"SOLVE_EQ")
    call parse_cmd_variable(plot3D       ,"PLOT3D")
    call parse_cmd_variable(fchi         ,"FCHI")
    call parse_cmd_variable(fupdate      ,"FUPDATE")
    !FILES&DIR:
    call parse_cmd_variable(irdSFILE      ,"IRDSFILE")
    call parse_cmd_variable(irdNkFILE     ,"IRDNKFILE")


    if(U==0.d0)Nloop=1
    Nfit=2*Nstep

    write(*,*)"CONTROL PARAMETERS"
    write(*,nml=variables)
    write(*,*)"--------------------------------------------"
    write(*,*)""
    call dump_input_file("used.")

  contains
    subroutine dump_input_file(prefix)
      character(len=*) :: prefix
      open(10,file=reg(prefix)//adjustl(trim(inputFILE)))
      write(10,nml=variables)
      close(10)
    end subroutine dump_input_file
  end subroutine read_input_init
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+----------------------------------------------------------------+
  !PURPOSE  : massive allocation of work array
  !+----------------------------------------------------------------+
  subroutine global_memory_allocation()
    integer          :: i
    real(8)          :: ex
    call msg("Allocating the arrays")
    !Weiss-fields:
    call allocate_keldysh_contour_gf(G0,Nstep*(Nstep+1)/2)
    !Interaction self-energies:
    call allocate_keldysh_contour_gf(Sigma,Nstep**2)
    !Local Green's functions:
    call allocate_keldysh_contour_gf(locG,Nstep*(Nstep+1)/2)
    !Bath self-energies:
    call allocate_keldysh_contour_gf(S0,Nfit**2)
    !Momentum-distribution:
    allocate(nk(Nstep,Lk),eq_nk(Lk))
    !Susceptibility/Optical response
    if(fchi)allocate(chi(2,2,Nstep,Nstep))
  end subroutine global_memory_allocation

  !******************************************************************
  !******************************************************************
  !******************************************************************


end module NEQ_VARS_GLOBAL

