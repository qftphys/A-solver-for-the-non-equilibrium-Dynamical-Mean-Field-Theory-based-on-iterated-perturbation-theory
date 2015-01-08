MODULE NEQ_INPUT_VARS
  USE SCIFOR_VERSION
  USE PARSE_INPUT
  implicit none

  !GIT VERSION
  include "revision.inc"

  !Gloabl  variables
  !=========================================================
  integer                                :: Ntime         !Number of Time steps
  integer                                :: Ntau          !Imaginary time slices
  integer                                :: Niw           !Number of frequencies
  real(8)                                :: dt            !real-time step
  real(8)                                :: dtau          !imag-time step
  real(8)                                :: beta          !inverse temperature
  integer                                :: nloop         !dmft loop variables
  real(8)                                :: U0            !equilibrium local interaction
  real(8)                                :: Uloc          !non-equilibrium local interaction
  real(8)                                :: Vbath         !coupling to the Thermostat
  integer                                :: Lbath         !number of frequency in the bash DOS
  real(8)                                :: Wbath         !Width of the BATH DOS
  real(8)                                :: Walpha        !exponent of the pseudo-gapped bath.
  real(8)                                :: Wgap          !gap of the gapped bath
  character(len=16)                      :: bath_type     !choose the shape of the BATH
  real(8)                                :: eps           !broadening
  real(8)                                :: dmft_error     !convergence error threshold
  character(len=16)                      :: field_type    !choose the profile of the electric field
  integer                                :: Nsuccess      !number of convergence success

  !ELECTRIC FIELD VARIABLES (& NML):
  !=========================================================  
  real(8)                                :: Efield        !Electric field strength
  real(8)                                :: Ex,Ey         !Electric field vectors as input
  real(8)                                :: t0,t1         !turn on/off time, t0 also center of the pulse
  integer                                :: Ncycles       !Number of cycles in pulsed light packet
  real(8)                                :: omega0        !parameter for the Oscilatting field and Pulsed light
  real(8)                                :: E1            !Electric field strenght for the AC+DC case (tune to resonate)


contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input_init(inputFILE)
    character(len=*)               :: inputFILE

    !GLOBAL
    call parse_input_variable(Ntime      , "NTIME" , inputFILE , default      =100 , comment="Number of Real-Time steps")
    call parse_input_variable(Ntau       , "NTAU" , inputFILE , default       =50 , comment="Number of Imag-Time steps")
    call parse_input_variable(Niw        , "Niw" , inputFILE , default        =2048 , comment="Number of Matsubara Frequencies")
    call parse_input_variable(dt         , "DT" , inputFILE , default         =0.1d0 , comment="Real-time step")
    call parse_input_variable(beta       , "BETA" , inputFILE , default       =10d0 , comment="Inverse temperature")
    call parse_input_variable(U0         , "U0" , inputFILE , default         =2d0 , comment="equilibrium local interaction")
    call parse_input_variable(Uloc       , "ULOC" , inputFILE , default       =2d0 , comment="non-equilibrium local interaction")
    call parse_input_variable(eps        , "EPS" , inputFILE , default        =0.01d0 , comment="broadening")
    call parse_input_variable(nloop      , "NLOOP" , inputFILE , default      =30 , comment="Max number of DMFT loop")
    call parse_input_variable(dmft_error , "DMFT_ERROR" , inputFILE , default =1.d-3 , comment="DMFT convergence threshold")
    call parse_input_variable(Nsuccess   , "NSUCCESS" , inputFILE , default   =1 , comment="number of consecutive successes on convergence")
    !BATH
    call parse_input_variable(Vbath      , "VBATH" , inputFILE , default      =0d0 , comment="coupling to the thermostat")
    call parse_input_variable(Lbath      , "LBATH" , inputFILE , default      =1000 , comment="number of frequencies in the thermostat DOS")
    call parse_input_variable(bath_type  , "BATH_TYPE" , inputFILE , default  ='flat' , comment="thermostat DOS type [flat,gauss,bethe,...]")
    call parse_input_variable(wbath      , "WBATH" , inputFILE , default      =20d0 , comment="Width of the thermostat DOS")
    call parse_input_variable(walpha     , "WALPHA" , inputFILE , default     =1d0 , comment="exponent of the pseudo-gapped thermostat")
    call parse_input_variable(wgap       , "WGAP" , inputFILE , default       =5d0 , comment="gap of the gapped thermostat")
    !EFIELD
    call parse_input_variable(field_type , "FIELD_TYPE" , inputFILE , default ='dc' , comment="profile type of the electric field ")
    call parse_input_variable(Efield     , "EFIELD" , inputFILE , default     =0d0 , comment="electric field strength")
    call parse_input_variable(Ex         , "EX" , inputFILE , default         =1d0 , comment="electric field direction along x-axis")
    call parse_input_variable(Ey         , "EY" , inputFILE , default         =0d0 , comment="electric field direction along y-axis")
    call parse_input_variable(t0         , "T0" , inputFILE , default         =0d0 , comment="turn on time or center of the pulse")
    call parse_input_variable(t1         , "T1" , inputFILE , default         =10000d0 , comment="turn off time")
    call parse_input_variable(ncycles    , "NCYCLES" , inputFILE , default    =1 , comment="number of cycles in pulsed light signal ")
    call parse_input_variable(omega0     , "OMEGA0" , inputFILE , default     =acos(-1d0) , comment="parameter for the Oscilatting field and Pulsed light")
    call parse_input_variable(E1         , "E1" , inputFILE , default         =0d0 , comment="Electric field strenght for the AC+DC case (tune to resonate)")
    call save_input_file(inputFILE)
    call sf_version(revision)
  end subroutine read_input_init

end module NEQ_INPUT_VARS

