!#####################################################################
!     Program  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Defines the global variables used thru all the code
!     AUTHORS  : Adriano Amaricci
!#####################################################################
MODULE VARS_GLOBAL
  USE LATTICE
  implicit none
  !number of time slices allowed first, number of energy points
  integer,parameter :: L=2*512, n=2*L, Ltau=128, Lmu=1000
  integer           :: Lk

  !Parameters
  !===============================================================
  real(8),parameter :: ts=1.d0  !hopping amplitude
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi=(0.d0,1.d0)
  complex(8),parameter :: one=(1.d0,0.d0)
  complex(8),parameter :: D=one
  real(8),parameter :: sqrt2 = 1.41421356237309504880169d0
  real(8),parameter :: sqrt3 = 1.73205080756887729352745d0
  real(8),parameter :: sqrt6 = 2.44948974278317809819728d0
  real(8),parameter :: pi=3.14159265358979d0
  real(8),parameter :: pi2=6.28318530717959d0
  real(8),parameter :: tmin=0.d0
  logical,parameter :: TT=.true., FF=.false.

  !Gloabl  variables
  !=========================================================
  integer :: iloop,nloop,nstep,neqloop
  real(8) :: Efield         !Electric field strength
  type(vect2D) :: Ek        !Electric field vector
  real(8) :: Vpd            !Hybridization between system and bath
  real(8) :: xmu            !chemical potential
  real(8) :: dt,dtau        !time step
  real(8) :: fmesh,rfmesh   !real frequency mesh 
  real(8) :: de             !energy step
  real(8) :: beta,temp,U    !inverse temp, temp, correlation
  real(8) :: beta0,xmu0,U0  !quench variables        
  real(8) :: eps            !broadening
  !Grid arrays:
  !=========================================================  
  real(8),dimension(:),allocatable :: wr,wrmini,wm,t,tau,e
  real(8)                          :: wmin,wmax
  real(8)                          :: emin,emax  

  !Total number of k-points
  integer :: Nx,Ny
  integer,allocatable,dimension(:) :: sorted_ik
  real(8),allocatable,dimension(:) :: sorted_epsik


  !Equiliubrium Green's function 
  !=========================================================  
  complex(8),dimension(2*L)    :: eqG0w,eqSw
  real(8),dimension(0:Ltau)    :: eqStau

  !non-equilibrium Green's function: 4 = G^<,G^>,G^\lceil,G^\rceil
  !=========================================================  
  !NON-INTERACTING
  complex(8),allocatable,dimension(:,:) :: G0gtr
  complex(8),allocatable,dimension(:,:) :: G0less
  complex(8),allocatable,dimension(:,:) :: G0lceil
  complex(8),allocatable,dimension(:,:) :: G0rceil


  !SELF-ENERGIES
  !=========================================================  
  !Sigma^V_k,mu(t,t`)
  complex(8),allocatable,dimension(:)   :: S0gtr,S0less,S0ret
  complex(8),allocatable,dimension(:,:) :: S0lceil,S0rceil
  !Sigma^U(t,t`)
  complex(8),allocatable,dimension(:,:) :: Sgtr,Sless,Sret
  complex(8),allocatable,dimension(:,:) :: Slceil,Srceil
  real(8),allocatable,dimension(:,:)    :: Smatsubara


  !INTERACTING
  !=========================================================  
  complex(8),allocatable,dimension(:,:) :: locGgtr,locGless
  complex(8),allocatable,dimension(:,:) :: locGlceil,locGrceil

  !INITIAL CONDITIONS
  !=========================================================  
  !saved initial conditions
  complex(8),allocatable,dimension(:)   :: icGkless
  complex(8),allocatable,dimension(:,:) :: icGklceil,icGkrceil
  real(8),allocatable,dimension(:,:)    :: Gmktau

  !MPI vars:
  !=========================================================
  integer :: mpiERR,mpiSIZE,mpiID


  !NAMELISTS:
  !=========================================================  
  namelist/variables/beta,U,Efield,Vpd,beta0,xmu0,U0
  namelist/parameters/nloop,nstep,neqloop,Nx,Ny,wmin,emin,eps



contains


  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input(inputFILE) 
    character(len=*) :: inputFILE
    !variables:
    beta  = 100.0
    beta0 = 100.0
    U     = 6.0
    U0    = 6.0
    xmu0  = 0.0
    Efield= 0.0
    Vpd   = 0.0000001
    !parameters:
    nloop  = 30
    neqloop= 30
    nstep  = 20
    Nx     = 5
    Ny     = 5
    wmin   = -10.0
    emin   = -10.0
    eps    =  0.05d0
    open(10,file=adjustl(trim(inputFILE)))
    read(10,nml=variables)
    read(10,nml=parameters)
    close(10)
    if(nstep > L)then
       print*,"# # # # # # # # # # # # "
       print*,"ERROR: nstep > L"
       print*,"# # # # # # # # # # # # "
       stop
    endif
    wmax=-wmin;emax=-emin

    write(*,*)"CONTROL PARAMETERS"
    write(*,*)"--------------------------------------------"
    write(*,nml=variables)
    write(*,*)"--------------------------------------------"
    write(*,nml=parameters)
    write(*,*)"--------------------------------------------"
    print*,''
    return
  end subroutine read_input
  !******************************************************************
  !******************************************************************
  !******************************************************************


end module VARS_GLOBAL

