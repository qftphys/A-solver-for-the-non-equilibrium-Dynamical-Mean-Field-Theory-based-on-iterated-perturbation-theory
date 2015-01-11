MODULE VARS_GLOBAL
  !#####################################################################
  !     PROGRAM  : VARS_GLOBAL
  !     TYPE     : Module
  !     PURPOSE  : Defines the global variables used thru all the code
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !#####################################################################
  implicit none
  !number of time slices allowed first, number of energy points
  integer,parameter :: L=128, Lepsi=64, n=2*L, Ltau=Lepsi
  integer :: iloop,nloop,nstep,neqloop
  integer :: mpiIERR,mpiSIZE,mpiMYID

  !Parameters
  !===============================================================
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
  real(8),parameter :: wmin=-10.d0,wmax=-wmin
  real(8),parameter :: emin=-4.d0,emax=-emin
  real(8),parameter :: ts=1.d0/sqrt2 !hopping amplitude


  !Gloabl  variables
  real(8) :: Efield         !external Electric field          
  real(8) :: Vpd            !Hybridization between system and bath
  real(8) :: xmu            !chemical potential
  real(8) :: dt             !time step
  real(8) :: fmesh          !real frequency mesh 
  real(8) :: de             !energy step
  real(8) :: beta,temp,U    !inverse temp, temp, correlation
  real(8) :: dtau
  real(8),dimension(0:L)     :: t   !Time
  real(8),dimension(0:Lepsi) :: e   !energy
  real(8),dimension(0:Ltau)  :: tau !Im Time 

  !Total number of k-points
  integer :: Lk

  !Lattice (real, reciprocal) unitary vector (direction in the lattice space)
  real(8),dimension(2)               :: ai,aj
  real(8),dimension(2)               :: bi,bj

  !k-points
  real(8),allocatable,dimension(:,:) :: k

  !Dispersion relations (\e(k), \e(k-A(t))
  real(8),allocatable,dimension(:)   :: epsik
  real(8),allocatable,dimension(:,:) :: epsikt,Eplus,Eminus

  !Non-interacting DOS
  real(8),dimension(L)   :: dens_lattice

  !Gloval character variables
  character(len=1) :: char

  !Equiliubrium Green's function 
  complex(8),dimension(2*L)   :: eqG0w
  real(8),dimension(0:Ltau)   :: eqG00tau
  complex(8),dimension(2*L)   :: eqSw
  real(8),dimension(0:Ltau)   :: eqStau

  !non-equilibrium Green's function
  !non-interacting:
  complex(8),dimension(0:L,0:L)    :: G0gtr,ppG0gtr
  complex(8),dimension(0:L,0:L)    :: G0less,ppG0less
  complex(8),dimension(0:L,0:L)    :: G0ret,ppG0ret,G0adv,ppG0adv
  complex(8),dimension(0:L,0:Ltau) :: G0lceil,ppG0lceil
  complex(8),dimension(0:Ltau,0:L) :: G0rceil,ppG0rceil

  !self-energies
  complex(8),dimension(0:L,0:L)      :: Sgtr,Sless
  complex(8),dimension(0:L,0:L)      :: Sret,Sadv
  complex(8),dimension(0:L,0:Ltau)   :: Slceil
  complex(8),dimension(0:Ltau,0:L)   :: Srceil
  real(8),dimension(0:Ltau,0:Ltau)   :: Smatsubara


  !initial condition (killed after use at iloop=1)
  real(8),dimension(0:Ltau,0:Ltau)  :: Gmktau

  !local interacting Green's functions
  complex(8),allocatable,dimension(:,:) :: locGret,locGadv
  complex(8),allocatable,dimension(:,:) :: locGgtr,locGless
  complex(8),allocatable,dimension(:,:) :: locGlceil,locGrceil

  !saved initial conditions
  complex(8),allocatable,dimension(:)   :: icGkless
  complex(8),allocatable,dimension(:,:) :: icGklceil,icGkrceil

end MODULE VARS_GLOBAL
