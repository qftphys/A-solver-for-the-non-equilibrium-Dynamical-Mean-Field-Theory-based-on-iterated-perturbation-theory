MODULE VARIABILI
  !########################################################################
  !     PROGRAM  : VARIABILI
  !     TYPE     : Module
  !     PURPOSE  : Defines the global variables used thru all the code
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !########################################################################
  implicit none
  integer :: L,N
  integer :: Lepsi,Nepsi
  integer :: iloop,nloop
  integer :: mpiIERR,mpiSIZE,mpiMYID

  !Gloabl real variables
  real(8) :: Efield
  real(8) :: ts
  real(8) :: xmu
  real(8) :: dt,tmin
  real(8) :: de,emin,emax
  real(8) :: beta,temp,U
  real(8) :: fmesh

  !Global complex variables
  complex(8) :: D

  !Gloval character variables
  character(len=1) :: char

  !Arrays: time, energy, etc 
  real(8),allocatable,dimension(:) :: t,e 
  character(len=1),dimension(4) :: vchar = (/ 't', '<', '>', 'a' /)


  !Green's and Sigma functions
  complex(8),allocatable,dimension(:,:,:) :: g0loc
  complex(8),allocatable,dimension(:,:)   :: g0hat
  complex(8),allocatable,dimension(:,:)   :: fgloc,ghat
  complex(8),allocatable,dimension(:,:)   :: sighat  
!  complex(8),allocatable,dimension(:,:)   :: g0epsiloc,glocal,glocal2,identity
  !Parameters
  !==========================================================================
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi=(0.d0,1.d0)
  complex(8),parameter :: one=(1.d0,0.d0)
  real(8),parameter :: sqrt2 = 1.41421356237309504880169d0
  real(8),parameter :: sqrt3 = 1.73205080756887729352745d0
  real(8),parameter :: sqrt6 = 2.44948974278317809819728d0
  real(8),parameter :: pi=3.14159265358979d0
  real(8),parameter :: pi2=6.28318530717959d0

end MODULE VARIABILI
