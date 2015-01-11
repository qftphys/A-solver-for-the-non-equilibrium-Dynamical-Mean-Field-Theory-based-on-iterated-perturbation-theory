  module GLOBAL_FUNCTIONS
  !########################################################################
  !     PROGRAM  : GLOBAL_FUNCTIONS
  !     TYPE     : Module
  !     PURPOSE  : Constructs some functions used in other places. 
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !########################################################################
  use GLOBAL_VARIABLES
  use MKL95_LAPACK    !f95 wrapper for MKL_LAPACK routinse
  use MKL95_PRECISION !f95 wrapper for MKL_LAPACK routines

  implicit none
  private
  integer :: i,j,k,ndim
  real(8) :: dens1,dens2
  real(8) :: sgn,deltat
  public heaviside,fermi,dens1D,dens2D,At,Bt,gather_blocks,InvMat,Xc

  save

  interface InvMat
     module procedure R_InvMat, C_InvMat
  end interface

contains
  !##################################################################
  !     FUNCTIONS
  !##################################################################
  !+-------------------------------------------------------------------+
  !PROGRAM  : HEAVISIDE
  !TYPE     : function
  !PURPOSE  : calculate the Heaviside or step function
  !+-------------------------------------------------------------------+
  function heaviside(x)
    real(8) :: x,heaviside
    if(x < 0.d0) then
       heaviside = 0.0d0
    elseif(x == 0.d0)then
       heaviside = 0.5d0
    elseif(x > 0.d0)then
       Heaviside = 1.0d0
    endif
  end function heaviside

  !+-------------------------------------------------------------------+
  !PROGRAM  : FERMI
  !TYPE     : function
  !PURPOSE  : calculate the Fermi distribution function
  !+-------------------------------------------------------------------+
  function fermi(x)
    real(8) :: x,fermi
    if(abs(x) <= 0.5d0)then
       fermi = 1.d0/(1.d0+exp(beta*(x)))      
    else
       if(x > 0.d0)fermi=0.d0
       if(x < 0.d0)fermi=1.d0
    endif
  end function fermi

  !+-------------------------------------------------------------------+
  !PROGRAM  : DENS_BETHE
  !TYPE     : function
  !PURPOSE  : calculate the non-interacting dos for BETHE lattice 
  !VERSION  : 31-01-2006
  !+-------------------------------------------------------------------+
  function dens_bethe(x)
    REAL(8) :: dens_bethe,x,t1,D
    complex(8):: root
    D=1.d0
    root=dcmplx((1.d0-1.d0*((x/D))**2),0.d0)
    root=sqrt(root)
    dens_bethe=(2.d0/(3.141592653589793238d0*D))*root
  end function dens_bethe
  
  !+-------------------------------------------------------------------+
  !PROGRAM  : DENS_HYPERC_1D
  !TYPE     : function
  !PURPOSE  : calculate the non-interacting dos for HYPERCUBIC lattice 1D
  !+-------------------------------------------------------------------+
  function dens_hyperc_1D(x)
    REAL(8):: dens_hyperc,x,t1,D
    complex(8):: root
    t1=1.d0/sqrt(2.d0)
    dens_hyperc = (1/(t1*sqrt(2*3.141592653589793238d0)))*exp(-(x**2)/(2*t1**2))
    return
  end function dens_hyperc_1D

  !+-------------------------------------------------------------------+
  !PROGRAM  : DENS_HYPERC_2D
  !TYPE     : function
  !PURPOSE  : calculate the non-interactin ghypercubic lattice dos
  !+-------------------------------------------------------------------+
  function dens_hyperc_2D(x1,x2)
    real(8) :: dens2D,x1,x2
    t1=1.d0/sqrt(2.d0)
    t2=1.d0/sqrt(2.d0)
    dens1 = (1.d0/(t1*sqrt(pi2)))*exp(-(x1**2)/(2.d0*t1**2))
    dens2 = (1.d0/(t2*sqrt(pi2)))*exp(-(x2**2)/(2.d0*t2**2))
    dens2D = dens1*dens2
  end function dens_hyperc_2D

  !+-------------------------------------------------------------------+
  !PROGRAM  : At
  !TYPE     : function
  !PURPOSE  : get A=int_t'^t cos(E*tau)*dtau
  !+-------------------------------------------------------------------+
  function At(i,j)
    integer :: i,j
    real(8) :: At
    At=0.d0
    if(i == j)then   !se sono uguali l'integrale e' zero
       At=0.d0

    elseif(j < i)then!se t' < t == t(j) < t(i) non cambi segno
       sgn=1.d0
       do k=j,i
          At=At+cos(Efield*t(k))*dt
       enddo

    elseif(i < j)then!se t' > t == t(j) > t(i) cambi segno
       sgn=-1.d0
       do k=i,j
          At=At+cos(Efield*t(k))*dt
       enddo
    endif
    At=At*sgn
  end function At

  !+-------------------------------------------------------------------+
  !PROGRAM  : Bt
  !TYPE     : function
  !PURPOSE  : get B=int_t'^t sin(E*tau)*dtau
  !+-------------------------------------------------------------------+
  function Bt(i,j)
    integer :: i,j
    real(8) :: Bt
    Bt=0.d0
    if(i == j)then   !se sono uguali l'integrale e' zero 
       Bt=0.d0

    elseif(j < i)then!se t' < t == t(j) < t(i) non cambi segno
       sgn=1.d0
       do k=j,i
          Bt=Bt+sin(Efield*t(k))*dt
       enddo

    elseif(i < j)then!se t' > t == t(j) > t(i) cambi segno
       sgn=-1.d0
       do k=i,j
          Bt=Bt+sin(Efield*t(k))*dt
       enddo
    endif
    Bt=Bt*sgn
  end function Bt

  !+-------------------------------------------------------------------+
  !PROGRAM  : Xc
  !TYPE     : function
  !PURPOSE  : build  the function X^c(e,t,t') equal to the term in front
  !           of the exponential in the expressions for G_0(t,t') 
  !COMMENT  : The function is defined WITHOUT the imaginary unit xi
  !+-------------------------------------------------------------------+
  function Xc(cc,e1,i,j)
    integer :: i,j
    real(8) :: Xc,e1
    character(len=1) :: cc
    Xc=0.d0
    if(cc == 't')then
       deltat=t(i)-t(j) !deltat=t-t'
       Xc=(fermi(e1-xmu) - heaviside(deltat))
    elseif(cc == '<')then
       Xc=fermi(e1-xmu)
    elseif(cc == '>')then
       Xc=(fermi(e1-xmu)-1.d0)
    elseif(cc == 'a')then
       deltat=t(j)-t(i) !tt=t'-t
       Xc=(fermi(e1-xmu)- heaviside(deltat))
    endif
  end function Xc



  !+-------------------------------------------------------------------+
  !PROGRAM  : gather_block
  !TYPE     : Subroutine
  !PURPOSE  : Build the Keldysh matrix (\ S/G_t, -S/G_< && 
  !                                       S/G_>, -S/G_tbar \)
  !           for the Sigma/Green's function, gathering the blocks.
  !+-------------------------------------------------------------------+
  subroutine gather_blocks(floc,f)
    complex(8),dimension(:,:,:) :: floc
    complex(8),dimension(:,:)   :: f
    print*,'Sono in GATHER_BLOCKS'
    ndim=size(floc,2)
    print*,ndim
    do i=1,ndim
       do j=1,ndim
          f(i,j)          = floc(1,i,j)!t
          f(i,j+ndim)     =-floc(2,i,j)!<
          f(i+ndim,j)     = floc(3,i,j)!>
          f(i+ndim,j+ndim)=-floc(4,i,j)!tbar
       enddo
    enddo
    print*,'sto uscendo da GATHER_BLOCKS'
  end subroutine gather_blocks

  !+-------------------------------------------------------------------+
  !PROGRAM  : R_InvMat(M,n)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a REAL*8 matrix M --> M^-1 using MKL_LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !+-------------------------------------------------------------------+
  subroutine R_InvMat(M,ndim)
    integer                      :: ndim
    real(8),dimension(ndim,ndim) :: M
    integer,dimension(ndim)      :: ipvt(ndim)
    call getrf(M,ipvt) !LU factorization
    call getri(M,ipvt) !Lapack Matrix Inversion
  end subroutine R_InvMat

  !+-------------------------------------------------------------------+
  !PROGRAM  : C_InvMat(M,n)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a COMPLEX*16 matrix M --> M^-1 using MKL_LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !+-------------------------------------------------------------------+
  subroutine C_InvMat(M,ndim)
    integer                         :: ndim
    complex(8),dimension(ndim,ndim) :: M
    integer,dimension(ndim)         :: ipvt(ndim)
    call getrf(M,ipvt) !LU factorization
    call getri(M,ipvt) !Lapack Matrix Inversion
  end subroutine C_InvMat

end module FUNZIONI
