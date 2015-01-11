module GLOC
  !########################################################################
  !     PROGRAM  : G0LOC
  !     TYPE     : Module
  !     PURPOSE  : Constructs or evaluate the G0(t,t') non-interacting GF. 
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !########################################################################
  use VARIABILI
  use FUNZIONI
  use MKL95_BLAS
  use MKL95_PRECISION

  implicit none
  include "mpif.h"
  private
  integer         :: i,j,k,k1,k2,ndim,id,mpiNT
  real(8)         :: AA,BB,I2,fug,rho,ee
  complex(8)      :: g,Chi



  public build_g0loc, extract_g0loc, get_g0epsiloc, get_gloc, ghat
  save

contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : build_g0loc
  !TYPE     : Subroutine
  !PURPOSE  : Construct the initial guess for G0_loc^c, c=t,<,>,tbar
  !+-------------------------------------------------------------------+
  subroutine build_g0loc(cc,g0c)
    complex(8),dimension(:,:) :: g0c
    character(len=1):: cc
    print*, 'Sono in BUILD_G0LOC'
    g0c=zero !set to zero the elements of the g0^c_loc matrix
    do i=1,N
       do j=1,N
          fug=exp(-xmu*(t(i)-t(j))) !exp(-mu(t-t'))
          I2=Bt(i,j)
          I2=I2**2
          I2=exp(-ts**2/4.d0*I2)!exp(-t*^2/4 I^2(t,t'))
          AA=At(i,j)
          do k=1,Nepsi
             ee=e(k)           !energy \e
             rho=dens1D(ee)    !\rho_o(\e)
             Chi=xi*Xc(cc,ee,i,j)!X^c(\e;t,t')
             g0c(i,j)=g0c(i,j) + de*rho*Chi*exp(-xi*ee*AA) !integrale su \e
          enddo
          g0c(i,j)=g0c(i,j)*fug*I2
       enddo
    enddo
    print*,'sto uscendo da BUILD_G0LOC'
  end subroutine build_g0loc




  !+-------------------------------------------------------------------+
  !PROGRAM  : extract_g0loc
  !TYPE     : Subroutine
  !PURPOSE  : Extract the blocks  G0_loc^c, c=t,<,>,tbar from \hat G0_loc
  !+-------------------------------------------------------------------+
  subroutine extract_g0loc(cc,g0,g0c)
    complex(8),dimension(:,:) :: g0c,g0
    character(len=1):: cc
    !get dim of g0loc = dim(g0hat)/2
    do i=1,N
       do j=1,N
          g=zero
          if (   cc=='t')then
             g= g0(i,j)

          elseif(cc=='<')then
             g=-g0(i,j+N)

          elseif(cc=='>')then
             g= g0(i+N,j)

          elseif(cc=='a')then
             g=-g0(i+N,j+N)
          endif
          g0c(i,j)=g
       enddo
    enddo
  end subroutine extract_g0loc


  !+-------------------------------------------------------------------+
  !PROGRAM  : g0epsiloc
  !TYPE     : Function
  !PURPOSE  : Construct the initial guess for G0_loc^c, c=t,<,>,tbar
  !+-------------------------------------------------------------------+
  subroutine get_g0epsiloc(e1,e2,g0c,M)
    integer :: M
    real(8) :: e1,e2
    complex(8),dimension(:,:) :: g0c
    character(len=1) :: cc

!    print*,'Entro in G0EPSILOC'

    do id=1,4      !loopa su char
       cc=vchar(id)
       do i=1,M    !loopa su t
          do j=1,M !loopa su t'
             fug=exp(-xmu*(t(i)-t(j))) !exp(-mu(t-t'))
             AA =e1*At(i,j)
             BB =e2*Bt(i,j)
             Chi=Xc(cc,e1,i,j)!X^c(\e;t,t')
             g=Chi*fug*exp(-xi*AA)*exp(-xi*BB)
             if(id == 1)then
                g0c(i,j)=g

             elseif(id == 2)then
                g0c(i,j+N)=-g

             elseif(id == 3)then
                g0c(i+N,j)=g

             elseif(id == 4)then
                g0c(i+N,j+N)=-g

             endif
          enddo
       enddo
    enddo
  end subroutine get_g0epsiloc
  
  !+-------------------------------------------------------------------+
  !PROGRAM  : GET_GLOC
  !TYPE     : Subroutine
  !PURPOSE  : Get the local GF \hat G_loc(t,t'), performing k-sum
  !+-------------------------------------------------------------------+
  subroutine get_gloc(g0c,M)
    integer :: M
    real(8) :: e1,e2
    complex(8),dimension(:,:) :: g0c
    complex(8),allocatable,dimension(:,:) :: g0epsiloc,glocal,glocal2,identity

    allocate(g0epsiloc(2*M,2*M),glocal(2*M,2*M),glocal2(2*M,2*M),identity(2*M,2*M))
    if(mpiMYID==0)print*,'Sono entrato in GET_GLOC'
    if(mpiMYID==1)print*,'Sono entrato in GET_GLOC'


    identity=zero
    do i=1,2*M
       identity(i,i)=one
    enddo
    !perhaps ID should be (\ ID, ID && ID, ID  \)
    
    g0c=zero
    !    do k1=1,Nepsi
    do k1=mpiMYID,Nepsi,mpiSIZE
       e1=e(k1)
!       if(mpiMYID==0)print*,k1
       do k2=1,Nepsi
          e2=e(k2)
          if(mpiMYID==0)print*,'get dens2D'
          rho=dens2D(e1,e2)
          if(mpiMYID==0)print*,'get g0epsiloc'
          call get_g0epsiloc(e1,e2,g0epsiloc,M)
          glocal=-identity
          if(mpiMYID==0)print*,'enters GEMM'
          call gemm(g0epsiloc,sighat,glocal) !G0loc*Sigma-\11
          glocal=-glocal                     !\11-G0loc*Sigma 
          if(mpiMYID==0)print*,'enters InvMat'
          call InvMat(glocal,2*N)            !(1-G0loc*Sigma)^-1
          glocal2=zero
          if(mpiMYID==0)print*,'enters GEMM'
          call gemm(glocal,g0epsiloc,glocal2)!{[(1-G0loc*Sigma)^-1]*G0loc}
          g0c=g0c+rho*glocal2*de**2          !rho(e1,e2)*{...}*de*de
       enddo
    enddo
  end subroutine get_gloc

end module GLOC
