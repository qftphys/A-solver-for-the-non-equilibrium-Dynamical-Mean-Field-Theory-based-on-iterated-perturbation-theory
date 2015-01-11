module FUNX_GLOBAL
  !###############################################################
  !     PROGRAM  : FUNCS_GLOBAL
  !     TYPE     : Module
  !     PURPOSE  : Constructs some functions used in other places. 
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !###############################################################
  use VARS_GLOBAL
  use FFTW
  use MKL95_LAPACK    !f95 wrapper for MKL_LAPACK routinse
  use MKL95_PRECISION !f95 wrapper for MKL_LAPACK routines

  implicit none
  private
  integer :: ndim
  real(8) :: dens1,dens2
  real(8) :: sgn,deltat
  public heaviside,&
       fermi,&
       dens_bethe,&
       dens_hyperc_1D,& 
       dens_hyperc_2D,&
       init_tgrid,&  
       init_egrid,&
       init_taugrid,&  
       extract, &
       getGmats, &
       shiftFW, &
       shiftBW, &
       plot_dislin,  &
       InvMat, &
       BuildLattice
  save

  interface InvMat
     module procedure R_InvMat, C_InvMat
  end interface

contains
  !+-----------------------------------------------------------------+
  !PROGRAM  : HEAVISIDE
  !TYPE     : function
  !PURPOSE  : calculate the Heaviside or step function
  !+-----------------------------------------------------------------+
  function heaviside(x)
    real(8) :: x,heaviside
    if(x < 0.d0) then
       heaviside = 0.0d0
    else
       Heaviside = 1.0d0
    endif
  end function heaviside



  !+-----------------------------------------------------------------+
  !PROGRAM  : FERMI
  !TYPE     : function
  !PURPOSE  : calculate the Fermi distribution function
  !+-----------------------------------------------------------------+
  function fermi(x)
    real(8) :: x,fermi
    fermi = 1.d0/(1.d0+exp(beta*x))
  end function fermi




  !+-----------------------------------------------------------------+
  !PROGRAM  : DENS_BETHE
  !TYPE     : function
  !PURPOSE  : calculate the non-interacting dos for BETHE lattice 
  !VERSION  : 31-01-2006
  !+-----------------------------------------------------------------+
  function dens_bethe(x)
    REAL(8) :: dens_bethe,x,t1,D
    complex(8):: root
    D=1.d0
    root=dcmplx((1.d0-1.d0*((x/D))**2),0.d0)
    root=sqrt(root)
    dens_bethe=(2.d0/(3.141592653589793238d0*D))*root
  end function dens_bethe




  !+-----------------------------------------------------------------+
  !PROGRAM  : DENS_HYPERC_1D
  !TYPE     : function
  !PURPOSE  : calculate non-interacting dos for 1d-HYPERCUBIC lattice
  !+-----------------------------------------------------------------+
  function dens_hyperc_1D(x)
    REAL(8):: dens_hyperc_1D,x,t1
    complex(8):: root
    t1=1.d0/sqrt(2.d0)
    dens_hyperc_1D= (1/(t1*sqrt(pi2)))*exp(-(x**2)/(2.*t1**2))
    return
  end function dens_hyperc_1D




  !+-----------------------------------------------------------------+
  !PROGRAM  : DENS_HYPERC_2D
  !TYPE     : function
  !PURPOSE  : calculate the non-interactin ghypercubic lattice dos
  !+-----------------------------------------------------------------+
  function dens_hyperc_2D(x1,x2)
    real(8) :: dens1,dens2,dens_hyperc_2D,x1,x2,t1,t2
    t1=sqrt2!1.d0/sqrt(2.d0)
    t2=sqrt2!1.d0/sqrt(2.d0)
    dens1 = (1.d0/(t1*sqrt(pi2)))*exp(-(x1**2)/(2.d0*t1**2))
    dens2 = (1.d0/(t2*sqrt(pi2)))*exp(-(x2**2)/(2.d0*t2**2))
    dens_hyperc_2D = dens1*dens2
  end function dens_hyperc_2D


  !+-----------------------------------------------------------------+
  !PROGRAM  : init_tgrid
  !TYPE     : Subroutine
  !PURPOSE  : initialize the real time axis grid of length given
  !+-----------------------------------------------------------------+
  subroutine init_tgrid(M)
    integer :: i,M
    do i=0,M
       t(i)=tmin+dble(i)*dt
    enddo
  end subroutine init_tgrid


  !+-----------------------------------------------------------------+
  !PROGRAM  : init_egrid
  !TYPE     : Subroutine
  !PURPOSE  : initialize the energy grid of length given
  !+-----------------------------------------------------------------+
  subroutine init_egrid(M)
    integer :: i,M
    de=(emax-emin)/dble(M)
    do i=0,M
       e(i)=emin+dble(i)*de
    enddo
  end subroutine init_egrid



  !+-----------------------------------------------------------------+
  !PROGRAM  : init_taugrid
  !TYPE     : Subroutine
  !PURPOSE  : initialize the real time axis grid of length given
  !+-----------------------------------------------------------------+
  subroutine init_taugrid(M,beta)
    integer :: i,M
    real(8) :: beta
    do i=0,M
       tau(i)=dble(i)*dtau
    enddo
  end subroutine init_taugrid



  !+-----------------------------------------------------------------+
  !PROGRAM  : EXTRACT
  !TYPE     : Subroutine
  !PURPOSE  : Sample a given function G(tau) over Nfak < N points.
  !COMMENTS : Incoming function is should have N+1 points (tau=0,beta)
  !this is modification with respect to the std extract routine used 
  !in HFqmc.
  !-g0 has N+1 points
  !-g00 will have Nfak+1 (tau_fak=0,beta)
  !+-----------------------------------------------------------------+
  subroutine extract(g0,g00)
    real(8),dimension(0:)  :: g0 !0:L
    real(8),dimension(0:) :: g00 !0:Lfak
    integer :: N,Nfak
    integer :: i,ip
    real(8) :: p,mismatch,beta

    N=size(g0)-1
    Nfak=size(g00)-1
    !Fix the end points
    g00(0)=g0(0)
    g00(Nfak)=1.d0-g0(0)
    mismatch=dble(N)/dble(Nfak)
    do i=1,Nfak-1
       p=dble(i)*mismatch
       ip=int(p)
       g00(i)=g0(ip)
    enddo
    return
  end subroutine extract



  !+-----------------------------------------------------------------+
  !PROGRAM  : getGmats
  !TYPE     : subroutine
  !PURPOSE  : calculate the Matsubara GF given the spectral density  
  !G^M(\iw)= int_\rrr d\e A(\e)/iw+\e=-\iw\int_\rrrd\e A(\e)/w^2+\e^2 
  !VERSION  : 10/2009
  !+-----------------------------------------------------------------+
  subroutine getGmats(gin,gout,beta)
    implicit none
    integer :: i,j
    real(8) :: w,wm,beta,A
    real(8) :: gmats_re,gmats_im
    complex(8),dimension(:) :: gin,gout
    if(size(gin) /= size(gout))stop "error in getGmats"
    do j=1,size(gin)
       wm=pi/beta*(2.d0*dble(j)-1.d0)
       gmats_re=zero;gmats_im=zero
       do i=1,n
          w=wmin+dble(i-1)*fmesh
          A=aimag(gin(i))/pi !DOS
          gmats_im=gmats_im + fmesh*A*wm/(wm**2+w**2)
          gmats_re=gmats_re + fmesh*A*w/(wm**2+w**2)
       enddo
       !gmats=-gmats*xi*wm
       gout(j)=gmats_re - xi*gmats_im
    enddo
  end subroutine getGmats




  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine shiftFW(Gin,Gout)
    complex(8),dimension(0:,0:) :: Gin
    complex(8),dimension(:,:)   :: Gout
    integer :: i,j,Ndim1,Ndim2
    Ndim1=size(Gout,1)
    Ndim2=size(Gout,2)
    do i=1,Ndim1
       do j=1,Ndim2
          Gout(i,j)=Gin(i-1,j-1)
       enddo
    enddo
  end subroutine shiftFW

  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine shiftBW(Gin,Gout)
    complex(8),dimension(0:,0:) :: Gout
    complex(8),dimension(:,:)   :: Gin
    integer :: i,j,Ndim1,Ndim2
    Ndim1=size(Gin,1)
    Ndim2=size(Gin,2)
    do i=1,Ndim1
       do j=1,Ndim2
          Gout(i-1,j-1)=Gin(i,j)
       enddo
    enddo
  end subroutine shiftBW


  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine plot_dislin(pname,xname,yname,zname, &
       X,Y,GF)
    character(len=160) :: cmd
    character(len=*) :: pname
    character(len=*) :: xname,yname,zname
    real(8),dimension(0:)    :: X,Y
    complex(8),dimension(0:,0:) :: GF

    !Intensity plot    
    call plot3D("3dPlot_Re"//trim(pname)//".png",&
         xname,yname,zname,real(X),real(Y),real(real(GF)))
    call plot3Dsurface("3dSurface_Re"//trim(pname)//".png",&
         xname,yname,zname,real(X),real(Y),real(real(GF)))

    call plot3D("3dPlot_Im"//trim(pname)//".png",&
         xname,yname,zname,real(X),real(Y),real(aimag(GF)))
    call plot3Dsurface("3dSurface_Im"//trim(pname)//".png",&
         xname,yname,zname,real(X),real(Y),real(aimag(GF)))

    print*,"Print "//trim(pname)
    cmd="rm -rf "//trim(pname)
    call system(trim(cmd))

    cmd="mkdir "//trim(pname)
    call system(trim(cmd))

    cmd="mv *"//trim(pname)//".png "//trim(pname)//"/"    
    call system(trim(cmd))
    print*,""

  end subroutine plot_dislin



  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine plot3D(pname,xname,yname,zname, &
       X,Y,GF)
    character(len=*) :: pname
    character(len=*) :: xname,yname,zname
    character(len=32):: xsuff,ysuff
    integer :: Nx,Ny,ixstep,iystep
    real(4) :: xstep,ystep
    real(4) :: xmin,xmax,dex
    real(4) :: ymin,ymax,dey
    real(4) :: zmin,zmax,dez
    real(4),dimension(0:)    :: X,Y
    real(4),dimension(0:,0:) :: GF

    Nx=size(GF,1)
    Ny=size(GF,2)

    xstep=real(dt);if(Nx == Ltau+1)xstep=1.
    ystep=real(dt);if(Ny == Ltau+1)ystep=1.
    X=X/xstep
    Y=Y/ystep

    xmin=minval(X);xmax=maxval(X)
    ixstep=nstep;if(Nx==Ltau+1)ixstep=int(beta)/2
    dex=abs(xmax-xmin); dex=dex/ixstep*5

    ymin=minval(Y);ymax=maxval(Y)
    iystep=nstep;if(Ny==Ltau+1)iystep=int(beta)/2
    dey=abs(ymax-ymin); dey=dey/iystep*5

    zmin=minval((GF));zmax=maxval((GF))
    GF=GF/abs(zmax-zmin) !normalization
    zmin=minval((GF))
    GF=GF+abs(zmin) !shift to 0

    xsuff="/$\Delta t$";if(Nx==Ltau+1)xsuff="/$\tau$"
    ysuff="/$\Delta t$";if(Ny==Ltau+1)ysuff="/$\tau$"


    CALL METAFL('PNG')
    CALL SETFIL(trim(pname))
    CALL SCRMOD('REVERSE')
    CALL WINSIZ(1920,1280)
    CALL ERRDEV('APPEND')
    CALL ERRFIL('dislin.out')
    CALL DISINI
    CALL ERRMOD('PROTOCOL','FILE')
    CALL PAGERA
    CALL HWFONT  
    CALL TITLIN(trim(pname),1)

    CALL TEXMOD('ON')
    CALL NAME(trim(xname)//trim(xsuff),'X')
    CALL NAME(trim(yname)//trim(ysuff),'Y')
    CALL NAME(trim(zname)//" ${\small\rm (normalized)}$",'Z')
    CALL LABDIG(3,'Z')
    CALL LABDIG(0,'XY')

    CALL INTAX
    CALL AUTRES(Nx,Ny)
    CALL SETVLT('RAIN') !TEMP
    CALL GRAF3(xmin,xmax,xmin,dex,  &
         ymin,ymax,ymin,dey,  &
         0.,1.05,0.,1.)
    CALL CRVMAT(GF,Nx,Ny,20,20)
    CALL TITLE

    CALL DISFIN
  end subroutine plot3D

  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine plot3Dsurface(pname,xname,yname,zname, &
       X,Y,GF)
    character(len=*) :: xname,yname,zname,pname
    character(len=7) :: char_temp
    character(len=32):: xsuff,ysuff
    integer :: Nx,Ny,ixstep,iystep,i
    real(4) :: xstep,ystep,angle
    real(4) :: xmin,xmax,dex
    real(4) :: ymin,ymax,dey
    real(4) :: zmin,zmax,dez
    real(4),dimension(0:) :: X,Y
    real(4),dimension(0:,0:) :: GF

    Nx=size(GF,1)
    Ny=size(GF,2)
    xstep=real(dt);if(Nx == Ltau+1)xstep=1.
    ystep=real(dt);if(Ny == Ltau+1)ystep=1.
    X=X/xstep
    Y=Y/ystep

    xmin=minval(X);xmax=maxval(X)
    ixstep=nstep;if(Nx==Ltau+1)ixstep=int(beta)/2
    dex=abs(xmax-xmin); dex=dex/ixstep*5

    ymin=minval(Y);ymax=maxval(Y)
    iystep=nstep;if(Ny==Ltau+1)iystep=int(beta)/2
    dey=abs(ymax-ymin); dey=dey/iystep*5

    zmin=minval((GF));zmax=maxval((GF))
    dez=max(abs(zmin),abs(zmax))    
    dez=dez/2.

    xsuff="/$\Delta t$";if(Nx==Ltau+1)xsuff="/$\tau$"
    ysuff="/$\Delta t$";if(Ny==Ltau+1)ysuff="/$\tau$"

    CALL METAFL('PNG')
    CALL SETFIL(trim(pname))
    CALL SCRMOD('REVERSE')
    CALL WINSIZ(1920,1280)

    CALL ERRDEV('APPEND')
    CALL ERRFIL('dislin.out')
    CALL DISINI
    CALL ERRMOD('PROTOCOL','FILE')

    CALL PAGERA
    CALL HWFONT  
    CALL TITLIN(trim(pname),2)

    CALL TEXMOD('ON')
    CALL NAME(trim(xname)//trim(xsuff),'X')
    CALL NAME(trim(yname)//trim(ysuff),'Y')
    CALL NAME(trim(zname),'Z')
    CALL LABDIG(6,'Z')
    CALL LABDIG(0,'XY')


    CALL VIEW3D(-3.75,-3.75,2.,'ABS')
    CALL SETVLT('RAIN') !TEMP
    CALL GRAF3D(xmin,xmax,xmin,dex,  &
         ymin,ymax,ymin,dey,  &
         zmin,zmax,zmin,dez)
    CALL SHDMOD("SMOOTH","SURFACE")
    CALL SURSHD(X,Nx,Y,Ny,GF)
    CALL SURMAT(GF,Nx,Ny,1,1)

    CALL TITLE
    CALL DISFIN

    return
  end subroutine plot3Dsurface

  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine plot3Dsurface_multi(pname,xname,yname,zname, &
       X,Y,GF)
    character(len=*) :: xname,yname,zname,pname
    character(len=7) :: char_temp
    character(len=32):: xsuff,ysuff
    integer :: Nx,Ny,ixstep,iystep,i
    real(4) :: xstep,ystep,angle
    real(4) :: xmin,xmax,dex
    real(4) :: ymin,ymax,dey
    real(4) :: zmin,zmax,dez
    real(4),dimension(0:) :: X,Y
    real(4),dimension(0:,0:) :: GF

    Nx=size(GF,1)
    Ny=size(GF,2)
    xstep=real(dt);if(Nx == Ltau+1)xstep=1.
    ystep=real(dt);if(Ny == Ltau+1)ystep=1.
    X=X/xstep
    Y=Y/ystep

    xmin=minval(X);xmax=maxval(X)
    ixstep=nstep;if(Nx==Ltau+1)ixstep=int(beta)/2
    dex=abs(xmax-xmin); dex=dex/ixstep*5

    ymin=minval(Y);ymax=maxval(Y)
    iystep=nstep;if(Ny==Ltau+1)iystep=int(beta)/2
    dey=abs(ymax-ymin); dey=dey/iystep*5

    zmin=minval((GF));zmax=maxval((GF))
    dez=max(abs(zmin),abs(zmax))    
    dez=dez/2.

    xsuff="/$\Delta t$";if(Nx==Ltau+1)xsuff="/$\tau$"
    ysuff="/$\Delta t$";if(Ny==Ltau+1)ysuff="/$\tau$"


    do i=1,10
       angle=360./10.*dble(i)
       print*,angle
       print*,adjustl(trim(char_temp)//trim(pname))
       write(char_temp,'(f7.2)')angle
       CALL METAFL('PNG')
       CALL SETFIL(trim(adjustl(trim(char_temp)//trim(pname))))
       CALL SCRMOD('REVERSE')
       CALL WINSIZ(1920,1280)
       CALL DISINI

       CALL PAGERA
       CALL HWFONT  
       CALL TITLIN(trim(pname),2)

       CALL TEXMOD('ON')
       CALL NAME(trim(xname)//trim(xsuff),'X')
       CALL NAME(trim(yname)//trim(ysuff),'Y')
       CALL NAME(trim(zname),'Z')
       CALL LABDIG(6,'Z')
       CALL LABDIG(0,'XY')


       CALL VIEW3D(angle,65.,6.,'ANGLE')
       CALL SETVLT('RAIN') !TEMP
       CALL GRAF3D(xmin,xmax,xmin,dex,  &
            ymin,ymax,ymin,dey,  &
            zmin,zmax,zmin,dez)
       CALL SHDMOD("SMOOTH","SURFACE")
       CALL SURSHD(X,Nx,Y,Ny,GF)
       CALL SURMAT(GF,Nx,Ny,0,0)

       CALL TITLE
       CALL DISFIN
    enddo
  end subroutine plot3Dsurface_multi




  !+-----------------------------------------------------------------+
  !PROGRAM  : R_InvMat(M,n)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a REAL*8 matrix M-->M^-1 using MKL_LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !+-----------------------------------------------------------------+
  subroutine R_InvMat(M,ndim)
    integer                      :: ndim
    real(8),dimension(ndim,ndim) :: M
    integer,dimension(ndim)      :: ipvt(ndim)
    call getrf(M,ipvt) !LU factorization
    call getri(M,ipvt) !Lapack Matrix Inversion
  end subroutine R_InvMat

  !+-----------------------------------------------------------------+
  !PROGRAM  : C_InvMat(M,n)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a COMPLEX*16 matrix M --> M^-1 using MKL_LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !+-----------------------------------------------------------------+
  subroutine C_InvMat(M,ndim)
    integer                         :: ndim
    complex(8),dimension(ndim,ndim) :: M
    integer,dimension(ndim)         :: ipvt(ndim)
    call getrf(M,ipvt) !LU factorization
    call getri(M,ipvt) !Lapack Matrix Inversion
  end subroutine C_InvMat



  !+-----------------------------------------------------------------+
  !PROGRAM  : LatticeStruct
  !TYPE     : Subroutine
  !PURPOSE  : Build the Lattice structure of the problem
  !COMMENT  : Square lattice for the time being
  !+-----------------------------------------------------------------+
  subroutine BuildLattice(Nx,Ny,t1,a)
    implicit none
    integer   :: Nx,Ny,i,j,it,iter
    real(8)   :: a,t1
    real(8)   :: c1,c2,dew,e
    real(8)   :: Kx,Ky,w,eps,Ed
    complex(8):: iw,gf

    eps   = 0.05d0
    dew=(wmax-wmin)/dble(L)

    Lk=(Nx+1)*(Ny+1)
    allocate(k(Lk,2),epsik(Lk),epsikt(Lk,0:L))

    !Get real lattice basis vectors 
    ai(1)=1.d0;ai(2)=0.d0;ai=a*ai
    aj(1)=0.d0;aj(2)=1.d0;aj=a*aj

    !Get reciprocal latice basis vectors
    bi(1)=1.d0;bi(2)=0.d0;bi=bi*pi2/a
    bj(1)=0.d0;bj(2)=1.d0;bj=bj*pi2/a

    call plot_realLat
    call plot_reciprocLat

    !Build the wave-vectors grid
    iter=0
    do i=0,Nx
       Kx=dble(i)/dble(Nx)
       do j=0,Ny
          Ky=dble(j)/dble(Ny)
          iter=iter+1
          k(iter,1)=Kx*bi(1)!Kx*i^^ !+Ky*bj(1)
          k(iter,2)=Ky*bj(2)!Ky^j^^ !Kx*bi(2)+Ky*bj(2)
       enddo
    enddo


    !Build dispersion relation arrays:
    !\epsilon(kx,ky)=epsk(i)
    !\epsilon(\ka - e\Ak(t))
    Ed=0.d0
    do i=1,Lk
       Kx=k(i,1)
       Ky=k(i,2)
       epsik(i)=-t1*(cos(Kx) + cos(Ky))            
       Ed=Ed+epsik(i)/dble(Lk)
       do it=0,L
          epsikt(i,it)=-t1*(cos(Kx+Efield*t(it))+cos(Ky+Efield*t(it)))
       enddo
    enddo
    call get_freeDOS
    call get_freeGtau

  contains 
    !=================================================================
    subroutine plot_realLat    
      implicit none
      open(30,file='RealLattice.dat')
      write(30,*)0,0
      write(30,*)''
      write(30,*)aj(1),aj(2)
      write(30,*)''
      write(30,*)ai(1),ai(2)
      write(30,*)''
      write(30,*)aj(1),-aj(2)
      write(30,*)''
      write(30,*)-aj(1),-aj(2)
      write(30,*)''
      write(30,*)-ai(1),ai(2)
      write(30,*)''
      write(30,*)-aj(1),aj(2)
      close(30)
      return
    end subroutine plot_realLat

    !=================================================================
    subroutine plot_reciprocLat
      implicit none
      open(40,file='ReciprocalLattice.dat')
      write(40,*)0,0
      write(40,*)''
      write(40,*)bj(1),bj(2)
      write(40,*)''
      write(40,*)bi(1)+bj(1),bi(2)+bj(2)
      write(40,*)''
      write(40,*)bi(1),bi(2)
      write(40,*)''
      write(40,*)-bj(1),-bj(2)
      write(40,*)''
      write(40,*)-bi(1),bi(2)
      write(40,*)''
      write(40,*)-bi(1),-bi(2)
      close(40)
    end subroutine plot_reciprocLat

    !=================================================================
    subroutine get_freeDOS
      implicit none
      open(70,file='DOSfree.dat')
      do j=1,L
         w=wmin+dble(j)*dew
         iw=cmplx(w,eps)
         gf=zero
         do i=1,Lk
            gf=gf+one/(iw+xmu-epsik(i))
         enddo
         gf=gf/dble(Lk)
         write(70,*)w,-1.d0/pi*aimag(gf)
         dens_lattice(j)=-aimag(gf)/pi
      enddo
    end subroutine get_freeDOS

    !=================================================================
    subroutine get_freeGtau
      implicit none
      complex(8),dimension(2*L) :: fg0iw
      real(8),dimension(0:L)    :: fg0tau
      fg0iw=zero
      do j=1,2*L
         w=pi/beta*dble(2*j-1)
         iw=cmplx(0.d0,w)
         gf=zero
         do i=1,Lk
            gf=gf+one/(iw+xmu-epsik(i))
         enddo
         gf=gf/dble(Lk)        
         fg0iw(j)=gf
      enddo
      call cfft_iw2it(fg0iw,fg0tau,beta)
      call extract(fg0tau,eqG00tau)
    end subroutine get_freeGtau
  end subroutine BuildLattice
end module FUNX_GLOBAL
