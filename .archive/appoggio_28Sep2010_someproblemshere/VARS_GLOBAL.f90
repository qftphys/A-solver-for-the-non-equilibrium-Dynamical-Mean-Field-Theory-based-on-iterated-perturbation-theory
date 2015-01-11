!#####################################################################
!     Program  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Defines the global variables used thru all the code
!     AUTHORS  : Adriano Amaricci
!#####################################################################
MODULE VARS_GLOBAL
  !LIBRARY:
  USE LATTICE
  USE TOOLS
  USE GRIDS
  implicit none
  integer,parameter :: Lmu=256  !# of bath energies
  integer,parameter :: Ltau=32  !
  integer,protected :: L        !a big number
  integer           :: Lk       !# of k-points

  !Parameters
  !===============================================================
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi=(0.d0,1.d0)
  complex(8),parameter :: one=(1.d0,0.d0)
  complex(8),parameter :: D=one
  real(8),parameter :: tiny=1.d-15
  real(8),parameter :: sqrt2 = 1.41421356237309504880169d0
  real(8),parameter :: sqrt3 = 1.73205080756887729352745d0
  real(8),parameter :: sqrt6 = 2.44948974278317809819728d0
  real(8),parameter :: pi=3.14159265358979d0
  real(8),parameter :: pi2=6.28318530717959d0
  real(8),parameter :: tmin=0.d0
  logical,parameter :: TT=.true., FF=.false.

  !Gloabl  variables
  !=========================================================
  integer      :: iloop,nloop    !DMFT loop variables
  integer      :: nstep          !Number of Time steps
  integer      :: eqnloop        !DMFT loop for equilibrium solution
  real(8)      :: Efield         !Electric field strength
  type(vect2D) :: Ek             !Electric field vector
  real(8)      :: Vpd            !Hybridization between system and bath
  real(8)      :: ts             !hopping amplitude
  real(8)      :: xmu            !chemical potential
  real(8)      :: dt,dtau        !time step
  real(8)      :: fmesh          !real frequency mesh 
  real(8)      :: de             !energy step
  real(8)      :: beta,temp,U    !inverse temp, temp, correlation
  real(8)      :: beta0,xmu0,U0  !quench variables        
  real(8)      :: eps            !broadening
  logical      :: iquench        !quench flag
  logical      :: eqflag,wfftw   !flag to fine tune the calculation
  !                              !eqflag=T work out associated equilibrium sol
  !                              !wfftw=T solve with fftw if Efield=0
  character(len=6):: method      !choose the perturbation theory method: IPT,SPT
  logical      :: plot3D         !flag to print or not 3D functions

  !Grid arrays:
  !=========================================================  
  real(8),dimension(:),allocatable :: wr,wm,t,tau,e
  real(8)                          :: wmin,wmax
  real(8)                          :: emin,emax  

  !Total number of k-points
  integer :: Nx,Ny

  !Dispersion relations arrays:
  !=========================================================  
  real(8),dimension(:),allocatable   :: epsimu,epsik,sorted_epsik
  integer,dimension(:),allocatable   :: sorted_ik              
  real(8),dimension(:,:),allocatable :: epsikt

  !EQUILIUBRIUM/WIGNER TRANSFORMED GREEN'S FUNCTION 
  !=========================================================  
  !correlated equilibrium initial conditions:
  real(8),allocatable,dimension(:)     :: icnk
  complex(8),allocatable,dimension(:)  :: icSiw,icGiw
  complex(8),allocatable,dimension(:)  :: icG0w,icG0tless,icG0tgtr
  real(8),allocatable,dimension(:)     :: icGtau,icStau

  !Frequency domain:
  complex(8),dimension(:),allocatable :: g0fret,g0fless,g0fgtr
  complex(8),dimension(:),allocatable :: gfret,gfless,gfgtr
  complex(8),dimension(:),allocatable :: sfret
  !Time domain:
  complex(8),dimension(:),allocatable :: g0tret,g0tless,g0tgtr
  complex(8),dimension(:),allocatable :: gtret,gtless,gtgtr
  complex(8),dimension(:),allocatable :: stret,stless,stgtr

  real(8),dimension(:),allocatable    :: exa

  !NON-EQUILIBRIUM GREEN'S FUNCTION: 4 = G^<,G^>,G^\lceil,G^\rceil
  !=========================================================  
  !NON-INTERACTING
  complex(8),allocatable,dimension(:,:) :: G0gtr
  complex(8),allocatable,dimension(:,:) :: G0less
  complex(8),allocatable,dimension(:,:) :: G0lceil
  complex(8),allocatable,dimension(:,:) :: G0rceil


  !SELF-ENERGIES
  !=========================================================  
  !Sigma^V_k,mu(t,t`)
  complex(8),allocatable,dimension(:)   :: S0gtr,S0less
  complex(8),allocatable,dimension(:,:) :: S0lceil,S0rceil
  !Sigma^U(t,t`)
  complex(8),allocatable,dimension(:,:) :: Sgtr,Sless
  complex(8),allocatable,dimension(:,:) :: Slceil,Srceil
  complex(8),allocatable,dimension(:,:) :: Smatsubara

  !INTERACTING
  !=========================================================  
  complex(8),allocatable,dimension(:,:) :: locGgtr,locGless
  complex(8),allocatable,dimension(:,:) :: locGlceil,locGrceil

  !INITIAL CONDITIONS
  !=========================================================  
  complex(8),allocatable,dimension(:)   :: icGkless
  complex(8),allocatable,dimension(:,:) :: icGktau

  !COMMON ARRAYS 4 KADANOFFBAYM module
  !=========================================================  
  complex(8),allocatable,dimension(:,:)   :: Gkless,Gkgtr
  complex(8),allocatable,dimension(:)     :: Ikless,Ikgtr
  complex(8),allocatable,dimension(:)     :: Ikless0,Ikgtr0
  complex(8),allocatable,dimension(:,:)   :: Gklceil,Gkrceil
  complex(8),allocatable,dimension(:)     :: Gktau
  complex(8),allocatable,dimension(:)     :: Iklceil,Ikrceil
  complex(8),allocatable,dimension(:)     :: Iklceil0,Ikrceil0
  complex(8)                              :: Ikdiag
  complex(8),allocatable,dimension(:,:)   :: Udelta,Vdelta
  real(8),allocatable,dimension(:,:)      :: nk




  !MPI vars:
  !=========================================================
  integer :: mpiERR,mpiSIZE,mpiID

  !NAMELISTS:
  !=========================================================  
  namelist/variables/dt,beta,U,Efield,Vpd,ts,nstep,nloop,eqnloop
  namelist/quench/iquench,beta0,xmu0,U0
  namelist/latticeN/Nx,Ny
  namelist/parameters/L,emin,eps,eqflag,method,wfftw,plot3D



contains


  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input(inputFILE) 
    character(len=*) :: inputFILE
    !variables:
    dt     = 0.157080
    beta   = 100.0
    U      = 6.0
    Efield = 0.0
    Vpd    = 0.0
    ts     = 1.0
    nstep  = 50
    nloop  = 30
    eqnloop= 20
    !LatticeN
    Nx     = 50
    Ny     = 50    
    !parameters:
    L      = 1024
    emin   = -10.0
    eps    =  0.05d0
    eqflag = .false.
    method = 'ipt'
    wfftw  = .false.
    plot3D = .true.
    !quench
    iquench= .false.
    beta0  = 100.0
    U0     = 6.0
    xmu0   = 0.0    
    open(10,file=adjustl(trim(inputFILE)))
    read(10,nml=variables)
    read(10,nml=LatticeN)
    read(10,nml=parameters)
    read(10,nml=quench)
    close(10)
    emax=-emin
    write(*,*)"CONTROL PARAMETERS"
    write(*,*)"--------------------------------------------"
    write(*,nml=variables)
    write(*,*)"--------------------------------------------"
    write(*,nml=quench)
    write(*,*)"--------------------------------------------"
    write(*,nml=LatticeN)
    write(*,*)"--------------------------------------------"
    write(*,nml=parameters)
    write(*,*)"--------------------------------------------"
    print*,''
    return
  end subroutine read_input
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+----------------------------------------------------------------+
  !PROGRAM  : INIT_CALC
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine init_calc
    xmu   = 0.d0  ;  temp=1.d0/beta
    wmax   = pi/dt ; wmin=-wmax    
    fmesh  = wmax/dble(nstep)
    de     = abs(emax-emin)/dble(Lmu)
    dtau   = beta/dble(Ltau)
    write(*,'(A,F12.6)')"dt   =",dt
    write(*,'(A,F12.6)')"dw   =",fmesh
    write(*,'(A,F12.6)')"wmax =",wmax
    call dump("")
    call dump("Init grids:")
    allocate(wr(2*nstep),wm(2*nstep),t(-nstep:nstep),e(0:Lmu),tau(0:Ltau))
    call init_wgrid(wr,wmin,fmesh)
    call init_wmgrid(wm,beta)
    call init_tgrid(t,dt,nstep)
    call init_taugrid(tau,dtau)
    call init_egrid(e,emin,de)
    call dump("")
  end subroutine init_calc
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+----------------------------------------------------------------+
  !PROGRAM  : MASSIVE_ALLOC
  !TYPE     : subroutine
  !PURPOSE  : massive allocation of work array
  !+----------------------------------------------------------------+
  subroutine alloc_memory(char)
    character(len=1) :: char
    integer          :: i
    real(8)          :: ex
    if(char=='a')then
       call dump("Allocating the memory:")
       allocate(G0gtr(0:nstep,0:nstep),G0less(0:nstep,0:nstep))
       allocate(Gkless(0:nstep,0:nstep),Gkgtr(0:nstep,0:nstep))

       allocate(S0gtr(-nstep:nstep),S0less(-nstep:nstep))
       allocate(Sgtr(0:nstep,0:nstep),Sless(0:nstep,0:nstep))

       allocate(locGless(0:nstep,0:nstep),locGgtr(0:nstep,0:nstep))

       allocate(Ikless(0:nstep),Ikgtr(0:nstep))
       allocate(Ikless0(0:nstep),Ikgtr0(0:nstep))

       allocate(icGkless(Lk))

#ifdef _mix
       allocate(G0lceil(0:nstep,0:Ltau),G0rceil(0:Ltau,0:nstep))
       allocate(Gklceil(0:nstep,0:Ltau),Gkrceil(0:Ltau,0:nstep))
       allocate(Gktau(-Ltau:Ltau))

       allocate(S0lceil(0:nstep,0:Ltau),S0rceil(0:Ltau,0:nstep))
       allocate(Slceil(0:nstep,0:Ltau),Srceil(0:Ltau,0:nstep)) 

       allocate(locGlceil(0:nstep,0:Ltau),locGrceil(0:Ltau,0:nstep))

       allocate(Iklceil(0:Ltau))
       allocate(Iklceil0(0:Ltau))

       allocate(icGktau(Lk,0:Ltau))
#endif

       allocate(Udelta(Lk,0:nstep),Vdelta(Lk,0:nstep))
       allocate(nk(0:nstep,Lk),icnk(Lk))

       allocate(g0fret(2*nstep),g0fless(2*nstep),g0fgtr(2*nstep))
       allocate(gfret(2*nstep),gfless(2*nstep),gfgtr(2*nstep))
       allocate(sfret(2*nstep))

       allocate(g0tret(-nstep:nstep),g0tless(-nstep:nstep),g0tgtr(-nstep:nstep))
       allocate(gtret(-nstep:nstep),gtless(-nstep:nstep),gtgtr(-nstep:nstep))
       allocate(stret(-nstep:nstep),stless(-nstep:nstep),stgtr(-nstep:nstep))

       allocate(icG0tless(-L:L),icG0tgtr(-L:L),icSiw(2*L),icGiw(2*L))
       allocate(icG0w(2*L))
       allocate(icGtau(-L:L),icStau(-L:L))

       allocate(exa(-nstep:nstep))
       ex=-1.d0       
       do i=-nstep,nstep
          ex=-ex
          exa(i)=ex
       enddo

       call dump("done")
       call dump("")

    elseif(char=='d')then
       call dump("Deallocating:")
       deallocate(G0gtr,G0less)
       deallocate(Gkless,Gkgtr)      
       deallocate(S0less,S0gtr)
       deallocate(Sgtr,Sless)
       deallocate(locGless,locGgtr)
       deallocate(Ikless,Ikgtr)
       deallocate(Ikless0,Ikgtr0)
       deallocate(icGkless)
#ifdef _mix
       deallocate(G0lceil,G0rceil)
       deallocate(Gklceil,Gkrceil,Gktau)       
       deallocate(S0lceil,S0rceil)
       deallocate(Slceil,Srceil)
       deallocate(locGlceil,locGrceil)
       deallocate(Iklceil,Iklceil0)
       deallocate(icGktau)
#endif

       deallocate(Udelta,Vdelta)
       deallocate(nk,icnk)

       deallocate(g0fret,g0fless,g0fgtr)
       deallocate(gfret,gfless,gfgtr)
       deallocate(sfret)
       deallocate(g0tret,g0tless,g0tgtr)
       deallocate(gtret,gtless,gtgtr)
       deallocate(stret,stless,stgtr)

       deallocate(icG0tless,icG0tgtr,icSiw,icGiw)
       deallocate(icG0w)
       deallocate(icGtau,icStau)
       deallocate(exa)
       call dump("done")
       call dump("")

    end if
  end subroutine alloc_memory
  !******************************************************************
  !******************************************************************
  !******************************************************************







  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !Build BATH dispersion arrays:
  !\epsilon_bath(\ka)          =epsimu(i)
  !+-----------------------------------------------------------------+
  subroutine get_epsimu(ebini,epsimu_)
    real(8),dimension(:),intent(inout) :: epsimu_
    real(8),intent(in)                 :: ebini
    integer                            :: ik,Lmu_
    real(8)                            :: ebfin,en
    Lmu_=size(epsimu_)
    ebfin=-ebini
    do ik=1,Lmu_     
       en=ebini + dble(ik)*abs(ebfin-ebini)/dble(Lmu_)
       epsimu_(ik)= en
    enddo
    call get_DOS(epsimu_,'DOSbath.lattice')
    call system("if [ ! -d LATTICEinfo ]; then mkdir LATTICEinfo; fi")
    call system("mv *.lattice LATTICEinfo/ 2>/dev/null")
  end subroutine get_epsimu
  !***************************************************************
  !***************************************************************
  !***************************************************************

end module VARS_GLOBAL

