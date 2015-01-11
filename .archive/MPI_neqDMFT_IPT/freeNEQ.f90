Program freeNEQ
  !###################################################################
  !PROGRAM  : freeNEQ
  !TYPE     : Main program
  !PURPOSE  : Solve system of free/interacting electrons driven 
  ! by electric field interacting with a resevoir/bath of free 
  ! electrons at temperature T
  !AUTHORS  : Adriano Amaricci
  !###################################################################
  !LIBRARY:  
  USE TOOLS
  USE GRIDS
  USE LATTICE
  USE DLPLOT
  USE MPI
  !LOCAL:
  USE VARS_GLOBAL 
  USE EQKELDYSH   
  USE FUNX_NEQ    
  USE KADANOFBAYM 
  implicit none

  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  print*,'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)

  !##################################################################
  !PREAMBLE TO THE CALCULATION:
  !Read inputs & set the grids:
  !comments: if _mpi this is done by everybody
  !===============================================================
  call read_input("inputFILE.in")
  xmu   = 0.d0  ;  temp=1.d0/beta
  fmesh = abs(wmax-wmin)/dble(2*L)
  rfmesh= abs(wmax-wmin)/dble(2*nstep)
  dt    = pi2/fmesh ; dt = dt/dble(2*L)
  de    = abs(emax-emin)/dble(Lmu)
  dtau  = beta/dble(Ltau)
  if(mpiID==0)then
     write(*,'(A,F9.6)')"dt=",dt
     write(*,'(A,F9.6)')"dw=",fmesh
     call dump("")
     call dump("Init grids:")
  end if
  allocate(wr(2*L),wrmini(2*nstep),wm(2*L),t(-L:L),e(0:Lmu),tau(0:Ltau))
  call init_wgrid(wr,wmin,fmesh)
  call init_wgrid(wrmini,wmin,rfmesh)
  call init_tgrid(t,dt,L)
  call init_wmgrid(wm,beta)
  call init_taugrid(tau,-dtau)
  call init_egrid(e,emin,de)
  if(mpiID==0)call dump("")

  !Build Lattice structure:
  !===============================================================
  if(mpiID==0)call dump("Get Lattice Structure:")
  call BuildLattice(Nx,Ny,Lk)
  call get_epsik(ts,Lk)
  call get_epsimu(emin,Lmu)
  Ek%x=1.d0;Ek%y=1.d0;Ek=(Efield*pi2/alat)*Ek!Define Electric field vector
  call get_epsikt(Lk,Ek,t(0:L))
  allocate(sorted_ik(Lk),sorted_epsik(Lk))
  sorted_epsik=epsik
  call sort_array(sorted_epsik,sorted_ik)
  if(mpiID==0)call dump("")


  !##################################################################
  !STARTS THE REAL WORK:
  !Solve the associated equilibrium problem:
  !===============================================================
  if(mpiID==0)call keldysheq(eqSw,eqG0w)
  call MPI_BCAST(eqSw,2*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  call MPI_BCAST(eqG0w,2*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)

  !Massive allocation of the working arrays:
  !===============================================================
  call allocate_locFunx('a')

  !Get the WF guess, Sigma guess, bath GF:
  !===============================================================
  call getG0guess
  call Build_bath_t
  call update_sigma

  !neq-DMFT iterative solution (NEARLY WORKING) 
  !===============================================================
  do iloop=1,nloop
     call kadanof_baym2Gloc()
     call dysonG0()
     call update_sigma
  enddo

  !Massive deallocation of the working arrays:
  !===============================================================
  call allocate_locFunx('d')
  CALL MPI_FINALIZE(mpiERR)
  print*,"BRAVO!"
end PROGRAM freeNEQ
!********************************************************************
!********************************************************************
!********************************************************************
