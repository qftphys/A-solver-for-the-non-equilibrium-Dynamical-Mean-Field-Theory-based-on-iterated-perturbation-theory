PROGRAM freeNEQ
  !###################################################################
  !     PROGRAM  : freeNEQ
  !     TYPE     : Main program
  !     PURPOSE  : Solve system of free/interacting electrons driven 
  ! by electric field interacting with a resevoir.bath of free 
  ! electrons at temperature T
  !     AUTHORS  : Adriano Amaricci
  !###################################################################
  USE VARS_GLOBAL   
  USE FUNX_GLOBAL   
  USE FUNX_NEQ      !funzioni caratteristiche del problema al NEQ
  USE FFTW          !MKL Fast Fourier Transforms 
  USE EQ_SOLUTION   !routine per la soluzione del modell all'EQ
  USE KADANOFBAYM   !calcolo dellle funzioni di Green locali 

  implicit none
  integer :: i,j

  namelist/variables/beta,U,Efield,Vpd
  namelist/parameters/nloop,nstep,neqloop

  open(10,file='inputFILE.in')
  read(10,nml=variables)
  read(10,nml=parameters)
  close(10)

  xmu=0.d0
  temp=1.d0/beta
  fmesh=abs(wmax-wmin)/dble(n)
  dt=pi2/fmesh          
  dt=dt/dble(n)

  de=abs(emax-emin)/dble(Lepsi)
  dtau=beta/dble(Ltau)

  !Build time and energy arrays
  print*,"dt,dw=",dt,fmesh
  print*,"Inizializzo le griglie"
  print*,"t"
  call init_tgrid(L)
  print*,"e"
  call init_egrid(Lepsi)
  print*,"tau"
  call init_taugrid(Ltau,beta)
  print*,"Build Lattice"
  call BuildLattice(40,40,1.d0,1.d0)

  !Get free functions (no V no U) and the free analytic solution
  !WARNING: the analytic solution still suffer some buggy bugs.
  !call getG0analytic
  call getG0free

  iloop=1
  call update_sigma
  allocate(locGret(0:nstep,0:nstep),  &         
       locGadv(0:nstep,0:nstep),  &
       locGless(0:nstep,0:nstep), &
       locGgtr(0:nstep,0:nstep), &
       locGlceil(0:nstep,0:Ltau), &
       locGrceil(0:Ltau,0:nstep))

  call kadanof_baym2Gloc(nstep)

  deallocate(locGret,locGadv,locGgtr,locGless,locGlceil,locGrceil)
  print*,"CIAO"
end PROGRAM freeNEQ

