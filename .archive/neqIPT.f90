PROGRAM neqIPT
  !###################################################################
  !     PROGRAM  : neqIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model with applied external 
  ! electric field E, on the hypercubic lattice, using neq. DMFT 
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !###################################################################
  USE VARS_GLOBAL   
  USE FUNX_GLOBAL   
  USE FUNX_NEQ      !funzioni caratteristiche del problema al NEQ
  USE FFTW          !MKL Fast Fourier Transforms 
  USE EQ_SOLUTION   !routine per la soluzione del modell all'EQ
  USE KADANOFBAYM   !calcolo dellle funzioni di Green locali 


  implicit none
  integer :: i,j

  namelist/variables/beta,U,Efield
  namelist/parameters/nloop,nstep,neqloop

  open(10,file='inputFILE.in')
  read(10,nml=variables)
  read(10,nml=parameters)
  close(10)

  xmu=0.d0 !half-filling
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

  !Solve the model at equilibrium to get initial condition
  !exported functions are:
  !- (eqSw)     Sigma(w)  : real freq. self-energy     
  !- (eqStau)   Sigma(tau): imaginary time self-energy
  !- (eqG0w)    G_0(w)    : real freq. bath function
  !- (eqG00tau) G_0(tau)  : imaginary time bath function
  call EQ_HMKELDYSH(neqloop,eqSw,eqStau,eqG0w,eqG00tau)

  !Obtain a "guess" for the local Weiss field \calG_0(t,t')
  !=\int d\e X^c*(\e) eq\calG_0(\e)
  !and all the related local function (>,<,R,\lceil)
  !This routine will build:
  !- G0gtr
  !- G0less
  !- G0ret 
  !- G0lceil
  print*,"Calcolo G0guess(t1,t2)"
  print*,"======================"
  print*,""
  call getG0guess()

  !Evaluate the first self-energies
  !this routine provides:
  !- Sret
  !- Sgtr
  !- Sless
  !- Slceil
  !used in the following to get local GF (see kadanof_baym2Gloc)
  print*,"Calcolo la nuova SIGMA"
  print*,"======================"
  print*,""
  call update_sigma()

  !Starts the DMFT iteration
  do iloop=1,nloop
     print*, 'iloop',iloop,'/',nloop

     !Build the local Green's function thru summmation over k
     !of the solution of the Kadanoff-Baym equations
     call kadanof_baym2Gloc(nstep)

     !use Dyson equation for  G0^r to update the Weiss field 
     call dyson4G0ret()

     !Update the Weiss fields and Sigma
     call update_G0gtr_less()
     call update_G0lceil()
     call update_sigma()
  enddo
  print*,"CIAO"
end PROGRAM neqIPT

