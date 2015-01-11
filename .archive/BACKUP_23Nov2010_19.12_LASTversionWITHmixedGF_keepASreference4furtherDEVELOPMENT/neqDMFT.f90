!###################################################################
!PROGRAM  : freeNEQ
!TYPE     : Main program
!PURPOSE  : Solve system of free/interacting electrons driven 
! by electric field interacting with a resevoir/bath of free 
! electrons at temperature T
!AUTHORS  : Adriano Amaricci && Cedric Weber
!###################################################################
program neqDMFT
  !LOCAL:
  USE VARS_GLOBAL               !global variables, calls to 3rd library 
  USE EQKELDYSH                 !solves the associated equilibrium pbm.
  USE OBSERVABLES               !evaluates and prints system observables
  USE FUNX_NEQ                  !contains all the neqDMFT functions
  USE KADANOFBAYM               !solves KB equations for k-summation
  implicit none

  !INIT THE CALCULATION: Read inputs & setup the grids:
  call read_input_init("inputFILE.in")

  !BUILD LATTICE STRUCTURE:
  call dump("Get Lattice Structure:")  
  call Build_2DSquareLattice(Nx,Ny,Lk)
  allocate(epsik(Lk),epsimu(Lmu),sorted_epsik(Lk),sorted_ik(Lk))
  call get_epsik(epsik,ts,0.d0,sorted_epsik,sorted_ik)
  call get_epsimu(emin,epsimu)

  !SET EXTERNAL ELECTRIC FIELD:
  Ek%x=1.d0;Ek%y=1.d0;Ek=(Efield*pi2/alat)*Ek !

  !STARTS THE REAL WORK:
  !Massive allocation of the working arrays:
  call dump("")
  call alloc_memory('a')

  !Solve the associated equilibrium problem if needed:
  call equilibium_ipt_keldysh(eqnloop)

  !Get the Guess // Build dissipative Bath :
  call guess_G0
  call get_Bath

  !Starts the DMFT-loop
  do iloop=1,nloop
     call dump("",lines=3)
     write(*,"(A,i5,A,i5)")"DMFT-loop",iloop,"/",nloop
     call dump("---------------------------------------")
     if(iloop>1)call update_G0_Dyson()
     call get_Sigma(method) 
     call get_Gloc_KadanoffBaym(iloop)
     !call get_Gloc_equilibrium
     call evaluate_print_observables(iloop)
  enddo

  !Massive deallocation of the working arrays:
  call alloc_memory('d')
  print*,"BRAVO!"
end PROGRAM neqDMFT
