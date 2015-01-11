Program neqDMFT
  !###################################################################
  !PROGRAM  : freeNEQ
  !TYPE     : Main program
  !PURPOSE  : Solve system of free/interacting electrons driven 
  ! by electric field interacting with a resevoir/bath of free 
  ! electrons at temperature T
  !AUTHORS  : Adriano Amaricci && Cedric Weber
  !###################################################################
  !LIBRARY:  
  USE TOOLS
  USE LATTICE
  !LOCAL:
  USE VARS_GLOBAL 
  USE EQKELDYSH   
  USE OBSERVABLES   
  USE FUNX_NEQ    
  USE KADANOFBAYM 
  implicit none

  !INIT THE CALCULATION: Read inputs & setup the grids:
  call read_input("inputFILE.in")
  call init_calc()

  !===========================================!
  !BUILD LATTICE STRUCTURE:
  call dump("Get Lattice Structure:")  
  call Build_2DSquareLattice(Nx,Ny,Lk)
  allocate(epsik(Lk),epsimu(Lmu),sorted_epsik(Lk),sorted_ik(Lk))
  call get_epsik(epsik,ts,0.d0,sorted_epsik,sorted_ik)
  call get_epsimu(emin,epsimu)
  !===========================================!

  !SET EXTERNAL ELECTRIC FIELD:
  Ek%x=1.d0;Ek%y=1.d0;Ek=(Efield*pi2/alat)*Ek !

  !STARTS THE REAL WORK:
  !===========================================!
  !Massive allocation of the working arrays:
  call dump("")
  call alloc_memory('a')

  !Solve the associated equilibrium problem if needed:
#ifdef _mix
  eqflag=.true.
#endif
  if(eqflag)call keldysheq(eqnloop)

  !Get the Guess // Build dissipative Bath // Get 1st Sigma:
  call guess_G0
  call get_Bath
  call get_Sigma(method)

  !Starts the DMFT-loop
  do iloop=1,nloop
     call dump("",lines=3)
     write(*,"(A,i5,A,i5)")"DMFT-loop",iloop,"/",nloop
     call dump("---------------------------------------")

     call get_Gloc_KadanoffBaym(iloop)
#ifdef _mix
     stop
#endif
     call update_G0_Dyson()
#ifdef _mix
     stop
#endif
     !-----------------------
     !this can (in principle) substituted with any Impurity Solver
     call get_Sigma(method)
     !----------------------

     call evaluate_print_observables(iloop)
     !call save_solution()
  enddo
  !===========================================!




  !Massive deallocation of the working arrays:
  !===========================================!
  call alloc_memory('d')
  print*,"BRAVO!"
  !===========================================!
end PROGRAM neqDMFT
!********************************************************************
!********************************************************************
!********************************************************************
