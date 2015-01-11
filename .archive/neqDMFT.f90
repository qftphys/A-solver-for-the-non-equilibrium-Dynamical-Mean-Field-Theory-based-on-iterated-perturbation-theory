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
  USE FUNX_NEQ                  !contains all the neqDMFT functions
  USE KADANOFBAYM               !solves KB equations for k-summation
  implicit none
  character(len=4) :: char

  !INIT THE CALCULATION: Read inputs & setup the grids:
  call read_input_init()

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

  !Creates a directory where to store data:
  call system("if [ ! -d DATAsrc ]; then mkdir DATAsrc; fi")

  !Get the Guess // Build dissipative Bath :
  call guess_G0
  call get_Bath

  !Start DMFT iteration
  do iloop=1,nloop
     call dump("",lines=3)
     write(*,"(A,i5,A,i5)")"DMFT-loop",iloop,"/",nloop
     call dump("---------------------------------------")
     if(iloop>1)call update_G0_Dyson(update_wfftw)
     call get_Sigma(method)
     call get_Gloc(iloop,solve_wfftw)
     if(trim(method)=="ipt")call obtain_Gimp
     call print_observables

     !Move the actual result to DIR_iloop
     write(char,'(i4)')iloop
     call system("if [ -d DATAsrc_iter"//adjustl(trim(char))//" ]; then rm -rf DATAsrc_iter"//adjustl(trim(char))//"; fi")
     call system("mv -v DATAsrc DATAsrc_iter"//adjustl(trim(char)))
     !If you succeed doing that you can remove the previous loop dir: DIR_(iloop-1)
     write(char,'(i4)')iloop-1
     call system("if [ -d DATAsrc_iter"//adjustl(trim(char))//" ]; then rm -rf DATAsrc_iter"//adjustl(trim(char))//"; fi")
  enddo
  !Rename the remaining data DIR
  write(char,'(i4)')nloop
  call system("mv -v DATAsrc_iter"//adjustl(trim(char))//" "//adjustl(trim(SRC)))  
  call system("tar czvf "//adjustl(trim(SRC))//".tgz "//adjustl(trim(SRC)))
  call system("rm -rf "//adjustl(trim(SRC)))
  call alloc_memory('d')
  print*,"BRAVO!"
end PROGRAM neqDMFT
