program getDATA
  USE VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE BATH
  USE RESULTS
  implicit none

  call read_input_init("used.inputFILE.in")
  include "grid_setup.f90"  
  Lk   = square_lattice_dimension(Nx,Ny)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx,Ny)
  epsik= square_lattice_dispersion_array(Lk,ts=ts)

  call set_efield_vector()

  call get_data()

  call plot_results()

contains

  !+-------------------------------------------------------------------+
  subroutine get_data()
    logical :: control

    call global_memory_allocation()

    call get_thermostat_bath()

    !Read the functions:
    call read_keldysh_contour_gf(G0,"G0")
    call read_keldysh_contour_gf(locG,"locG")
    call read_keldysh_contour_gf(Sigma,"Sigma")

    call sread("nk.neq",nk(:,1:Lk))

    if(fchi)then
       call sread("locChi_11.neq",chi(1,1,:,:))
       call sread("locChi_12.neq",chi(1,2,:,:))
       call sread("locChi_21.neq",chi(2,1,:,:))
       call sread("locChi_22.neq",chi(2,2,:,:))
    endif

    return
  end subroutine get_data
  !+-------------------------------------------------------------------+

end program getDATA
