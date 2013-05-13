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
    call read_keldysh_contour_gf(G0,trim(data_dir)//"/G0")
    call read_keldysh_contour_gf(locG,trim(data_dir)//"/locG")
    call read_keldysh_contour_gf(Sigma,trim(data_dir)//"/Sigma")

    call sread(trim(data_dir)//"/nk.data",nk(0:nstep,1:Lk))

    if(fchi)then
       call sread(trim(data_dir)//"/locChi_11.data",chi(1,1,0:nstep,0:nstep))
       call sread(trim(data_dir)//"/locChi_12.data",chi(1,2,0:nstep,0:nstep))
       call sread(trim(data_dir)//"/locChi_21.data",chi(2,1,0:nstep,0:nstep))
       call sread(trim(data_dir)//"/locChi_22.data",chi(2,2,0:nstep,0:nstep))
    endif

    return
  end subroutine get_data
  !+-------------------------------------------------------------------+

end program getDATA
