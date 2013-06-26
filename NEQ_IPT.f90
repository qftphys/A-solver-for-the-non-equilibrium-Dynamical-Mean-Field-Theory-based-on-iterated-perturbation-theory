!###############################################################
!     PURPOSE  : A non-equilibrium IPT solver module. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module NEQ_IPT
  USE NEQ_VARS_GLOBAL
  implicit none
  private

  public  :: neq_init_run
  public  :: neq_solve_ipt



contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading self-energy
  !+-------------------------------------------------------------------+
  subroutine neq_init_run()
    logical :: init,checknk
    integer :: i,j,ik

    init = inquire_keldysh_contour_gf(irdSFILE)
    if(init)then
       call msg(bold("Reading components of the input Self-energy and nk"))
       call read_keldysh_contour_gf(Sigma,irdSFILE)
    else
       call msg(bold("Start from the Hartree-Fock self-energy"))
       Sigma=zero
    endif

    inquire(file=reg(irdnkfile),exist=checkNk)
    if(.not.checkNk)inquire(file=reg(irdNkfile)//".gz",exist=checkNk)
    if(checkNk)then
       call read_nkfile(eq_nk,trim(irdnkfile))
    else
       !Get non-interacting n(k):
       do ik=1,Lk
          eq_nk(ik)=fermi((epsik(ik)),beta)
       enddo
    endif
    call splot("init_nkVSepsk.ipt",epsik,eq_nk)

  contains

    subroutine read_nkfile(irdnk,file)
      character(len=*)     :: file
      real(8),dimension(Lk):: irdnk
      integer              :: redLk
      real(8),allocatable  :: rednk(:),redek(:)
      integer,allocatable  :: orderk(:)
      real(8),allocatable  :: uniq_rednk(:),uniq_redek(:)
      logical,allocatable  :: maskk(:)
      !n(k): A lot of work here to reshape the array
      redLk=file_length(file)
      allocate(rednk(redLk),redek(redLk),orderk(redLk))
      call sread(file,redek,rednk)
      !work on the read arrays:
      !1 - sorting: sort the energies (X-axis), mirror on occupation (Y-axis) 
      !2 - delete duplicates energies (X-axis), mirror on occupation (Y-axis) 
      !3 - interpolate to the actual lattice structure (epsik,nk)
      call sort_array(redek,orderk)
      call reshuffle(rednk,orderk)
      call uniq(redek,uniq_redek,maskk)
      allocate(uniq_rednk(size(uniq_redek)))
      uniq_rednk = pack(rednk,maskk)
      call linear_spline(uniq_rednk,uniq_redek,irdnk,epsik)
    end subroutine read_nkfile
  end subroutine neq_init_run




  !+-------------------------------------------------------------------+
  !PURPOSE  : BUild the 2^nd IPT sigma functions:
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt()
    integer    :: i,j,k
    complex(8) :: G0less,G0gtr
    call msg("Get Sigma(t,t')")
    do i=1,nstep
       do j=1,nstep
          k =pack_index(i,j,nstep)
          G0less = pack_less_tri(i,j,G0)
          G0gtr  = pack_gtr_tri(j,i,G0)
          Sigma%less(k) = U**2*G0less**2*G0gtr
          !
          G0gtr = pack_gtr_tri(i,j,G0)
          G0less= pack_less_tri(j,i,G0)
          Sigma%gtr(k)  = U**2*G0gtr**2*G0less
       enddo
    enddo
    !Save data:
    call msg("Saving Sigma function:")
    call write_keldysh_contour_gf(Sigma,"Sigma")
  end subroutine neq_solve_ipt



  !********************************************************************
  !********************************************************************
  !********************************************************************






end module NEQ_IPT
