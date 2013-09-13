!###############################################################
!     PURPOSE  : A non-equilibrium IPT solver module. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module NEQ_IPT
  USE NEQ_VARS_GLOBAL
  implicit none
  private

  public  :: neq_solve_ipt



contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve with the 2^nd IPT sigma functions
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt()
    integer      :: i,j,itau
    real(8),dimension(nstep)            :: nt             !occupation(time)
    complex(8),dimension(Nstep,Nstep)   :: G0gtr,Sigmagtr

    call msg("Get Sigma(t,t')")
    G0gtr = G0%less + G0%ret - conjg(transpose(G0%ret))
    forall(i=1:nstep,j=1:nstep)
       Sigma%less(i,j) = (U**2)*(G0%less(i,j)**2)*G0gtr(j,i)    !get Sigma^<(t,t`)
       Sigmagtr(i,j)   = (U**2)*(G0gtr(i,j)**2)*G0%less(j,i)    !get Sigma^>(t,t`)
    end forall
    forall(i=1:nstep,j=1:nstep)Sigma%ret(i,j)  = heaviside(time(i)-time(j))*(Sigmagtr(i,j)-Sigma%less(i,j))

    !Save data:
    call write_keldysh_contour_gf(Sigma,reg(data_dir)//"/Sigma")
    if(plot3D)call plot_keldysh_contour_gf(Sigma,time,reg(plot_dir)//"/Sigma")
    forall(i=1:nstep)nt(i)=-xi*Sigma%less(i,i)
    call splot("nsVStime.ipt",time,nt,append=.true.)
  end subroutine neq_solve_ipt



end module NEQ_IPT
