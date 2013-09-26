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
    integer      :: i,j,itau,ints,intg
    real(8),dimension(nstep)            :: nt             !occupation(time)
    complex(8),dimension(Nstep,Nstep)   :: G0gtr,Sigmagtr

    call msg("Get Sigma(t,t')")
    G0gtr = zero
    forall(i=1:Nstep,j=1:Nstep,i>=j)G0gtr(i,j) = G0%less(i,j) + G0%ret(i,j)
    forall(i=1:Nstep,j=1:Nstep,i<j)G0gtr(i,j)=-conjg(G0gtr(j,i))


    Sigma   = zero
    Sigmagtr= zero
    forall(i=1:nstep,j=1:nstep,i>=j)
       Sigma%less(i,j) = U*U*G0%less(i,j)*G0%less(i,j)*G0gtr(j,i) !get Sigma^<(t,t`)
       Sigmagtr(i,j)   = U*U*G0gtr(i,j)*G0gtr(i,j)*G0%less(j,i)   !get Sigma^>(t,t`)
    end forall                                                    !
    Sigma%ret = Sigmagtr - Sigma%less                             ! i>=j
    forall(i=1:Nstep,j=1:Nstep,i<j)                               !
       Sigma%less(i,j)=-conjg(Sigma%less(j,i))
       Sigmagtr(i,j)=-conjg(Sigmagtr(j,i))
    end forall

    !Save data:
    call write_keldysh_contour_gf(Sigma,reg(data_dir)//"/Sigma")
    if(plot3D)call plot_keldysh_contour_gf(Sigma,time,reg(plot_dir)//"/Sigma")

    !<<<DEBUG
    forall(i=1:nstep)nt(i)=-xi*Sigma%less(i,i)
    call splot("nsVStime.ipt",time,2d0*nt,append=.true.)
    forall(i=1:nstep)nt(i)=-xi*G0%less(i,i)
    call splot("nt0VStime.ipt",time,2d0*nt,append=.true.)
    ints=500;intg=100
    do j=1,Nstep,Nstep/10
       ints=ints+1
       intg=intg+1
       rewind(intg);rewind(ints)
       do i=1,Nstep
          write(ints,"(7F26.16)")time(i),dimag(Sigma%less(i,j)),dreal(Sigma%less(i,j)),&
               dimag(Sigmagtr(i,j)),dreal(Sigmagtr(i,j)),&
               dimag(Sigma%ret(i,j)),dreal(Sigma%ret(i,j))
          write(intg,"(7F26.16)")time(i),dimag(G0%less(i,j)),dreal(G0%less(i,j)),&
               dimag(G0gtr(i,j)),dreal(G0gtr(i,j)),&
               dimag(G0%ret(i,j)),dreal(G0%ret(i,j))
       enddo
    enddo
    !>>>DEBUG

  end subroutine neq_solve_ipt



end module NEQ_IPT
