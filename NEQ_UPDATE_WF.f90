!###############################################################
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module NEQ_UPDATE_WF
  USE MATRIX
  USE NEQ_VARS_GLOBAL
  implicit none
  private

  public                           :: neq_update_weiss_field

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine neq_update_weiss_field
    integer :: i,j,k
    !These are used in the default inversion
    complex(8),dimension(2*Nstep,2*Nstep) :: mat_locG
    complex(8),dimension(2*Nstep,2*Nstep) :: mat_Sigma
    complex(8),dimension(2*Nstep,2*Nstep) :: mat_G0
    !These are used on call for case(2)
    complex(8),dimension(nstep,nstep) :: locGret,Sret
    complex(8),dimension(nstep,nstep) :: locGless,locGgtr
    complex(8),dimension(nstep,nstep) :: Sless,Sgtr
    complex(8),dimension(nstep,nstep) :: G0less,G0gtr
    complex(8),dimension(nstep,nstep) :: Uno,GammaRet,G0ret
    !These are used on call for case(3)
    complex(8),dimension(:,:),allocatable :: mat_Delta
    complex(8),dimension(:,:),allocatable :: mat_Gamma
    !
    type(keldysh_contour_gf),save             :: G0_old

    if(G0_old%status.EQV..false.)call allocate_keldysh_contour_gf(G0_old,Nstep*(Nstep+1)/2)
    G0_old=G0    

    call msg("Update WF")

    select case(fupdate)
    case default
       call msg("update with direct inversion of keldysh matrix GF")
       mat_locG = build_keldysh_matrix_gf(locG,nstep)   !Build Gloc matrix
       mat_Sigma = build_keldysh_matrix_gf(Sigma,nstep) !Build Sigma matrix
       mat_locG = mat_locG*dt**2                        !Invert Gloc
       call matrix_inverse(mat_locG)                    !
       mat_G0 = mat_locG + mat_Sigma                    !Get G0^-1 matrix:
       mat_G0 = mat_G0*dt**2                            !Invert G0^-1
       call matrix_inverse(mat_G0)
       do i=1,Nstep
          do j=1,Nstep
             k = pack_index_tri(i,j)
             G0%less(k) = -mat_G0(i,Nstep+j)
             G0%gtr(k)  =  mat_G0(Nstep+i,j)
          enddo
       enddo

    case(2)
       call msg("update with equations for <,>")
       do i=1,nstep
          do j=1,nstep
             locGless(i,j) = pack_less_tri(i,j,locG)
             locGgtr(i,j)  = pack_gtr_tri(i,j,locG)
             Sless(i,j)    = pack_less(i,j,Nstep,Sigma)
             Sgtr(i,j)     = pack_gtr(i,j,Nstep,Sigma)
             locGret(i,j)  = heaviside(t(i)-t(j))*(locGgtr(i,j)-locGless(i,j))
             Sret(i,j)     = heaviside(t(i)-t(j))*(Sgtr(i,j)-Sless(i,j))
          enddo
       enddo

       ! !G0ret = [\11 + Gret * Sret]^-1 * Gret
       Uno=zero  ; forall(i=1:nstep)Uno(i,i)=One/dt
       GammaRet = Uno+matmul(locGret,Sret)*dt
       GammaRet = GammaRet*dt**2
       call matrix_inverse(GammaRet)
       G0ret = matmul(GammaRet,locGret)*dt
       !### COMMENTING THIS LINE THE RESULTS ARE IDENTICAL WITH THE THREE METHODS OF UPDATE ###
       !forall(i=1:nstep)G0ret(i,i)=-xi !???
       !G0less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sless * gA
       G0less = matmul(GammaRet,matmul(locGless,conjg(transpose(GammaRet)))*dt)*dt -&
            matmul(G0ret,matmul(Sless,conjg(transpose(G0ret)))*dt)*dt
       !G0gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sgtr * gA
       G0gtr  = matmul(GammaRet,matmul(locGgtr,conjg(transpose(GammaRet)))*dt)*dt  -&
            matmul(G0ret,matmul(Sgtr,conjg(transpose(G0ret)))*dt)*dt
       do i=1,Nstep
          do j=1,Nstep
             k = pack_index_tri(i,j)
             G0%less(k) = G0less(i,j)
             G0%gtr(k)  = G0gtr(i,j)
          enddo
       enddo


    case(3)
       call msg("update with inversion and matrix-multiplication")
       !Build Gloc matrix
       mat_locG = build_keldysh_matrix_gf(locG,nstep)
       !Build Sigma matrix
       mat_Sigma= build_keldysh_matrix_gf(Sigma,nstep)
       !Allocate space for other matrices:
       allocate(mat_Delta(2*nstep,2*nstep))
       allocate(mat_Gamma(2*nstep,2*nstep))
       mat_Delta = zero ; forall(i=1:2*nstep)mat_Delta(i,i)=One/dt
       mat_Gamma = mat_Delta + matmul(mat_Sigma,mat_locG)*dt
       mat_Gamma = mat_Gamma*dt**2
       call matrix_inverse(mat_Gamma)
       mat_G0  = matmul(mat_locG,mat_Gamma)*dt
       do i=1,Nstep
          do j=1,Nstep             
             k = pack_index_tri(i,j)
             G0%less(k) = -mat_G0(i,Nstep+j)
             G0%gtr(k)  =  mat_G0(Nstep+i,j)
          enddo
       enddo
       deallocate(mat_Delta,mat_Gamma)

    end select


    G0%less = weight*G0%less + (1.d0-weight)*G0_old%less
    G0%gtr  = weight*G0%gtr  + (1.d0-weight)*G0_old%gtr

    !Save data:
    call write_keldysh_contour_gf(G0,"G0")
  end subroutine neq_update_weiss_field



  !********************************************************************
  !********************************************************************
  !********************************************************************



  function build_keldysh_matrix_gf(G,N) result(matG)
    type(keldysh_contour_gf)      :: G
    complex(8),dimension(2*N,2*N) :: matG
    integer                       :: i,j,k,N
    complex(8)                    :: Ggtr,Gless
    do i=1,N
       do j=1,N
          Gless = pack_less_tri(i,j,G)
          Ggtr  = pack_gtr_tri(i,j,G)
          matG(i,j)     = step(t(i)-t(j))*Ggtr + step(t(j)-t(i))*Gless
          matG(i,N+j)   =-Gless
          matG(N+i,j)   = Ggtr
          matG(N+i,N+j) =-( step(t(i)-t(j))*Gless + step(t(j)-t(i))*Ggtr )
       enddo
    enddo
  end function build_keldysh_matrix_gf



  !********************************************************************
  !********************************************************************
  !********************************************************************



end module NEQ_UPDATE_WF
