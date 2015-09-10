module NEQ_IPT
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_INPUT_VARS
  USE NEQ_MEASURE
  implicit none
  private

  interface neq_solve_ipt
     module procedure neq_solve_ipt_normal
  end interface neq_solve_ipt

  public  :: neq_solve_ipt

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: Solve with the 2^nd IPT sigma functions
  !  - _normal: solve the 2nd order IPT for the normal (repulsive case)
  ! at half-filling.
  !TODO 1: work out the case with interaction in the form Unn.
  !TODO 2: extend the routines to the non-half-filling.
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt_normal(G0,G,Sigma,params)
    type(kb_contour_gf)                   :: G0
    type(kb_contour_gf)                   :: G
    type(kb_contour_sigma)                :: Sigma
    type(kb_contour_params)               :: params
    integer                               :: N,L
    complex(8),dimension(:,:),allocatable :: G0_gtr,G0_rmix,Sigma_reg_gtr
    real(8)                               :: dens
    integer                               :: i,j,itau
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    !<DEBUG add test on G0,G,Sigma
    !
    !>DEBUG
    !
    !Get the HF contribution:
    ! in this formulation this quantity is zero
    ! because is [or it should be] n=0.5, while Sigma_HF=U*(n-0.5).
    ! In the future we may want to implement something where the HF 
    ! contribution is actually evaluated, either because n!=0.5 or
    ! interaction is such that this term is of the form Sigma_HF=U*n
    dens = measure_dens(g,params)
    Sigma%hf(N) = zero          !U*(dens-0.5d0)
    !
    !Get the Regular contribution:
    allocate( G0_gtr(N,N), G0_rmix(0:L,N), Sigma_reg_gtr(N,N) )
    do j=1,N
       G0_gtr(N,j)=G0%less(N,j)+G0%ret(N,j)
    end do
    do i=1,N-1
       G0_gtr(i,N)=G0%less(i,n)-conjg(G0%ret(N,i))
    end do
    do j=0,L
       G0_rmix(j,N)  = conjg(G0%lmix(N,L-j))
    enddo
    !
    !Vertical edge
    do j=1,N
       Sigma%reg%less(N,j)= U*U*G0%less(N,j)*G0_gtr(j,N)*G0%less(N,j)
       Sigma_reg_gtr(N,j) = U*U*G0_gtr(N,j)*G0%less(j,N)*G0_gtr(N,j)
    end do
    !Horizontal edge
    do i=1,N-1
       Sigma%reg%less(i,N)= U*U*G0%less(i,N)*G0_gtr(N,i)*G0%less(i,N)
       Sigma_reg_gtr(i,N) = U*U*G0_gtr(i,N)*G0%less(N,i)*G0_gtr(i,N)
    end do
    !Imaginary time edge:
    forall(i=0:L)Sigma%reg%lmix(N,i) = U*Ui*G0%lmix(N,i)*G0_rmix(i,N)*G0%lmix(N,i)
    forall(j=1:N)Sigma%reg%ret(N,j)    = Sigma_reg_gtr(N,j) - Sigma%reg%less(N,j)

  end subroutine neq_solve_ipt_normal












end module NEQ_IPT
