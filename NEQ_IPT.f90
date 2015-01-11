!###############################################################
!     PURPOSE  : A non-equilibrium IPT solver module. 
!###############################################################
module NEQ_IPT
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_VARS_GLOBAL
  USE NEQ_INPUT_VARS
  USE CONSTANTS
  USE DMFT_TOOLS
  implicit none
  private

  public  :: neq_solve_ipt_first_step
  public  :: neq_solve_ipt

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve with the 2^nd IPT sigma functions
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt_first_step(g0,self,params)
    type(kb_contour_gf)     :: g0
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    complex(8)              :: Self_gtr
    integer                 :: i,j,k,ik,unit,len,N,L,Lf
    real(8),dimension(0:1)  :: Scoeff
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.self%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    N = params%Itime
    L = params%Ntau
    Lf= params%Niw
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    !(this step depends on the imp. solv.)
    ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
    ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
    ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
    ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
    ! self^R(0,0) = self^> - self^<
    do j=1,L
       Self%mats(j) = U0*U0*g0%mats(j)*g0%mats(L-j+1)*g0%mats(j)
    end do
    Scoeff  = tail_coeff_sigma(uloc,0.5d0)
    call fft_sigma_tau2iw(Self%iw,Self%mats,beta,Scoeff)
    Self%iw = xi*dimag(self%iw) !!ACTHUNG: imposing half-filling symmetry
    Self%less(1,1)=(xi**3)*Uloc*Uloc*g0%mats(L)*g0%mats(1)*g0%mats(L)
    Self_gtr      =(xi**3)*Uloc*Uloc*g0%mats(1)*g0%mats(L)*g0%mats(1)
    do j=1,L
       Self%lmix(1,j)=(xi**3)*Uloc*U0*g0%mats(L-j+1)*g0%mats(j)*g0%mats(L-j+1)
    end do
    Self%ret(1,1) = Self_gtr - Self%less(1,1)
  end subroutine neq_solve_ipt_first_step




  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve with the 2^nd IPT sigma functions
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt(G0,self,params)
    type(kb_contour_gf)                   :: G0
    type(kb_contour_gf)                   :: self
    type(kb_contour_params)               :: params
    integer                               :: N,L
    complex(8),dimension(:,:),allocatable :: G0_gtr,self_gtr,G0_rmix
    integer                               :: i,j,itau
    !
    N   = params%Itime                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    allocate(G0_gtr(N,N),self_gtr(N,N),G0_rmix(L,N))
    do j=1,N
       G0_gtr(N,j)=G0%less(N,j)+G0%ret(N,j)
    end do
    do i=1,N-1
       G0_gtr(i,N)=G0%less(i,n)-conjg(G0%ret(N,i))
    end do
    do j=1,L
       G0_rmix(j,N)  = conjg(G0%lmix(N,L-j+1))
    enddo

    !Vertical edge
    do j=1,N
       self%less(N,j)= Uloc*Uloc*G0%less(N,j)*G0_gtr(j,N)*G0%less(N,j)
       self_gtr(N,j) = Uloc*Uloc*G0_gtr(N,j)*G0%less(j,N)*G0_gtr(N,j)
    end do
    !Horizontal edge
    do i=1,N-1
       self%less(i,N)= Uloc*Uloc*G0%less(i,N)*G0_gtr(N,i)*G0%less(i,N)
       self_gtr(i,N) = Uloc*Uloc*G0_gtr(i,N)*G0%less(N,i)*G0_gtr(i,N)
    end do
    !Imaginary time edge:
    forall(i=1:L)self%lmix(N,i)  = Uloc*U0*G0%lmix(N,i)*G0_rmix(i,N)*G0%lmix(N,i)
    forall(j=1:N)self%ret(N,j) = self_gtr(N,j) - self%less(N,j)


    ! !#################################
    ! if(ifourth)then
    !    print*,"Entering 4th order..."

    !    Ntot = 2*N+L+1

    !    do i=1,4
    !       call allocate_kb_contour_gf(sigma4(i),params)
    !       sigma4(i) = zero
    !    enddo

    !    allocate(G0mat(Ntot,Ntot),G2mat(Ntot,Ntot),G3mat(Ntot,Ntot))
    !    allocate(Chimat(Ntot,Ntot),Smat(Ntot,Ntot),Sfoo(Ntot,Ntot))
    !    allocate(tdiff(Ntot),Utime(Ntot))

    !    tdiff(1:N)           =  params%dt
    !    tdiff(N+1:2*N)       = -params%dt
    !    tdiff(2*N+1:2*N+L+1) = -xi*params%dtau

    !    Utime(1:N)           = U 
    !    Utime(N+1:2*N)       = U
    !    Utime(2*N+1:2*N+L+1) = Ui


    !    !Sigma^(4a):
    !    !================================================================
    !    print*,"get 4a..."
    !    call kb_contour_gf2kb_matrix(G0,N,L,G0mat)
    !    do i=1,Ntot
    !       do j=1,Ntot
    !          G2mat(i,j) = G0mat(i,j)*G0mat(i,j)
    !          G3mat(i,j) = G0mat(i,j)*G0mat(i,j)*G0mat(i,j)
    !       enddo
    !    enddo
    !    print*,"convolute [G0^2*G0^2]:"
    !    call convolute_kb_matrix_gf(G2mat,G2mat,N,L,params,Chimat)
    !    print*,"convolute [G0^2*Chi]:"
    !    call convolute_kb_matrix_gf(G2mat,Chimat,N,L,params,Sfoo)
    !    do i=1,Ntot
    !       do j=1,Ntot
    !          Smat(i,j) = G0mat(i,j)*Sfoo(i,j)
    !       enddo
    !    enddo
    !    print*,"extract Sigma^a:"
    !    call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(1))





    !    !Sigma^(4b):
    !    print*,"Get 4b..."
    !    print*,"convolute [G0^3*G0]:"
    !    call convolute_kb_matrix_gf(G3mat,G0mat,N,L,params,Chimat)
    !    print*,"convolute [G0*Chi]:"
    !    call convolute_kb_matrix_gf(G0mat,Chimat,N,L,params,Sfoo)
    !    do i=1,Ntot
    !       do j=1,Ntot
    !          Smat(i,j) = G2mat(i,j)*Sfoo(i,j)
    !       enddo
    !    enddo
    !    print*,"extract Sigma^b:"
    !    call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(2))



    !    !Sigma^(4c):
    !    print*,"Get 4c"
    !    Smat  = zero
    !    do i=1,Ntot
    !       do j=1,Ntot
    !          int_ij=zero
    !          do ia=1,Ntot
    !             intb=zero
    !             do ib=1,Ntot
    !                intb = intb + G0mat(i,ia)*G0mat(ia,ib)*G0mat(ib,j)*G2mat(i,ib)*G2mat(ia,j)*tdiff(ib)
    !             enddo
    !             int_ij = int_ij + intb*tdiff(ia)
    !          enddo

    !          Smat(i,j) = int_ij
    !       enddo
    !    enddo
    !    print*,"extract Sigma^c:"
    !    call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(3))





    !    !Sigma^(4d):
    !    print*,"Get 4d"
    !    Smat  = zero
    !    do i=1,Ntot
    !       do j=1,Ntot
    !          int_ij=zero
    !          do ia=1,Ntot
    !             intb=zero
    !             do ib=1,Ntot
    !                intb = intb + G0mat(i,ia)*G0mat(ia,j)*G0mat(i,ib)*G0mat(ib,j)*G2mat(ia,ib)*tdiff(ib)
    !             enddo
    !             int_ij = int_ij + intb*tdiff(ia)
    !          enddo
    !          Smat(i,j) = G0mat(i,j)*int_ij
    !       enddo
    !    enddo
    !    print*,"extract Sigma^d:"
    !    call kb_matrix2kb_contour_gf(Smat,N,L,Sigma4(4))


    !    !call plot_kb_contour_gf("Sigma4a_test",Sigma4(1),params)
    !    !call plot_kb_contour_gf("Sigma4b_test",Sigma4(2),params)
    !    !call plot_kb_contour_gf("Sigma4c",Sigma4(3),params)
    !    !call plot_kb_contour_gf("Sigma4d",Sigma4(4),params)

    !    !============================================





    !    print*,"Sum up Sigma2+Sigma4:"
    !    do j=1,N
    !       Sigma%less(N,j) = Sigma%less(N,j) + &
    !            3.d0*U*U*U*U*(Sigma4(1)%less(N,j) + Sigma4(2)%less(N,j) + Sigma4(3)%less(N,j) - Sigma4(4)%less(N,j))
    !       Sigma%ret(N,j) = Sigma%ret(N,j) + &
    !            3.d0*U*U*U*U*(Sigma4(1)%ret(N,j) + Sigma4(2)%ret(N,j) + Sigma4(3)%ret(N,j) - Sigma4(4)%ret(N,j))
    !    enddo
    !    do i=1,N-1
    !       Sigma%less(i,N) = Sigma%less(i,N) + &
    !            3.d0*U*U*U*U*(Sigma4(1)%less(i,N) + Sigma4(2)%less(i,N) + Sigma4(3)%less(i,N) - Sigma4(4)%less(i,N))          
    !    enddo
    !    !!ACTHUNG!! U_i might not be correct in this expression!!
    !    do j=0,L
    !       Sigma%lmix(N,j) = Sigma%lmix(N,j) + &
    !            3.d0*U*U*U*Ui*(Sigma4(1)%lmix(N,j) + Sigma4(2)%lmix(N,j) + Sigma4(3)%lmix(N,j) - Sigma4(4)%lmix(N,j))
    !    enddo


    !    print*,"Deallocating:"
    !    do i=1,4
    !       call deallocate_kb_contour_gf(sigma4(i))
    !    enddo
    !    deallocate(G0mat,G2mat,G3mat)
    !    deallocate(Chimat,Smat,Sfoo)
    !    deallocate(tdiff)
    ! endif

  end subroutine neq_solve_ipt













end module NEQ_IPT
