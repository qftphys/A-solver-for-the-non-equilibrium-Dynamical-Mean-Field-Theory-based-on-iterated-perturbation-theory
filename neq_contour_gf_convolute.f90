!C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max ; t_max_index==N
subroutine convolute_kb_contour_gf_gf(A,B,C,params,dcoeff,ccoeff)
  type(kb_contour_gf)                 :: A,B,C
  type(kb_contour_params)             :: params
  real(8),optional                    :: dcoeff
  complex(8),optional                 :: ccoeff
  integer                             :: N,L
  real(8)                             :: dt,dtau
  complex(8),dimension(:),allocatable :: AxB    
  integer                             :: i,j,k,itau,jtau
  logical                             :: checkA,checkB,checkC
  !
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  if(N==1) stop "convolute_kb_contour_gf_gf error: called with N=1"
  if(  (.not.A%status).OR.&
       (.not.B%status).OR.&
       (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: A,B,C not allocated"
  checkA=check_dimension_kb_contour(A,params%Ntime,L) 
  checkB=check_dimension_kb_contour(B,params%Ntime,L)
  checkC=check_dimension_kb_contour(C,params%Ntime,L)
  !
  allocate(AxB(0:max(L,N)))
  !
  !Ret. component
  !C^R(t,t')=\int_{t'}^t ds A^R(t,s)*B^R(s,t')
  C%ret(N,1:N)=zero
  do j=1,N           !for all t' in 0:t_max
     AxB(0:)  = zero !AxB should be set to zero before integration
     do k=j,N        !store the convolution between t'{=j} and t{=N}
        AxB(k) = A%ret(N,k)*B%ret(k,j)
     enddo
     C%ret(n,j) = C%ret(n,j) + dt*kb_trapz(AxB(0:),j,N)
  enddo
  !
  !Lmix. component
  !C^\lmix(t,tau')=\int_0^{beta} ds A^\lmix(t,s)*B^M(s,tau')
  !                  +\int_0^{t} ds A^R(t,s)*B^\lmix(s,tau')
  !               = I1 + I2
  C%lmix(N,0:L)=zero
  do jtau=0,L
     !I1:
     AxB(0:) = zero
     do k=0,jtau
        AxB(k)=A%lmix(N,k)*B%mats(L+k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
     do k=jtau,L
        AxB(k)=A%lmix(N,k)*B%mats(k-jtau)
     end do
     C%lmix(n,jtau)=C%lmix(n,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
     !I2:
     AxB(0:) = zero
     do k=1,N
        AxB(k) = A%ret(N,k)*B%lmix(k,jtau)
     enddo
     C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
  enddo
  !
  !Less component
  !C^<(t,t')=-i\int_0^{beta} ds A^\lmix(t,s)*B^\rmix(s,t')
  !             +\int_0^{t'} ds A^<(t,s)*B^A(s,t')
  !             +\int_0^{t} ds A^R(t,s)*B^<(s,t')
  ! (t,t')=>(N,j) <==> Vertical side, with tip (j=1,N) !!no tip (j=1,N-1)
  do j=1,N                      !-1
     C%less(N,j)=zero
     do k=0,L
        AxB(k)=A%lmix(N,k)*get_rmix(B,k,j,L)
     end do
     C%less(N,j)=C%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
     do k=1,j
        AxB(k)=A%less(N,k)*get_adv(B,k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
     do k=1,N
        AxB(k)=A%ret(N,k)*B%less(k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
  end do
  !
  ! (t,t`)=>(i,N) <==> Horizontal side, w/ no tip (i=1,N-1)
  do i=1,N-1
     C%less(i,N)=zero
  enddo
  ! do i=1,N-1
  !    C%less(i,N)=zero
  !    do k=0,L
  !       !         AxB(k)=A%lmix(i,k)*conjg(B%lmix(n,L-k))
  !       AxB(k)=A%lmix(i,k)*get_rmix(B,k,N,L)
  !    end do
  !    C%less(i,N)=C%less(i,N)-xi*dtau*kb_trapz(AxB(0:),0,L)
  !    !
  !    do k=1,N
  !       !         AxB(k)=A%less(i,k)*conjg(B%ret(N,k))
  !       AxB(k)=A%less(i,k)*get_adv(B,k,N)
  !    end do
  !    C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,N)
  !    !
  !    do k=1,i
  !       AxB(k)=A%ret(i,k)*B%less(k,N)
  !    end do
  !    C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,i)
  ! end do
  !    
  if(present(dcoeff))then
     C%lmix(N,0:L) = dcoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=dcoeff*C%less(N,1:N-1)
     C%less(1:N,N)=dcoeff*C%less(1:N,N)
     C%ret(N,1:N)=dcoeff*C%ret(N,1:N)
  endif
  if(present(ccoeff))then
     C%lmix(N,0:L) = ccoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=ccoeff*C%less(N,1:N-1)
     C%less(1:N,N)=ccoeff*C%less(1:N,N)
     C%ret(N,1:N)=ccoeff*C%ret(N,1:N)
  endif
  deallocate(AxB)
end subroutine convolute_kb_contour_gf_gf


!C(t,t')=[A*(Bhf*delta + Breg)](t,t'), with t=t_max && t'=0,t_max ; t_max_index==N
subroutine convolute_kb_contour_gf_sigma(A,B,C,params,dcoeff,ccoeff)
  type(kb_contour_gf)                 :: A
  type(kb_contour_sigma)              :: B
  type(kb_contour_gf)                 :: C
  type(kb_contour_params)             :: params
  real(8),optional                    :: dcoeff
  complex(8),optional                 :: ccoeff
  integer                             :: N,L
  real(8)                             :: dt,dtau
  complex(8),dimension(:),allocatable :: AxB    
  integer                             :: i,j,k,itau,jtau
  logical                             :: checkA,checkB,checkC
  !
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  if(N==1) stop "convolute_kb_contour_gf_sigma error: called with N=1"
  if(  (.not.A%status).OR.&
       (.not.B%status).OR.&
       (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: A,B,C not allocated"
  checkA=check_dimension_kb_contour(A,params) 
  checkB=check_dimension_kb_contour(B,params)
  checkC=check_dimension_kb_contour(C,params)
  !
  allocate(AxB(0:max(L,N)))
  !
  !Ret. component
  !C^R(t,t')=A^R(t,t')*Bhf(t') + \int_{t'}^t ds A^R(t,s)*Breg^R(s,t')
  C%ret(N,1:N)=zero
  do j=1,N
     C%ret(N,j) = A%ret(N,j)*B%hf(j)
     AxB(0:)=zero
     do k=j,N
        AxB(k) = A%ret(N,k)*B%reg%ret(k,j)
     enddo
     C%ret(N,j) = C%ret(N,j) + dt*kb_trapz(AxB(0:),j,N)
  enddo
  !
  !Lmix. component
  !C^\lmix(t,tau')= A^\lmix(t,tau')*Bhf(0+) +
  !                 \int_0^{beta} ds A^\lmix(t,s)*Breg^M(s,tau') + 
  !                 \int_0^{t} ds A^R(t,s)*Breg^\lmix(s,tau') 
  !               = A^\lmix*Bhf + I1 + I2
  C%lmix(N,0:L)=zero
  do jtau=0,L
     !A^\lmix*Bhf:
     C%lmix(N,jtau) = A%lmix(N,jtau)*B%hf(1)
     !I1:
     AxB(0:)=zero
     do k=0,jtau
        AxB(k)=A%lmix(N,k)*B%reg%mats(L+k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
     AxB(0:)=zero
     do k=jtau,L
        AxB(k)=A%lmix(N,k)*B%reg%mats(k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
     !I2:
     AxB(0:)=zero
     do k=1,N
        AxB(k) = A%ret(N,k)*B%reg%lmix(k,jtau)
     enddo
     C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
  enddo
  !
  !Less component
  !C^<(t,t') =  A^<(t,t')*conjg(Bhf(t'))
  !            -i\int_0^{beta} ds A^\lmix(t,s)*Breg^\rmix(s,t')
  !            + \int_0^{t'} ds A^<(t,s)*Breg^A(s,t')
  !            + \int_0^{t} ds A^R(t,s)*Breg^<(s,t')
  !          =  A^<*Bhf* + I1 + I2 + I3
  ! (t,t')=>(N,j) <==> Vertical side, with tip (j=1,N)
  C%less(N,1:N)=zero
  do j=1,N
     !A^<*Bhf*:
     C%less(N,j)=A%less(N,j)*conjg(B%hf(j))
     !I1:
     AxB(0:)=zero
     do k=0,L
        AxB(k)=A%lmix(N,k)*get_rmix(B%reg,k,j,L)
     end do
     C%less(N,j)=C%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
     !I2:
     AxB(0:)=zero
     do k=1,j
        AxB(k)=A%less(N,k)*get_adv(B%reg,k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
     !I3:
     AxB(0:)=zero
     do k=1,N
        AxB(k)=A%ret(N,k)*B%reg%less(k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
  end do
  !
  ! (t,t`)=>(i,N) <==> Horizontal side, w/ no tip (i=1,N-1)
  C%less(1:N-1,N)=zero
  !    
  if(present(dcoeff))then
     C%lmix(N,0:L) = dcoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=dcoeff*C%less(N,1:N-1)
     C%less(1:N,N)=dcoeff*C%less(1:N,N)
     C%ret(N,1:N)=dcoeff*C%ret(N,1:N)
  endif
  if(present(ccoeff))then
     C%lmix(N,0:L) = ccoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=ccoeff*C%less(N,1:N-1)
     C%less(1:N,N)=ccoeff*C%less(1:N,N)
     C%ret(N,1:N)=ccoeff*C%ret(N,1:N)
  endif
  deallocate(AxB)
end subroutine convolute_kb_contour_gf_sigma


!C(t,t')=[(Ahf*delta + Areg)*B](t,t'), with t=t_max && t'=0,t_max ; t_max_index==N
subroutine convolute_kb_contour_sigma_gf(A,B,C,params,dcoeff,ccoeff)
  type(kb_contour_sigma)              :: A
  type(kb_contour_gf)                 :: B
  type(kb_contour_gf)                 :: C
  type(kb_contour_params)             :: params
  real(8),optional                    :: dcoeff
  complex(8),optional                 :: ccoeff
  integer                             :: N,L
  real(8)                             :: dt,dtau
  complex(8),dimension(:),allocatable :: AxB    
  integer                             :: i,j,k,itau,jtau
  logical                             :: checkA,checkB,checkC
  !
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  if(N==1) stop "convolute_kb_contour_sigma_gf error: called with N=1"
  if(  (.not.A%status).OR.&
       (.not.B%status).OR.&
       (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: A,B,C not allocated"
  checkA=check_dimension_kb_contour(A,params) 
  checkB=check_dimension_kb_contour(B,params)
  checkC=check_dimension_kb_contour(C,params)
  !
  allocate(AxB(0:max(L,N)))
  !
  !Ret. component
  !C^R(t,t')=Ahf(t)*B^R(t,t') + \int_{t'}^t ds Areg^R(t,s)*B ^R(s,t')
  C%ret(N,1:N)=zero
  do j=1,N
     C%ret(N,j) = A%hf(N)*B%ret(N,j)
     AxB(0:)=zero
     do k=j,N
        AxB(k) = A%reg%ret(N,k)*B%ret(k,j)
     enddo
     C%ret(N,j) = C%ret(N,j) + dt*kb_trapz(AxB(0:),j,N)
  enddo
  !
  !Lmix. component
  !C^\lmix(t,tau')= Ahf(t)*B^\lmix(t,tau') +
  !                 \int_0^{beta} ds A^\lmix(t,s)*Breg^M(s,tau') + 
  !                 \int_0^{t} ds A^R(t,s)*Breg^\lmix(s,tau') 
  !               = Ahf*B^\lmix + I1 + I2
  C%lmix(N,0:L)=zero
  do jtau=0,L
     !Ahf*B^\lmix:
     C%lmix(N,jtau) = A%hf(N)*B%lmix(N,jtau)
     !I1: take care of the sign of (tau-tau').
     AxB(0:) = zero
     do k=0,jtau
        AxB(k)=A%reg%lmix(N,k)*B%mats(L+k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
     AxB(0:) = zero
     do k=jtau,L
        AxB(k)=A%reg%lmix(N,k)*B%mats(k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
     !I2:
     AxB(0:) = zero
     do k=1,N
        AxB(k) = A%reg%ret(N,k)*B%lmix(k,jtau)
     enddo
     C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
  enddo
  !
  !Less component
  !C^<(t,t') =  Ahf(t)*B^<(t,t')
  !            -i\int_0^{beta} ds A^\lmix(t,s)*Breg^\rmix(s,t')
  !            + \int_0^{t'} ds A^<(t,s)*Breg^A(s,t')
  !            + \int_0^{t} ds A^R(t,s)*Breg^<(s,t')
  !          = Ahf*B^< + I1 + I2 + I3
  ! (t,t')=>(N,j) <==> Vertical side, with tip (j=1,N)
  C%less(N,1:N)=zero
  do j=1,N
     !Ahf*B^<:
     C%less(N,j)=A%hf(N)*B%less(N,j)
     !I1:
     AxB(0:) = zero
     do k=0,L
        AxB(k)=A%reg%lmix(N,k)*get_rmix(B,k,j,L)
     end do
     C%less(N,j)=C%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
     !I2:
     AxB(0:) = zero
     do k=1,j
        AxB(k)=A%reg%less(N,k)*get_adv(B,k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
     !I3:
     AxB(0:) = zero
     do k=1,N
        AxB(k)=A%reg%ret(N,k)*B%less(k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
  end do
  !
  ! (t,t`)=>(i,N) <==> Horizontal side, w/ no tip (i=1,N-1)
  C%less(1:N-1,N)=zero
  !    
  if(present(dcoeff))then
     C%lmix(N,0:L) = dcoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=dcoeff*C%less(N,1:N-1)
     C%less(1:N,N)=dcoeff*C%less(1:N,N)
     C%ret(N,1:N)=dcoeff*C%ret(N,1:N)
  endif
  if(present(ccoeff))then
     C%lmix(N,0:L) = ccoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=ccoeff*C%less(N,1:N-1)
     C%less(1:N,N)=ccoeff*C%less(1:N,N)
     C%ret(N,1:N)=ccoeff*C%ret(N,1:N)
  endif
  deallocate(AxB)
end subroutine convolute_kb_contour_sigma_gf


!C(t,t') = (a_1*a_2*...*a_n)(t,t')with t=t_max && t'=0,t_max; t_max_index==N
! a_i = contour_gf forall i=1,...,n
subroutine convolute_kb_contour_gf_gf_recursive(A,C,params,dcoeff,ccoeff)
  type(kb_contour_gf)                 :: A(:)
  type(kb_contour_gf)                 :: C
  type(kb_contour_gf)                 :: Knew,Kold
  type(kb_contour_params)             :: params
  real(8),optional                    :: dcoeff
  complex(8),optional                 :: ccoeff
  integer                             :: N,L
  integer                             :: i,j,k,itau,jtau,Na
  logical                             :: checkA,checkC
  !
  Na = size(A)
  !
  if(  (.not.C%status))stop "contour_gf/convolute_kb_contour_gf: C not allocated"
  do i=1,Na
     if(  (.not.A(i)%status))stop "contour_gf/convolute_kb_contour_gf: A(i) not allocated"
  enddo
  !
  call allocate_kb_contour_gf(Knew,params)
  call allocate_kb_contour_gf(Kold,params)
  !
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  if(N==1) stop "convolute_kb_contour_gf_gf_recursive error: called with N=1"
  checkC=check_dimension_kb_contour(C,params%Ntime,L)
  do i=1,Na
     checkA=check_dimension_kb_contour(A(i),params%Ntime,L) 
  enddo
  !
  !Perform the recursive convolution:
  Kold = A(1)
  do i=2,Na-1
     call convolute_kb_contour_gf(Kold,A(i),Knew,params)
     Kold=Knew
  enddo
  call convolute_kb_contour_gf(Knew,A(Na),C,params)
  !
  if(present(dcoeff))then
     C%lmix(N,0:L) = dcoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=dcoeff*C%less(N,1:N-1)
     C%less(1:N,N)=dcoeff*C%less(1:N,N)
     C%ret(N,1:N)=dcoeff*C%ret(N,1:N)
  endif
  if(present(ccoeff))then
     C%lmix(N,0:L) = ccoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=ccoeff*C%less(N,1:N-1)
     C%less(1:N,N)=ccoeff*C%less(1:N,N)
     C%ret(N,1:N)=ccoeff*C%ret(N,1:N)
  endif
  call deallocate_kb_contour_gf(Knew)
  call deallocate_kb_contour_gf(Kold)
  !    
end subroutine convolute_kb_contour_gf_gf_recursive


!C(t,t') = (a_1*a_2*...*a_n)(t,t')with t=t_max && t'=0,t_max; t_max_index==N
! a_i = contour_gf    for i=j_1,...,j_ng
! a_i = contour_sigma for i=k_1,...,k_ns
subroutine convolute_kb_contour_gf_sigma_recursive(G,S,mask,P,params,dcoeff,ccoeff)
  type(kb_contour_gf),dimension(:)    :: G     !array of contour_GF
  type(kb_contour_sigma),dimension(:) :: S     !array of contour_Sigma
  integer,dimension(:)                :: mask  !mask array to select contour_GF (0) over contour_Sigma (1)
  type(kb_contour_gf)                 :: P     !result
  type(kb_contour_params)             :: params
  real(8),optional                    :: dcoeff
  complex(8),optional                 :: ccoeff
  integer                             :: Ng,Ns,Ntot
  type(kb_contour_gf)                 :: Knew,Kold
  integer                             :: N,L
  integer                             :: i,j,k,ig,is
  logical                             :: checkG,checkP,checkS
  !
  Ng   = size(G)
  Ns   = size(S)
  Ntot = size(Mask)
  !
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  if(Ntot/=Ng+Ns)stop "convolute_kb_contour_gf_sigma_recursive error: Ntot != Ns + Ng"
  if(N==1)stop "convolute_kb_contour_gf_sigma_recursive error: call with N=1"
  do i=1,Ng
     if(  (.not.G(i)%status) )stop "convolute_kb_contour_gf_sigma_recursive error: G(i) not allocated"
  enddo
  do i=1,Ns
     if(  (.not.S(i)%status) )stop "convolute_kb_contour_gf_sigma_recursive error: S(i) not allocated"
  enddo
  if(  (.not.P%status) )stop "convolute_kb_contour_gf_sigma_recursive error: P not allocated"
  do i=1,Ng
     checkG=check_dimension_kb_contour(G(i),params) 
  enddo
  do i=1,Ns
     checkS=check_dimension_kb_contour(S(i),params) 
  enddo
  checkP=check_dimension_kb_contour(P,params)
  !
  call allocate_kb_contour_gf(Knew,params)
  call allocate_kb_contour_gf(Kold,params)
  !
  !check that mask does not contain 2 consecutive ones:
  j=mask(1)
  do i=2,Ntot
     if(j*mask(i)/=0)then
        print*,"convolute_kb_contour_gf_sigma_recursive error: mask has 2 consecutives ones at",i
        stop
     endif
     j=mask(i)
  enddo
  !
  !Perform the recursive convolution:
  !get the first product P1 = g1/s1 * g2/s2 (can not be s1*s2 though)
  !this step features 3 cases: i) g1*g2, ii) g1*s2, iii) s1*g2
  is=1
  ig=1
  if(mask(1)==0.AND.mask(2)==0)then
     ig=ig+1
     call convolute_kb_contour_gf(G(1),G(ig),Kold,params)
  elseif(mask(1)==1.AND.mask(2)==0)then
     ig=ig+1
     call convolute_kb_contour_gf(S(1),G(ig),Kold,params)
  elseif(mask(1)==0.AND.mask(2)==1)then
     is=is+1
     call convolute_kb_contour_gf(G(1),S(is),Kold,params)
  else
     print*,"convolute_kb_contour_gf_sigma_recursive error: mask(1)==mask(2)==1"
     stop
  endif
  !
  do i=2,Ntot-1
     if(mask(i)==0)then
        ig=ig+1
        call convolute_kb_contour_gf(Kold,G(ig),Knew,params)
     elseif(mask(i)==1)then
        is=is+1
        call convolute_kb_contour_gf(Kold,S(is),Knew,params)
     else
        print*,"convolute_kb_contour_gf_sigma_recursive error: mask(i)!=1.OR.0"
        stop
     endif
     Kold = Knew
  enddo
  if(mask(Ntot)==0)then
     ig=ig+1
     call convolute_kb_contour_gf(Kold,G(ig),P,params)
  elseif(mask(Ntot)==1)then
     is=is+1
     call convolute_kb_contour_gf(Kold,S(is),P,params)
  else
     print*,"convolute_kb_contour_gf_sigma_recursive error: mask(Ntot)!=1.OR.0"
     stop
  endif
  !
  if(ig/=Ng.OR.is/=Ns)stop "convolute_kb_contour_gf_sigma_recursive error: ig!=Ng OR is!=Ns"
  !
  if(present(dcoeff))then
     P%lmix(N,0:L) = dcoeff*P%lmix(N,0:L)
     P%less(N,1:N-1)=dcoeff*P%less(N,1:N-1)
     P%less(1:N,N)=dcoeff*P%less(1:N,N)
     P%ret(N,1:N)=dcoeff*P%ret(N,1:N)
  endif
  if(present(ccoeff))then
     P%lmix(N,0:L) = ccoeff*P%lmix(N,0:L)
     P%less(N,1:N-1)=ccoeff*P%less(N,1:N-1)
     P%less(1:N,N)=ccoeff*P%less(1:N,N)
     P%ret(N,1:N)=ccoeff*P%ret(N,1:N)
  endif
  call deallocate_kb_contour_gf(Knew)
  call deallocate_kb_contour_gf(Kold)
  !    
end subroutine convolute_kb_contour_gf_sigma_recursive



!C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max; t_max_index==N
! The case A(t,t')=A(t)delta(t,t')
subroutine convolute_kb_contour_gf_delta_left(A,B,C,params)
  type(kb_contour_gf)                 :: B,C
  complex(8),dimension(:)             :: A
  type(kb_contour_params)             :: params
  integer                             :: N,L
  integer                             :: i,j,k,itau,jtau
  logical                             :: checkB,checkC
  if(  (.not.B%status).OR.&
       (.not.C%status))stop "convolute_kb_contour_gf_delta_left error: B,C not allocated"
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  if(size(A)<N)stop "convolute_kb_contour_gf_delta_left error: size(A) < N"
  checkB=check_dimension_kb_contour(B,params%Ntime,L) 
  checkC=check_dimension_kb_contour(C,params%Ntime,L)
  !
  ! Matsubara
  if (N==1)C%mats = A(1)*B%mats
  !Ret. component
  do j=1,N
     C%ret(n,j) = A(n)*B%ret(n,j)
  enddo
  !Less. component
  do j=1,N
     C%less(n,j) = A(n)*B%less(n,j)
  enddo
  !Lmix. component
  do jtau=0,L
     C%lmix(n,jtau) = A(n)*B%lmix(n,jtau)
  enddo
end subroutine convolute_kb_contour_gf_delta_left


!C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max; t_max_index==N
! The case B(t,t')=B(t)delta(t,t')
subroutine convolute_kb_contour_gf_delta_right(A,B,C,params)
  type(kb_contour_gf)                 :: A,C
  complex(8),dimension(:)             :: B
  type(kb_contour_params)             :: params
  integer                             :: N,L
  integer                             :: i,j,k,itau,jtau
  logical                             :: checkA,checkC
  if(  (.not.A%status).OR.&
       (.not.C%status))stop "convolute_kb_contour_gf_delta_right error: A,C not allocated"
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  checkA=check_dimension_kb_contour(A,params%Ntime,L)
  if(size(B)<N)stop "convolute_kb_contour_gf_delta_right error: size(B) < N"
  checkC=check_dimension_kb_contour(C,params%Ntime,L)
  !
  ! Matsubara
  if (N==1)C%mats(0:) = A%mats(0:)*B(1)
  !Ret. component
  do j=1,N
     C%ret(n,j) = A%ret(n,j)*B(j)
  enddo
  !Less. component
  do j=1,N
     C%less(n,j) = A%less(n,j)*B(j)
  enddo
  !Lmix. component
  do jtau=0,L
     C%lmix(N,jtau) = A%lmix(N,jtau)*B(1)
  enddo
end subroutine convolute_kb_contour_gf_delta_right
