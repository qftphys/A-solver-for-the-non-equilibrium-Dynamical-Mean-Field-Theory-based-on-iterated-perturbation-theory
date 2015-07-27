!C(t,t')=[A*(Bhf*delta + Breg)](t,t'), with t=t_max && t'=0,t_max ; t_max_index==N
subroutine convolute_kb_contour_sigma_simple(A,B,C,params,dcoeff,ccoeff)
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
     AxB(0:)  = zero
     do k=j,N
        AxB(k) = A%ret(N,k)*B%reg%ret(k,j)
     enddo
     C%ret(N,j) = C%ret(N,j) + dt*kb_trapz(AxB(0:),j,N) + A%ret(N,j)*B%hf(j)
  enddo
  !
  !Lmix. component
  !C^\lmix(t,tau')= A^\lmix(t,tau')*Bhf(0+) +
  !                 \int_0^{beta} ds A^\lmix(t,s)*Breg^M(s,tau') + 
  !                 \int_0^{t} ds A^R(t,s)*Breg^\lmix(s,tau') 
  !               = A*Bhf + I1 + I2
  C%lmix(N,0:L)=zero
  do jtau=0,L
     AxB(0:) = zero
     !I1: take care of the sign of (tau-tau').
     do k=0,jtau
        AxB(k)=A%lmix(N,k)*B%reg%mats(L+k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
     do k=jtau,L
        AxB(k)=A%lmix(N,k)*B%reg%mats(k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
     !I2:
     AxB(0:) = zero
     do k=1,N
        AxB(k) = A%ret(N,k)*B%reg%lmix(k,jtau)
     enddo
     C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
     !A*Bhf:
     C%lmix(N,jtau) = C%lmix(N,jtau) + A%lmix(N,jtau)*B%hf(1)
  enddo
  !
  !Less component
  !C^<(t,t') =  A^<(t,t')*conjg(Bhf(t'))
  !            -i\int_0^{beta} ds A^\lmix(t,s)*Breg^\rmix(s,t')
  !            + \int_0^{t'} ds A^<(t,s)*Breg^A(s,t')
  !            + \int_0^{t} ds A^R(t,s)*Breg^<(s,t')
  !          = A^<*Bhf* + I1 + I2 + I3
  ! (t,t')=>(N,j) <==> Vertical side, with tip (j=1,N)
  do j=1,N
     C%less(N,j)=zero
     !I1:
     do k=0,L
        AxB(k)=A%lmix(N,k)*get_rmix(B%reg,k,j,L)
     end do
     C%less(N,j)=C%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
     !I2:
     do k=1,j
        AxB(k)=A%less(N,k)*get_adv(B%reg,k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
     !I3:
     do k=1,N
        AxB(k)=A%ret(N,k)*B%reg%less(k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
     !A^<*Bhf*:
     C%less(N,j)=C%less(N,j)+A%less(N,j)*conjg(B%hf(j))
  end do
  !
  ! (t,t`)=>(i,N) <==> Horizontal side, w/ no tip (i=1,N-1)
  do i=1,N-1
     C%less(i,N)=zero
  enddo
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
end subroutine convolute_kb_contour_sigma_simple
