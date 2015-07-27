subroutine allocate_kb_contour_gf_main(G,params)
  type(kb_contour_gf)     :: G
  type(kb_contour_params) :: params
  integer                 :: i,j,N,L,Lf
  if(allocated(G%less))deallocate(G%less)
  if(allocated(G%ret)) deallocate(G%ret)
  if(allocated(G%lmix))deallocate(G%lmix)
  if(allocated(G%mats))deallocate(G%mats)
  if(allocated(G%iw))deallocate(G%iw)
  N = params%Ntime            !<== allocate at maximum time
  L = params%Ntau
  Lf= params%Niw
  allocate(G%less(N,N))  ; G%less=zero
  allocate(G%ret(N,N))   ; G%ret=zero
  allocate(G%lmix(N,0:L)); G%lmix=zero
  allocate(G%mats(0:L))  ; G%mats=0d0
  allocate(G%iw(Lf))     ; G%iw=zero
  G%status=.true.
end subroutine allocate_kb_contour_gf_main
!
subroutine allocate_kb_contour_gf_nambu_redux(G,params)
  type(kb_contour_gf),dimension(2) :: G
  type(kb_contour_params)          :: params
  call allocate_kb_contour_gf_main(G(1),params)
  call allocate_kb_contour_gf_main(G(2),params)
  G(1)%anomalous=.false.
  G(2)%anomalous=.true.
end subroutine allocate_kb_contour_gf_nambu_redux
!
subroutine allocate_kb_contour_gf_nambu(G,params)
  type(kb_contour_gf),dimension(2,2) :: G
  type(kb_contour_params)            :: params
  integer                            :: i,j
  do i=1,2
     do j=1,2
        call allocate_kb_contour_gf_main(G(i,j),params)
     enddo
  enddo
  G(1,1)%anomalous=.false.
  G(1,2)%anomalous=.true.
  G(2,1)%anomalous=.true.
  G(2,2)%anomalous=.false.
end subroutine allocate_kb_contour_gf_nambu
!
!
!
subroutine allocate_kb_contour_sigma_main(S,params)
  type(kb_contour_sigma)  :: S
  type(kb_contour_params) :: params
  integer                 :: N
  if(S%status)call deallocate_kb_contour_sigma(S)
  N = params%Ntime              !<== allocate at maximum time
  call allocate_kb_contour_gf_main(S%reg,params)
  allocate(S%hf(N))
  S%status=.true.
end subroutine allocate_kb_contour_sigma_main
!
subroutine allocate_kb_contour_sigma_nambu_redux(S,params)
  type(kb_contour_sigma),dimension(2) :: S
  type(kb_contour_params)             :: params
  call allocate_kb_contour_sigma_main(S(1),params)
  call allocate_kb_contour_sigma_main(S(2),params)
  S(1)%reg%anomalous=.false.
  S(2)%reg%anomalous=.true.
end subroutine allocate_kb_contour_sigma_nambu_redux
!
subroutine allocate_kb_contour_sigma_nambu(S,params)
  type(kb_contour_sigma),dimension(2,2) :: S
  type(kb_contour_params)               :: params
  integer                               :: i,j
  do i=1,2
     do j=1,2
        call allocate_kb_contour_sigma_main(S(i,j),params)
     enddo
  enddo
  S(1,1)%reg%anomalous=.false.
  S(1,2)%reg%anomalous=.true.
  S(2,1)%reg%anomalous=.true.
  S(2,2)%reg%anomalous=.false.
end subroutine allocate_kb_contour_sigma_nambu
!
!
!
subroutine allocate_kb_contour_dgf_main(dG,params,wgtr)
  type(kb_contour_dgf)    :: dG
  type(kb_contour_params) :: params
  integer                 :: i,j,N,L
  logical,optional        :: wgtr
  if(allocated(dG%less))deallocate(dG%less)
  if(allocated(dG%ret)) deallocate(dG%ret)
  if(allocated(dG%lmix))deallocate(dG%lmix)
  N=params%Ntime           !<== allocate at maximum time
  L=params%Ntau
  allocate(dG%less(N))  ; dG%less=zero
  allocate(dG%ret(N))   ; dG%ret=zero
  allocate(dG%lmix(0:L)); dG%lmix=zero
  if(present(wgtr).AND.wgtr)then
     allocate(dG%gtr(N))
     dG%gtr=zero
  endif
  dG%status=.true.
end subroutine allocate_kb_contour_dgf_main
!
subroutine allocate_kb_contour_dgf_nambu_redux(dG,params)
  type(kb_contour_dgf),dimension(2)     :: dG
  type(kb_contour_params) :: params
  call allocate_kb_contour_dgf_main(dG(1),params)
  call allocate_kb_contour_dgf_main(dG(2),params)
  dG(1)%anomalous=.false.
  dG(2)%anomalous=.true.
end subroutine allocate_kb_contour_dgf_nambu_redux
!
subroutine allocate_kb_contour_dgf_nambu(dG,params)
  type(kb_contour_dgf),dimension(2,2)     :: dG
  type(kb_contour_params) :: params
  integer :: i,j
  do i=1,2
     do j=1,2
        call allocate_kb_contour_dgf_main(dG(1,j),params)
     enddo
  enddo
  dG(1,1)%anomalous=.false.
  dG(1,2)%anomalous=.true.
  dG(2,1)%anomalous=.true.
  dG(2,2)%anomalous=.false.
end subroutine allocate_kb_contour_dgf_nambu





! DEALLOCATE
subroutine deallocate_kb_contour_gf_main(G)
  type(kb_contour_gf) :: G
  if(.not.G%status)stop "contour_gf/deallocate_kb_contour_gf: G not allocated"
  deallocate(G%less,G%ret,G%lmix,G%mats,G%iw)
  G%status=.false.
end subroutine deallocate_kb_contour_gf_main
!
subroutine deallocate_kb_contour_gf_nambu(G)
  type(kb_contour_gf),dimension(2,2) :: G
  integer :: i,j
  do i=1,2
     do j=1,2
        call deallocate_kb_contour_gf_main(G(i,j))
        G(i,j)%anomalous=.false.
     enddo
  enddo
end subroutine deallocate_kb_contour_gf_nambu
!
subroutine deallocate_kb_contour_gf_nambu_redux(G)
  type(kb_contour_gf),dimension(2) :: G
  integer :: i,j
  do i=1,2
     call deallocate_kb_contour_gf_main(G(i))
     G(i)%anomalous=.false.
  enddo
end subroutine deallocate_kb_contour_gf_nambu_redux
!
!
!
subroutine deallocate_kb_contour_sigma_main(S)
  type(kb_contour_sigma) :: S
  if(.not.S%status)stop "contour_gf/deallocate_kb_contour_sigma: S not allocated"
  if(S%reg%status)call deallocate_kb_contour_gf_main(S%reg)
  if(allocated(S%hf)) deallocate(S%hf)
  S%status=.false.
end subroutine deallocate_kb_contour_sigma_main
!
subroutine deallocate_kb_contour_sigma_nambu(S)
  type(kb_contour_sigma),dimension(2,2) :: S
  integer                               :: i,j
  do i=1,2
     do j=1,2
        call deallocate_kb_contour_sigma_main(S(i,j))
        S(i,j)%reg%anomalous=.false.
     enddo
  enddo
end subroutine deallocate_kb_contour_sigma_nambu
!
subroutine deallocate_kb_contour_sigma_nambu_redux(S)
  type(kb_contour_sigma),dimension(2) :: S
  integer                             :: i
  do i=1,2
     call deallocate_kb_contour_sigma_main(S(i))
     S(i)%reg%anomalous=.false.
  enddo
end subroutine deallocate_kb_contour_sigma_nambu_redux
!
!
!
subroutine deallocate_kb_contour_dgf_main(dG)
  type(kb_contour_dgf) :: dG
  if(.not.dG%status)stop "contour_gf/deallocate_kb_contour_dgf: dG not allocated"
  deallocate(dG%less,dG%ret,dG%lmix)
  if(allocated(dG%gtr))deallocate(dG%gtr)
  dG%status=.false.
end subroutine deallocate_kb_contour_dgf_main
!
subroutine deallocate_kb_contour_dgf_nambu(dG)
  type(kb_contour_dgf),dimension(2,2) :: dG
  integer                             :: i,j
  do i=1,2
     do j=1,2
        call deallocate_kb_contour_dgf_main(dG(i,j))
        dG(i,j)%anomalous=.false.
     enddo
  enddo
end subroutine deallocate_kb_contour_dgf_nambu
!
subroutine deallocate_kb_contour_dgf_nambu_redux(dG)
  type(kb_contour_dgf),dimension(2) :: dG
  integer                           :: i
  do i=1,2
     call deallocate_kb_contour_dgf_main(dG(i))
     dG(i)%anomalous=.false.
  enddo
end subroutine deallocate_kb_contour_dgf_nambu_redux
