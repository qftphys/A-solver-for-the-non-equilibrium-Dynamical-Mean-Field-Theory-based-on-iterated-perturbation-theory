MODULE NEQ_MEASURE
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_INPUT_VARS
  USE ELECTRIC_FIELD
  !
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only:free_unit
  implicit none

  private

  interface measure_ekin
     module procedure measure_ekin_hk
     module procedure measure_ekin_bethe
  end interface measure_ekin

  interface measure_epot
     module procedure measure_epot_normal
     module procedure measure_epot_superc
     module procedure measure_epot_ntdocc
  end interface measure_epot

  interface measure_docc
     module procedure measure_docc_normal
     module procedure measure_docc_superc
  end interface measure_docc

  interface measure_observables
     module procedure measure_observables_normal
     module procedure measure_observables_normal_bethe
     module procedure measure_observables_superc
  end interface measure_observables

  public :: measure_observables
  public :: measure_dens
  public :: measure_docc
  public :: measure_delta
  public :: measure_ekin
  public :: measure_epot
  public :: measure_current

contains



  !+-------------------------------------------------------------------+
  !PURPOSE: measure some observables and print them
  !+-------------------------------------------------------------------+
  subroutine measure_observables_normal(g,self,Gk,Hk,Wtk,params)
    type(kb_contour_gf)       :: g
    type(kb_contour_gf)       :: self
    type(kb_contour_gf)       :: gk(:)
    complex(8),dimension(:,:) :: Hk
    real(8)                   :: wtk(size(gk))
    type(kb_contour_params)   :: params
    integer                   :: unit,itime,Lk
    real(8)                   :: dens,docc,ekin,epot,etot
    Lk=size(Hk,2)
    if(size(Hk,1)<params%Ntime)stop "measure_observables: size(Hk),1 is not Ntime"
    if(size(Gk)/=Lk)           stop "measure_observables: size(Gk) is not Lk"
    itime = params%Nt
    unit = free_unit()
    open(unit,file="observables.info")
    write(unit,"(8A20)")"time","N","Docc","Epot","Ekin","Etot"
    close(unit)
    open(unit,file="observables.neqipt",position="append")
    dens = measure_dens(g,params)
    docc = measure_docc_normal(g,self,params)
    epot = measure_epot(dens,docc,params)
    ekin = measure_ekin_hk(Gk,Hk,Wtk,params)
    write(unit,"(6F20.12)")params%t(itime),dens,docc,epot,ekin,ekin+epot
    close(unit)
  end subroutine measure_observables_normal

  subroutine measure_observables_normal_bethe(g,self,params)
    type(kb_contour_gf)                :: g
    type(kb_contour_gf)                :: self
    type(kb_contour_params)            :: params
    integer                            :: unit,itime,Lk
    real(8)                            :: dens,docc,ekin,epot,etot
    itime = params%Nt
    unit = free_unit()
    open(unit,file="observables.info")
    write(unit,"(8A20)")"time","N","Docc","Epot","Ekin","Etot"
    close(unit)
    open(unit,file="observables.neqipt",position="append")
    dens = measure_dens(g,params)
    docc = measure_docc_normal(g,self,params)
    epot = measure_epot(dens,docc,params)
    ekin = measure_ekin_bethe(G,params)
    write(unit,"(6F20.12)")params%t(itime),dens,docc,epot,ekin,ekin+epot
    close(unit)
  end subroutine measure_observables_normal_bethe

  subroutine measure_observables_superc(g,self,Gk,Hk,Wtk,params)
    type(kb_contour_gf)       :: g(2,2)
    type(kb_contour_gf)       :: self(2,2)
    type(kb_contour_gf)       :: gk(:,:)             ![2][Lk]
    complex(8),dimension(:,:) :: Hk                  ![Ntime][Lk]
    real(8)                   :: wtk(size(gk,2))     ![Lk]
    type(kb_contour_params)   :: params
    integer                   :: unit,itime,Lk
    real(8)                   :: dens,docc,delta,ekin,epot,etot
    !
    Lk = size(Hk,2)
    if(size(Hk,1)<params%Ntime)stop "measure_observables_superc: size(Hk),1 is not Ntime"
    if(size(Gk,1)/=2)          stop "measure_observables_superc: size(Gk,2) is not 2"
    if(size(Gk,2)/=Lk)         stop "measure_observables_superc: size(Gk,2) is not 2"
    !
    itime = params%Nt
    unit = free_unit()
    open(unit,file="observables.info")
    write(unit,"(8A20)")"time","N","Docc","Delta","Epot","Ekin","Etot"
    close(unit)
    open(unit,file="observables.neqipt",position="append")
    dens = measure_dens(g(1,1),params)
    docc = measure_docc_superc(g,self,params)
    delta= measure_delta(g(1,2),params)
    epot = measure_epot(dens,docc,params)
    ekin = measure_ekin_hk(Gk(1,1),Hk,Wtk,params)
    write(unit,"(6F20.12)")params%t(itime),dens,delta,docc,epot,ekin,ekin+epot
    close(unit)
  end subroutine measure_observables_superc










  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  !Lattice: -i\sum_ks epsik(k,t)G_k(t,t)= 2*\sum_k e(k,t)n(k,t)
  !Bethe: E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_ekin_hk(Gk,Hk,Wtk,params) result(ekin)
    type(kb_contour_gf),dimension(:)    :: gk
    complex(8),dimension(:,:)           :: Hk
    real(8),dimension(:)                :: Wtk
    type(kb_contour_params)             :: params
    real(8)                             :: Ekin
    integer                             :: i,ik,Lk,N
    complex(8),dimension(:),allocatable :: Ker
    real(8)                             :: nkt
    N  = params%Nt
    Lk = size(Gk)
    if(size(Hk,1)<params%Ntime)stop "measure_ekin: size(Hk),1 is not Ntime"
    if(size(Hk,2)<Lk)          stop "measure_ekin: size(Hk),2 is not Lk"
    if(size(Wtk)<Lk)           stop "measure_ekin: size(Wt),1 is not Lk"
    !
    Ekin=0d0
    do ik=1,Lk
       Nkt  = dimag(Gk(ik)%less(N,N))
       Ekin = Ekin + 2d0*Wtk(ik)*Nkt*Hk(N,ik)
    enddo
  end function measure_ekin_hk

  function measure_ekin_bethe(g,params) result(ekin)
    type(kb_contour_gf)                 :: g
    type(kb_contour_params)             :: params
    real(8)                             :: ekin
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: Ker
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    allocate(Ker(0:max(N,L)))
    if(N==1)then
       do k=0,L
          Ker(k)=G%mats(L-k)*G%mats(k)
       end do
       ekin = -2.d0*params%dtau*kb_trapz(Ker(0:),0,L)
    else
       do k=0,L
          Ker(k)=G%lmix(N,k)*conjg(G%lmix(N,L-k))
       end do
       ekin=2.d0*params%dtau*dimag( (-xi)*kb_trapz(Ker(0:),0,L) )
       do k=1,N
          Ker(k)=G%ret(N,k)*G%less(k,N)
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
       do k=1,N
          Ker(k)=G%less(N,k)*conjg(G%ret(N,k))
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
    endif
    deallocate(Ker)
  end function measure_ekin_bethe






  !+-------------------------------------------------------------------+
  !PURPOSE: measure current
  !+-------------------------------------------------------------------+
  subroutine measure_current(Gk,Vkt,Wtk,params)
    type(kb_contour_gf)            :: gk(:)
    real(8),dimension(:,:,:)       :: Vkt
    real(8),dimension(:)           :: wtk
    type(kb_contour_params)        :: params
    integer                        :: unit,itime,Lk,ik,i
    real(8),dimension(size(Vkt,3)) :: Jloc
    real(8)                        :: nkt
    !
    Lk=size(gk)
    itime = params%Nt
    !
    if(size(Vkt,1)<params%Ntime)stop "neq_measure_current: dim(Vkt,2) < Ntime"
    if(size(Vkt,2)/=Lk)stop "neq_measure_current: dim(Vkt,3) != Lk"
    !
    unit = free_unit()
    open(unit,file="current.info")
    write(unit,"(8A20)")"time","Jx","Jy","Jz"
    close(unit)
    !
    open(unit,file="current.neqipt",position="append")
    Jloc=0d0
    do ik=1,Lk
       nkt  = dimag(Gk(ik)%less(itime,itime))
       Jloc = Jloc + Wtk(ik)*Nkt*Vkt(itime,ik,:)
    enddo
    write(unit,"(4F20.12)")params%t(itime),(Jloc(i),i=1,size(Jloc))
    close(unit)
  end subroutine measure_current




  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the density at a given istant of time
  ! n(t)=-xi*G^<(t,t)
  !+-------------------------------------------------------------------+
  function measure_dens(g,params) result(dens)
    type(kb_contour_gf)                 :: g
    type(kb_contour_params)             :: params
    real(8)                             :: dens
    integer                             :: N
    N = params%Nt
    dens = dimag(G%less(N,N))
  end function measure_dens


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the order parameter at a given istant of time
  ! \Delta(t)=-xi*F^<(t,t)
  !+-------------------------------------------------------------------+
  function measure_delta(f,params) result(delta)
    type(kb_contour_gf)                 :: f
    type(kb_contour_params)             :: params
    real(8)                             :: delta
    integer                             :: N
    N = params%Nt
    delta = dimag(F%less(N,N))
  end function measure_delta


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the double occupancy at a given istant of time
  ! d(t)=n_up(t)*n_do(t)-1/U0*[S^M*G^M]
  !      n_up(t)*n_do(t)-i/U*[S^R*G^< + S^<*G^A + S^\lmix*G^\rmix](t,t)
  !
  ! d(t)=n_up(t)*n_do(t)-1/U0*[S_11^M*G_11^M + S_12^M*G_21^M ]
  !      n_up(t)*n_do(t)-i/U*[S_11^R*G_11^< + S_11^<*G_11^A + S_11^\lmix*G_11^\rmix](t,t)
  !                     -i/U*[S_12^R*G_21^< + S_12^<*G_21^A + S_12^\lmix*G_21^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_docc_normal(g,self,params) result(docc)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: docc
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: SxG
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    nt   = dimag(G%less(N,N))
    allocate(SxG(0:max(N,L)))
    docc = nt**2
    if(N==1)then
       if(ui/=0.d0)then
          do k=0,L
             SxG(k)=Self%mats(L-k)*G%mats(k)
          end do
          docc=docc-1.d0/Ui*params%dtau*kb_trapz(SxG(0:),0,L)
       endif
    else
       if(u/=0.d0)then
          do k=0,L
             SxG(k)=Self%lmix(N,k)*conjg(G%lmix(N,L-k))
          end do
          docc=docc + 1.d0/U*params%dtau*dimag( (-xi)*kb_trapz(SxG(0:),0,L) )
          do k=1,N
             SxG(k)=Self%ret(N,k)*G%less(k,N)
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
          do k=1,N
             SxG(k)=Self%less(N,k)*conjg(G%ret(N,k))
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
       endif
    endif
    deallocate(SxG)
  end function measure_docc_normal

  function measure_docc_superc(g,self,params) result(docc)
    type(kb_contour_gf)                 :: g(2,2)
    type(kb_contour_gf)                 :: self(2,2)
    type(kb_contour_params)             :: params
    real(8)                             :: docc
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: SxG
    real(8)                             :: nt
    !
    N = params%Nt
    L = params%Ntau
    !
    nt   = dimag(G(1,1)%less(N,N))
    !
    allocate(SxG(0:max(N,L)))
    docc = nt**2
    if(N==1)then
       if(ui/=0.d0)then
          do k=0,L
             SxG(k)=Self(1,1)%mats(L-k)*G(1,1)%mats(k) + Self(1,2)%mats(L-k)*G(2,1)%mats(k)
          end do
          docc=docc-1.d0/Ui*params%dtau*kb_trapz(SxG(0:),0,L)
       endif
    else
       if(u/=0.d0)then
          do k=0,L
             SxG(k)=Self(1,1)%lmix(N,k)*get_rmix(G(1,1),k,N,L) + Self(1,2)%lmix(N,k)*get_rmix(G(2,1),k,N,L)
          end do
          docc=docc + 1.d0/U*params%dtau*dimag( (-xi)*kb_trapz(SxG(0:),0,L) )
          do k=1,N
             SxG(k)=Self(1,1)%ret(N,k)*G(1,1)%less(k,N) + Self(1,2)%ret(N,k)*G(2,1)%less(k,N)
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
          do k=1,N
             SxG(k)=Self(1,1)%less(N,k)*get_adv(G(1,1),k,N) + Self(1,2)%less(N,k)*get_adv(G(2,1),k,N)
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
       endif
    endif
    deallocate(SxG)
  end function measure_docc_superc





  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! U(t)= U*docc(t) - n(t) + 1/4
  !+-------------------------------------------------------------------+
  function measure_epot_normal(g,self,params) result(epot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: epot,docc,nt
    integer                             :: i,k,j,N,L
    N = params%Nt
    L = params%Ntau
    !
    if(N==1)then
       nt   = measure_dens(g,self,params)
       docc = measure_docc(g,self,params)
       epot = abs(Ui)*(docc - nt + 0.25d0)
    else
       nt   = measure_dens(g,self,params)
       docc = measure_docc(g,self,params)
       epot = abs(U)*(docc - nt + 0.25d0)
    endif
  end function measure_epot_normal

  function measure_epot_superc(g,self,params) result(epot)
    type(kb_contour_gf)                 :: g(2,2)
    type(kb_contour_gf)                 :: self(2,2)
    type(kb_contour_params)             :: params
    real(8)                             :: epot,docc,nt
    integer                             :: i,k,j,N,L
    !
    N = params%Nt
    L = params%Ntau
    !
    if(N==1)then
       nt   = measure_dens(g(1,1),self(1,1),params)
       docc = measure_docc(g,self,params)
       epot = abs(Ui)*(docc - nt + 0.25d0)
    else
       nt   = measure_dens(g(1,1),self(1,1),params)
       docc = measure_docc(g,self,params)
       epot = abs(U)*(docc - nt + 0.25d0)
    endif
  end function measure_epot_superc

  function measure_epot_ntdocc(nt,docc,params) result(epot)
    real(8)                             :: nt,docc
    type(kb_contour_params)             :: params
    real(8)                             :: epot
    integer                             :: i,k,j,N,L
    !
    N = params%Nt
    L = params%Ntau
    !
    if(N==1)then
       epot = abs(Ui)*(docc - nt + 0.25d0)
    else
       epot = abs(U)*(docc - nt + 0.25d0)
    endif
  end function measure_epot_ntdocc






END MODULE NEQ_MEASURE
