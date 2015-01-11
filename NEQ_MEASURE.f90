MODULE NEQ_MEASURE
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_INPUT_VARS
  USE ELECTRIC_FIELD
  USE CONSTANTS
  USE IOTOOLS, only:free_unit
  USE TIGHT_BINDING
  implicit none

  private

  public :: neq_measure_observables
  public :: neq_measure_current
  public :: neq_measure_dens
  public :: neq_measure_docc
  public :: neq_measure_ekin
  public :: neq_measure_epot
  public :: neq_measure_etot


contains


  !+-------------------------------------------------------------------+
  !PURPOSE: measure some observables and print them
  !+-------------------------------------------------------------------+
  subroutine neq_measure_observables(g,self,params)
    type(kb_contour_gf)     :: g
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    integer                 :: unit,itime
    real(8)                 :: dens,docc,ekin,epot,etot
    itime = params%Itime
    unit = free_unit()
    open(unit,file="observables.info")
    write(unit,"(8A20)")"time","N","Docc","Ekin","Epot","Etot"
    close(unit)
    open(unit,file="observables.ipt",position="append")
    dens = neq_measure_dens(g,self,params)
    docc = neq_measure_docc(g,self,params)
    ekin = neq_measure_ekin(g,self,params)
    epot = neq_measure_epot(g,self,params)
    etot = ekin + epot
    write(unit,"(6F20.12)")params%t(itime),dens,docc,ekin,epot,etot
    close(unit)
  end subroutine neq_measure_observables



  !+-------------------------------------------------------------------+
  !PURPOSE: measure current
  !+-------------------------------------------------------------------+
  subroutine neq_measure_current(Gk,Vkt,Wtk,params)
    type(kb_contour_gf)            :: gk(:)
    real(8),dimension(:,:,:)       :: Vkt
    real(8),dimension(:)           :: wtk
    type(kb_contour_params)        :: params
    integer                        :: unit,itime,Lk,ik,i
    real(8),dimension(size(Vkt,1)) :: Kt,Ak,Jloc
    real(8)                        :: nkt
    !
    Lk=size(gk)
    itime = params%Itime
    !
    if(size(Vkt,2)<params%Ntime)stop "neq_measure_current: dim(Vkt,2) < Ntime"
    if(size(Vkt,3)/=Lk)stop "neq_measure_current: dim(Vkt,3) != Lk"
    !
    unit = free_unit()
    open(unit,file="current.info")
    write(unit,"(8A20)")"time","Jx","Jy","Jz"
    close(unit)
    !
    open(unit,file="current.ipt",position="append")
    Jloc=0d0
    do ik=1,Lk
       nkt  = dimag(Gk(ik)%less(itime,itime))
       Jloc(:)  = Jloc(:) + Wtk(ik)*Nkt*Vkt(:,itime,ik)
    enddo
    write(unit,"(4F20.12)")params%t(itime),(Jloc(i),i=1,size(Jloc))
    close(unit)
  end subroutine neq_measure_current



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the density at a given istant of time
  ! n(t)=-xi*G^<(t,t)
  !+-------------------------------------------------------------------+
  function neq_measure_dens(g,self,params) result(dens)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: dens
    integer                             :: N
    N = params%Itime
    dens = dimag(G%less(N,N))
  end function neq_measure_dens


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the double occupancy at a given istant of time
  ! d(t)=n_up(t)*n_do(t)-1/U0*[Self^M*G^M]
  !      n_up(t)*n_do(t)-i/U*[Self^R*G^< + Self^<*G^A + Self^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function neq_measure_docc(g,self,params) result(docc)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: docc
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: SxG
    real(8)                             :: nt
    N = params%Itime
    L = params%Ntau
    !
    nt   = dimag(G%less(N,N))
    allocate(SxG(max(N,L)))
    docc = nt**2
    if(N==1)then
       if(U0/=0d0)then
          do k=1,L
             SxG(k)=Self%mats(L-k+1)*G%mats(k)
          end do
          docc=docc-1.d0/U0*params%dtau*kb_trapz(SxG,1,L)
       endif
    else
       if(uloc/=0.d0)then
          do k=1,L
             SxG(k)=Self%lmix(N,k)*conjg(G%lmix(N,L-k+1))
          end do
          docc=docc + 1.d0/Uloc*params%dtau*dimag( (-xi)*kb_trapz(SxG,1,L) )
          do k=1,N
             SxG(k)=Self%ret(N,k)*G%less(k,N)
          end do
          docc=docc + 1.d0/Uloc*params%dt*dimag(kb_trapz(SxG,1,N))
          do k=1,N
             SxG(k)=Self%less(N,k)*conjg(G%ret(N,k))
          end do
          docc=docc + 1.d0/Uloc*params%dt*dimag(kb_trapz(SxG,1,N))
       endif
    endif
    deallocate(SxG)
  end function neq_measure_docc


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function neq_measure_ekin(g,self,params) result(ekin)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: ekin
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: Ker
    real(8)                             :: nt
    N = params%Itime
    L = params%Ntau
    !
    allocate(Ker(max(N,L)))
    if(N==1)then
       do k=1,L
          Ker(k)=G%mats(L-k+1)*G%mats(k)
       end do
       ekin = -2.d0*params%dtau*kb_trapz(Ker,1,L)
    else
       do k=1,L
          Ker(k)=G%lmix(N,k)*conjg(G%lmix(N,L-k+1))
       end do
       ekin=2.d0*params%dtau*dimag( (-xi)*kb_trapz(Ker,1,L) )
       do k=1,N
          Ker(k)=G%ret(N,k)*G%less(k,N)
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker,1,N))
       do k=1,N
          Ker(k)=G%less(N,k)*conjg(G%ret(N,k))
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker,1,N))
    endif
    deallocate(Ker)
  end function neq_measure_ekin



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! U(t)= U*docc(t) - n(t) + 1/4
  !+-------------------------------------------------------------------+
  function neq_measure_epot(g,self,params) result(epot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: epot,docc,nt
    integer                             :: i,k,j,N,L
    N = params%Itime
    L = params%Ntau
    !
    if(N==1)then
       nt   = neq_measure_dens(g,self,params)
       docc = neq_measure_docc(g,self,params)
       epot = U0*(docc - nt + 0.25d0)
    else
       nt   = neq_measure_dens(g,self,params)
       docc = neq_measure_docc(g,self,params)
       epot = Uloc*(docc - nt + 0.25d0)
    endif
  end function neq_measure_epot



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function neq_measure_etot(g,self,params) result(etot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: etot,ekin,epot
    ekin = neq_measure_ekin(g,self,params)
    epot = neq_measure_epot(g,self,params)
    etot = ekin + epot
  end function neq_measure_etot


END MODULE NEQ_MEASURE
