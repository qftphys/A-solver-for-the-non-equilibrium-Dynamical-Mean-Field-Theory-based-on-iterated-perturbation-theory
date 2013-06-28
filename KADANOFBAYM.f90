!########################################################################
!PROGRAM  : KADANOFFBAYM
!TYPE     : module
!PURPOSE  : Evolve the Green's functions G^>,< on the real time axis
!given the initial condition at t=t'=0.
!AUTHORS  : Adriano Amaricci
!########################################################################
module KADANOFBAYM
  !LOCAL:
  USE NEQ_VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE INTEGRATE
  implicit none
  private

  !k-dependent GF:
  complex(8),allocatable,dimension(:,:)   :: Gkless,Gkgtr

  !Vector for KB propagation solution
  complex(8),allocatable,dimension(:)     :: Ikless,Ikgtr
  complex(8),allocatable,dimension(:)     :: Ikless0,Ikgtr0
  real(8)                                 :: Ikdiag
  !Auxiliary operators:
  complex(8),allocatable,dimension(:,:)   :: Udelta,Vdelta
  !Chi local:
  real(8),allocatable,dimension(:,:,:,:)  :: chik
  !
  public                                  :: neq_get_localgf_pc

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the k-sum to get the local Green's functions (R,>,<)
  !+-------------------------------------------------------------------+
  subroutine neq_get_localgf_pc()
    call kadanoff_baym_to_localgf()
    call print_out_Gloc()
  end subroutine neq_get_localgf_pc


  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the local Green's functions G^<,>
  !along the real time axis, using a discretized verion of the Kadanoff-Baym 
  !equations. The evolution is performed by a two-pass procedure, following:
  ! H.S. K\"ohler, N.H.Kwong, Hashim A. Yousif
  !"A Fortran code for solving the Kadanoff-Baym equations for a homogeneous 
  !fermion system" Computer Physics Communications,123(1999),123-142
  !+-------------------------------------------------------------------+
  subroutine kadanoff_baym_to_localgf()
    integer                  :: istep,i,j,ik,k

    call allocate_funx

    call buildUV

    call msg("Entering Kadanoff-Baym / Predictor-Corrector")

    !Set to Zero GF and nk:
    locG   = zero
    nk     = 0.d0

    !=============START K-POINTS LOOP======================
    call start_timer
    do ik=1,Lk,mpiSIZE
       Gkless=zero
       Gkgtr =zero

       !t-step loop
       do istep=1,nstep-1
          call step_keldysh_contour_gf(ik,istep)
       enddo

       !sum over k-point
       !sum over k-point
       do i=1,nstep
          do j=1,nstep
             k=pack_index_tri(i,j)
             locG%less(k) = locG%less(k) + Gkless(i,j)*wt(ik)
             locG%gtr(k)  = locG%gtr(k)  + Gkgtr(i,j)*wt(ik)
          enddo
       enddo
       do istep=1,nstep
          nk(istep,ik)=-xi*Gkless(istep,istep)
       enddo
       ! if(fchi)then
       !    call get_chi(ik)
       !    chi = chi + chik*wt(ik)
       ! endif
       call eta(ik,Lk)
    enddo
    call stop_timer
    !=============END K-POINTS LOOP======================


    !Deallocate tmp arrays:
    call deallocate_funx
  end subroutine kadanoff_baym_to_localgf



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+----------------------------------------------------------------+
  !PURPOSE  :  allocation of working array
  !+----------------------------------------------------------------+
  subroutine allocate_funx()
    call msg("Allocating KB memory:")
    !Predictor-corrector solver arrays: store the time-step
    allocate(Ikless(nstep),Ikgtr(nstep))
    allocate(Ikless0(nstep),Ikgtr0(nstep))
    !Aux. operators
    allocate(Udelta(Lk,nstep),Vdelta(Lk,nstep))
    if(fchi)allocate(chik(2,2,nstep,nstep))
    !Allocate k-dependent GF:
    allocate(Gkless(Nstep,Nstep),Gkgtr(Nstep,Nstep))
  end subroutine allocate_funx

  subroutine deallocate_funx()
    call msg("Deallocating KB memory:")
    deallocate(Ikless,Ikgtr)
    deallocate(Ikless0,Ikgtr0)
    deallocate(Udelta,Vdelta)
    if(fchi)deallocate(chik)
    deallocate(Gkless,Gkgtr)
  end subroutine deallocate_funx


  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : Execute the 2pass procedure to solve KB equations:
  !+-------------------------------------------------------------------+
  subroutine step_keldysh_contour_gf(ik,istep)
    integer :: ips,istep,ik
    integer :: i,j,itau
    integer :: it
    !istep==0: provide initial conditions to the real-time part of the contour:
    if(istep==1)then
       Gkless(1,1)=xi*eq_nk(ik)
       Gkgtr(1,1) =xi*(eq_nk(ik)-1.d0)
    endif

    do ips=1,2
       select case(ips)
       case(1)
          !First step: get collision integrals up to t=T=istep
          call get_Ikcollision(istep)
          Ikless0  = Ikless ; Ikgtr0   = Ikgtr
          !Ikdiag   = -2.d0*dreal(Ikless(istep))
          Ikdiag   = dreal(Ikgtr(istep))-dreal(Ikless(istep))
       case(2)
          !Second Pass: get collision integrals up to t=T+\Delta=istep+1
          call get_Ikcollision(istep+1)
          Ikless   = (Ikless  + Ikless0)/2.d0
          Ikgtr    = (Ikgtr   + Ikgtr0)/2.d0
          !Ikdiag   = (-2.d0*dreal(Ikless(istep+1)) + Ikdiag)/2.d0
          Ikdiag   = (dreal(Ikgtr(istep+1))-dreal(Ikless(istep+1)) + Ikdiag)/2.d0 
       end select

       !Evolve the solution of KB equations for all the k-points:
       forall(it=1:istep)
          Gkless(it,istep+1) = Gkless(it,istep)*conjg(Udelta(ik,istep))-Ikless(it)*conjg(Vdelta(ik,istep))
          Gkgtr(istep+1,it)  = Gkgtr(istep,it)*Udelta(ik,istep)-Ikgtr(it)*Vdelta(ik,istep)
       end forall
       Gkgtr(istep+1,istep)=(Gkless(istep,istep)-xi)*Udelta(ik,istep)-Ikgtr(istep)*Vdelta(ik,istep)
       Gkless(istep+1,istep+1)= Gkless(istep,istep)-xi*dt*Ikdiag
       Gkgtr(istep+1,istep+1) = Gkless(istep+1,istep+1)-xi
       do i=1,istep+1
          Gkless(istep+1,i)=-conjg(Gkless(i,istep+1))
          Gkgtr(i,istep+1) =-conjg(Gkgtr(istep+1,i))
       enddo
    enddo
  end subroutine step_keldysh_contour_gf



  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !COMMENTS : VER3.0_Apr2010. This version is based on the homogeneity
  ! with repects to k' of the matrices F_{k,k'}(t,t'), related to the 
  ! constant expression of Vhyb{k,k'}\== Vpd.
  ! A more general form has been tested in the past, using direct
  ! expression of the matrices in {k,k'}. See repo.
  !+-------------------------------------------------------------------+
  subroutine get_Ikcollision(Nt)
    integer,intent(in)          :: Nt
    integer                     :: i,itau,it,itp
    complex(8)                  :: I1,I2,Ib
    complex(8),dimension(0:Nt)  :: Vadv,Vret,Vless,Vgtr
    complex(8),dimension(0:Nt)  :: Fadv,Fret

    Ikless=zero; Ikgtr=zero
    if(Nt==1)return

    do i=1,Nt
       Vless(i)= pack_less(i,Nt,Nstep,Sigma) + pack_less(i,Nt,Nstep,S0)!Sigma%less(i,Nt)    + S0%less(i,Nt)
       Vgtr(i) = pack_gtr(i,Nt,Nstep,Sigma)  + pack_gtr(i,Nt,Nstep,S0)!Sigma%gtr(Nt,i)     + S0%gtr(Nt,i)
       Vret(i) = SretF(Nt,i)                 + S0retF(Nt,i)
       Vadv(i) = conjg(Vret(i))
    end do

    select case(int_method)
    case default             !trapz
       !Get I^<(t=it,t'=Nt) it:0,...,Nt == t=0,...,T+\Delta
       !I1 = \int_0^{t=it}  G^R*Sigma^<
       !I2 = \int_0^{t'=Nt} G^<*Sigma^A
       do it=1,Nt
          do i=1,it
             Fret(i) = GkretF(it,i)
          enddo
          I1=trapz(dt,Fret(1:it)*Vless(1:it))
          I2=trapz(dt,Gkless(it,1:Nt)*Vadv(1:Nt))
          Ikless(it)=I1 + I2
       enddo
       !
       !Get I^>(t=Nt,t'=itp) itp:0,...,Nt == t' = 0,...,T+\Delta=Nt
       !I1 = \int_0^{t=Nt}   S^R*G^>
       !I2 = \int_0^{t'=itp} S^>*G^A
       do itp=1,Nt
          do i=1,itp
             Fadv(i) = conjg(GkretF(itp,i))
          enddo
          I1=trapz(dt,Vret(1:Nt)*Gkgtr(1:Nt,itp))
          I2=trapz(dt,Vgtr(1:itp)*Fadv(1:itp))
          Ikgtr(itp)=I1 + I2
       enddo



    case ("simps")
       do it=1,Nt
          do i=1,it
             Fret(i) = GkretF(it,i)
          enddo
          I1=simps(dt,Fret(1:it)*Vless(1:it))
          I2=simps(dt,Gkless(it,1:Nt)*Vadv(1:Nt))
          Ikless(it)=I1 + I2
       enddo
       !
       do itp=1,Nt
          do i=1,itp
             Fadv(i) = conjg(GkretF(itp,i))
          enddo
          I1=simps(dt,Vret(1:Nt)*Gkgtr(1:Nt,itp))
          I2=simps(dt,Vgtr(1:itp)*Fadv(1:itp))
          Ikgtr(itp)=I1 + I2
       enddo

    case ("rect")
       do it=1,Nt
          do i=1,it
             Fret(i) = GkretF(it,i)
          enddo
          I1=sum(Fret(1:it)*Vless(1:it))*dt
          I2=sum(Gkless(it,1:Nt)*Vadv(1:Nt))*dt
          Ikless(it)=I1 + I2
       enddo
       !
       do itp=1,Nt
          do i=1,itp
             Fadv(i) = conjg(GkretF(itp,i))
          enddo
          I1=sum(Vret(1:Nt)*Gkgtr(1:Nt,itp))*dt
          I2=sum(Vgtr(1:itp)*Fadv(1:itp))*dt
          Ikgtr(itp)=I1 + I2
       enddo


    end select

  end subroutine get_Ikcollision



  !******************************************************************
  !******************************************************************
  !******************************************************************


  function GkretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: GkretF
    GkretF = heaviside(t(i)-t(j))*(Gkgtr(i,j)-Gkless(i,j))
  end function GkretF

  function S0retF(i,j)
    integer,intent(in) :: i,j
    complex(8)         :: S0retF
    S0retF = heaviside(t(i)-t(j))*(pack_gtr(i,j,Nstep,S0)-pack_less(i,j,Nstep,S0))
  end function S0retF

  function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i)-t(j))*(pack_gtr(i,j,Nstep,Sigma)-pack_less(i,j,Nstep,Sigma))
  end function SretF


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine buildUV
    integer :: ik,i
    do  ik=1,Lk
       do i=1,nstep
          Udelta(ik,i)=UdeltaF(ik,i)
          Vdelta(ik,i)=VdeltaF(ik,i)
       enddo
    enddo
  end subroutine buildUV

  function UdeltaF(ik,istep) 
    integer,intent(in)    :: ik,istep
    complex(8) :: UdeltaF
    real(8) :: arg
    arg=Hbar(ik,istep)
    UdeltaF=exp(-xi*arg*dt)
  end function UdeltaF

  function VdeltaF(ik,istep)
    integer,intent(in)    :: ik,istep
    complex(8) :: VdeltaF
    real(8) :: arg
    arg=Hbar(ik,istep)
    VdeltaF=exp(-xi*arg*dt)
    if(abs(arg*dt) <= 1.d-5)then
       VdeltaF=xi*dt
    else
       VdeltaF=(1.d0-VdeltaF)/arg
    endif
  end function VdeltaF



  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PURPOSE  : This is the time-dependent hamiltonian
  ! a more general code would require this quantity to be definied 
  ! externally as an array for every k-point (ik) and time (istep) but
  ! I am lazy...
  !+-------------------------------------------------------------------+
  function Hbar(ik,istep)
    integer,intent(in) :: ik,istep  
    integer      :: i,j
    real(8)      :: Hbar
    real(8)      :: tbar
    type(vect2D) :: kt,Ak
    tbar=t(istep) + dt/2.d0
    i=ik2ix(ik)
    j=ik2iy(ik)
    Ak=Afield(tbar,Ek)
    kt=kgrid(i,j) - Ak
    Hbar=square_lattice_dispersion(kt)
  end function Hbar


  !******************************************************************
  !******************************************************************
  !******************************************************************





  ! subroutine get_chi(ik)    
  !   integer      :: ik
  !   integer      :: i,j,ix,iy
  !   type(vect2D) :: Ak,kt,vel
  !   real(8)      :: sx,sy,ex,ey
  !   !
  !   chik=zero
  !   !
  !   ix=ik2ix(ik)
  !   iy=ik2iy(ik)
  !   do i=1,nstep
  !      Ak = Afield(t(i),Ek)
  !      kt = kgrid(ix,iy)-Ak
  !      !vel= square_lattice_velocity(kt)
  !      sx = 2.d0*sin(kt%x)
  !      sy = 2.d0*sin(kt%y)
  !      ex = 2.d0*cos(kt%x)
  !      ey = 2.d0*cos(kt%y)
  !      do j=1,nstep
  !         chik(1,1,i,j)=-2.d0*sx**2*dimag(GkretF(i,j)*Gkless(j,i))
  !         chik(1,2,i,j)=-2.d0*sx*sy*dimag(GkretF(i,j)*Gkless(j,i))
  !         chik(2,1,i,j)=-2.d0*sx*sy*dimag(GkretF(i,j)*Gkless(j,i))
  !         chik(2,2,i,j)=-2.d0*sy**2*dimag(GkretF(i,j)*Gkless(j,i))
  !      enddo
  !      chik(1,1,i,i)=chik(1,1,i,i) + 2.d0*ex*xi*Gkless(i,i)
  !      chik(2,2,i,i)=chik(2,2,i,i) + 2.d0*ey*xi*Gkless(i,i)
  !   enddo
  ! end subroutine get_chi


  !******************************************************************
  !******************************************************************
  !******************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : evaluate the susceptibility (para-magnetic contribution)
  ! !+-------------------------------------------------------------------+
  ! subroutine get_chi_pm(ik)
  !   integer      :: i,j,ik,ix,iy
  !   type(vect2D) :: Ak,kt,vel
  !   real(8)      :: sx,sy
  !   ix=ik2ix(ik)
  !   iy=ik2iy(ik)
  !   do i=1,nstep
  !      Ak = Afield(t(i),Ek)
  !      kt = kgrid(ix,iy)-Ak
  !      !vel= square_lattice_velocity(kt)
  !      sx = 2.d0*sin(kt%x)
  !      sy = 2.d0*sin(kt%y)
  !      do j=1,nstep
  !         chik_pm(1,1,i,j)=-2.d0*sx**2*wt(ik)*dimag(GkretF(i,j)*Gkless(j,i))
  !         chik_pm(1,2,i,j)=-2.d0*sx*sy*wt(ik)*dimag(GkretF(i,j)*Gkless(j,i))
  !         chik_pm(2,1,i,j)=-2.d0*sx*sy*wt(ik)*dimag(GkretF(i,j)*Gkless(j,i))
  !         chik_pm(2,2,i,j)=-2.d0*sy**2*wt(ik)*dimag(GkretF(i,j)*Gkless(j,i))
  !      enddo
  !   enddo
  !   if(plot3d)call splot3d(trim(plot_dir)//"/Chi_PM",t(0:),t(0:),chi_pm(1,1,0:,0:))
  ! end subroutine get_chi_pm



  !******************************************************************
  !******************************************************************
  !******************************************************************




  !   !+-------------------------------------------------------------------+
  !   !PURPOSE  : evaluate the susceptibility (dia-magnetic contribution)
  !   !+-------------------------------------------------------------------+
  !   subroutine get_chi_dia(ik)
  !     integer       :: i,j,ik,ix,iy
  !     type(vect2D)  :: vel,kt,Ak
  !     real(8)       :: eab(2)
  !     ix=ik2ix(ik)
  !     iy=ik2iy(ik)
  !     do i=1,nstep
  !        Ak = Afield(t(i),Ek)
  !        kt = kgrid(ix,iy)-Ak
  !        eab(1)=2.d0*cos(kt%x)
  !        eab(2)=2.d0*cos(kt%y)
  !        chik_dia(1,1,i)=2.d0*wt(ik)*eab(1)*xi*Gkless(i,i)
  !        chik_dia(2,2,i)=2.d0*wt(ik)*eab(2)*xi*Gkless(i,i)
  !     enddo
  !  enddo
  !  call splot("ChiDIA_t.ipt",t(0:),chi_dia(1,1,0:))
  ! end subroutine get_chi_dia

  !******************************************************************
  !******************************************************************
  !******************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : Solve the equilibrium case
  ! !+-------------------------------------------------------------------+
  ! subroutine get_equilibrium_localgf()
  !   integer    :: i,j,ik
  !   complex(8) :: A,zetan
  !   real(8)    :: w,n
  !   complex(8) :: funcM(L),sigmaM(L)
  !   real(8)    :: funcT(0:L) 
  !   if(mpiID==0)then
  !      !Get Sret(w) = FFT(Sret(t-t'))
  !      forall(i=0:nstep,j=0:nstep) sf%ret%t(i-j)=heaviside(t(i-j))*(Sigma%gtr(i,j)-Sigma%less(i,j))
  !      sf%ret%t=exa*sf%ret%t ; call fftgf_rt2rw(sf%ret%t,sf%ret%w,nstep) ; sf%ret%w=dt*sf%ret%w

  !      !Get locGret(w)
  !      gf%ret%w=zero
  !      do i=1,2*nstep
  !         w=wr(i)
  !         zetan=cmplx(w,eps,8)-sf%ret%w(i) !-eqsbfret(i)
  !         do ik=1,Lk
  !            gf%ret%w(i)=gf%ret%w(i)+wt(ik)/(zetan-epsik(ik))
  !         enddo
  !      enddo

  !      !Get locG<(w/t),locG>(w/t)
  !      gf%less%w=less_component_w(gf%ret%w,wr,beta)
  !      gf%gtr%w=gtr_component_w(gf%ret%w,wr,beta)
  !      call fftgf_rw2rt(gf%less%w,gf%less%t,nstep)  ; gf%less%t=exa*fmesh/pi2*gf%less%t
  !      call fftgf_rw2rt(gf%gtr%w,gf%gtr%t,nstep)    ; gf%gtr%t=exa*fmesh/pi2*gf%gtr%t


  !      forall(i=0:nstep,j=0:nstep)
  !         locG%less(i,j) = gf%less%t(i-j)
  !         locG%gtr(i,j)  = gf%gtr%t(i-j)
  !         gf%ret%t(i-j) = heaviside(t(i-j))*(locG%gtr(i,j)-locG%less(i,j))
  !      end forall

  !      !This is just to get n(k)
  !      call get_matsubara_gf_from_dos(wr,sf%ret%w,sigmaM,beta)
  !      do ik=1,Lk
  !         funcM=zero
  !         do i=1,L
  !            w=pi/beta*dble(2*i-1) ; zetan=cmplx(0.d0,w,8) - sigmaM(i)
  !            funcM(i)=one/(zetan - epsik(ik))
  !         enddo
  !         call fftgf_iw2tau(funcM,funcT,beta)
  !         n=-funcT(L)
  !         nk(:,ik)=n
  !      enddo
  !   endif
  !   call MPI_BCAST(locG%less,(nstep+1)**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  !   call MPI_BCAST(locG%gtr,(nstep+1)**2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  !   call MPI_BCAST(nk,(nstep+1)*Lk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  !   call splot('nkVSepsk.ipt',epsik,nk(nstep/2,:),append=.true.)
  !   call splot('locSM_iw.ipt',wm,sigmaM,append=.true.)
  !   call splot("eqG_w.ipt",wr,gf%ret%w,append=.true.)
  !   call splot("eqSigma_w.ipt",wr,sf%ret%w,append=.true.)
  !   return
  ! end subroutine get_equilibrium_localgf



  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : print out useful information
  !+-------------------------------------------------------------------+
  subroutine print_out_Gloc()
    integer                         :: i,j,ix,iy,ik
    type(vect2D)                    :: Jk,Ak
    type(vect2D),dimension(0:nstep) :: Jloc                   !local Current 
    real(8),dimension(0:nstep)      :: nt,modJloc             !occupation(time)
    !call write_keldysh_contour_gf(locG,trim(data_dir)//"/locG")
    call store_data("nk.neq",nk)

    ! forall(i=0:nstep,j=0:nstep)
    !    gf%less%t(i-j) = locG%less(i,j)
    !    gf%gtr%t(i-j)  = locG%gtr(i,j)
    !    gf%ret%t(i-j)  = heaviside(t(i-j))*(locG%gtr(i,j)-locG%less(i,j))
    ! end forall
    ! !if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
    ! call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ;  gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
    ! call splot("locGless_t.ipt",t,gf%less%t,append=.true.)
    ! call splot("locGgtr_t.ipt",t,gf%gtr%t,append=.true.)
    ! call splot("locGret_t.ipt",t,gf%ret%t,append=.true.)
    ! call splot("locGret_realw.ipt",wr,gf%ret%w,append=.true.)
    ! call splot("locDOS.ipt",wr,-aimag(gf%ret%w)/pi,append=.true.)

    ! if(fchi)then
    !    call store_data(trim(data_dir)//"/locChi_11.data",chi(1,1,0:,0:))
    !    call store_data(trim(data_dir)//"/locChi_12.data",chi(1,2,0:,0:))
    !    call store_data(trim(data_dir)//"/locChi_21.data",chi(2,1,0:,0:))
    !    call store_data(trim(data_dir)//"/locChi_22.data",chi(2,2,0:,0:))
    ! endif

    do i=1,nstep
       nt(i)=-xi*pack_less_tri(i,i,locG)
    enddo
    Jloc=Vzero    
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=1,nstep
          Ak= Afield(t(i),Ek)
          Jk= nk(i,ik)*square_lattice_velocity(kgrid(ix,iy) - Ak)
          Jloc(i) = Jloc(i) +  wt(ik)*Jk
       enddo
    enddo
    call splot("nVStime.neq",t,2.d0*nt,append=.true.)
    if(Efield/=0.d0)call splot("JlocVStime.neq",t,Jloc(:)%x,Jloc(:)%y,append=.true.)
  end subroutine print_out_Gloc



  !******************************************************************
  !******************************************************************
  !******************************************************************

end module KADANOFBAYM
