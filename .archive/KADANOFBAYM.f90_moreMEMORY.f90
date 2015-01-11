!########################################################################
!PROGRAM  : KADANOFF_BAYM
!TYPE     : module
!PURPOSE  : Evolve the Green's functions G^>,< on the real time axis
!given the initial condition at t=t'=0.
!AUTHORS  : Adriano Amaricci
!########################################################################
module KADANOFBAYM
  !LIBRARY:
  !===================
  USE TOOLS
  USE LATTICE
  USE PROGRESS
  USE FFTW
  USE GRIDS
  USE DLPLOT
  !LOCAL:
  !===================
  USE VARS_GLOBAL
  USE FUNX_NEQ
  implicit none
  private
  public kadanof_baym2Gloc
  save

  !These are in common to the different routines in the module
  complex(8),allocatable,dimension(:,:,:) :: Gkless,Gkgtr
  complex(8),allocatable,dimension(:,:)   :: Ikless,Ikgtr
  complex(8),allocatable,dimension(:,:)   :: Ikless0,Ikgtr0
  complex(8),allocatable,dimension(:)     :: Ikdiag
  complex(8),allocatable,dimension(:,:)   :: Udelta(:,:),Vdelta(:,:)
contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : KADANOF_BAYM2GLOC
  !TYPE     : Subroutine
  !PURPOSE  : Build the local Green's functions G^<,> [G^\l/rceil]
  !along the real time axis, using a discretized verion of the Kadanoff-Baym 
  !equations. The evolution is performed by a two-pass procedure, following:
  ! H.S. K\"ohler, N.H.Kwong, Hashim A. Yousif
  !"A Fortran code for solving the Kadanoff-Baym equations for a homogeneous 
  !fermion system" Computer Physics Communications,123(1999),123-142
  !+-------------------------------------------------------------------+
  subroutine kadanof_baym2Gloc
    integer    :: istep,i,j,ik,jk,it,itau,ix,iy
    logical    :: last
    do i=1,10
       call dump("")
    enddo
    call dump("Entering Kadanoff-Baym")
    call dump("----------------------")

    !Allocate functions:
    !===============================================================
    call alloc_kfunx('a')

    !First loop ever, build the Initial Condition
    !===============================================================
    if(iloop==1)call build_ic

    !Set to Zero loc & k-dep GF:
    !===============================================================
    Gkless=zero;   Gkgtr=zero
    locGless=zero; locGgtr=zero

    !Set the evolution operators
    !===============================================================
    forall(ik=1:Lk,i=0:nstep) 
       Udelta(ik,i)=UdeltaF(ik,i)
       Vdelta(ik,i)=VdeltaF(ik,i)
    end forall

    !Start the time-stepping procedure: solve KB equations
    !===============================================================
    call dump("Solving KB equations:")
    do istep=0,nstep-1
       write(*,'(A,I4,A)',advance="no")"T-step:",istep+1," |"
       last=FF;if(istep==nstep-1)last=TT
       !Recover Initial Condition for Gk^X=<,>,\lceil,\rceil
       call recover_ic()
       !Start growing the functions, 2-pass procedure: 
       call GFstep(1,istep) !1st-pass
       call GFstep(2,istep) !2nd-pass
       !if last t-loop: sum on the fly to get local-GF and local-observables:
       if(last)then
          do ik=1,Lk
             locGless =locGless + Gkless(ik,:,:)*wt(ik)
             locGgtr  =locGgtr  + Gkgtr(ik,:,:)*wt(ik)
          enddo
          call save_initial_condition()
       endif
    enddo
    call dump("")

    !Gloc^>(t',t)= - Gloc^>(t,t')^T
    !===============================================================
    forall(i=0:nstep,j=0:nstep,i>j)locGless(i,j)=-conjg(locGless(j,i))
    forall(i=0:nstep,j=0:nstep,i<j)locGgtr(i,j)=-conjg(locGgtr(j,i))

    !Plot local Green's functions
    !===============================================================
    call plot_dislin3D("locGless_t1t2","X/$\Delta t$","Y/$\Delta t$","Z",&
         t(0:nstep)/dt,t(0:nstep)/dt,locGless)
    call plot_dislin3D("locGgtr_t1t2","X/$\Delta t$","Y/$\Delta t$","Z",&
         t(0:nstep)/dt,t(0:nstep)/dt,locGgtr)

    !Get the (local) current and the (local) occupation
    !===============================================================
    call get_Observables

    !Deallocate to free space:
    !===============================================================
    call alloc_kfunx('d')
    do i=1,10
       call dump("")
    enddo
    return
  end subroutine kadanof_baym2gloc
  !******************************************************************
  !******************************************************************
  !******************************************************************















  !+-------------------------------------------------------------------+
  !PROGRAM  : GFSTEP
  !TYPE     : Subroutine
  !PURPOSE  : Execute the 2pass procedure to solve KB equations:
  !+-------------------------------------------------------------------+
  subroutine GFstep(ips,istep)
    integer :: ips,istep,ik,jk
    integer :: i,j,itau,jtau
    integer :: it
    real(8) :: nk
    complex(8),allocatable,dimension(:,:) :: dummy



    if(ips == 1)then
       !First Pass: get collision integrals up to t=T=istep
       !===============================================================
       Ikless0 = zero;Ikgtr0  = zero;Ikdiag  = zero
       call get_Ikcollision(istep)
       Ikless0  = Ikless
       Ikgtr0   = Ikgtr
       Ikdiag(:)= Ikgtr(:,istep)-Ikless(:,istep)

    elseif(ips==2)then
       !Second Pass: get collision integrals up to t=T+\Delta=istep+1
       !===============================================================
       call get_Ikcollision(istep+1)
       Ikless0  = (Ikless  + Ikless0)/2.d0
       Ikgtr0   = (Ikgtr   + Ikgtr0)/2.d0
       Ikdiag(:)= (Ikgtr(:,istep+1)-Ikless(:,istep+1) + Ikdiag(:))/2.d0
    endif


    !Evolve the solution of KB equations for all the k-points:
    !===============================================================
    do it=0,istep
       Gkless(:,it,istep+1) = Gkless(:,it,istep)*Udelta(:,istep)+Ikless0(:,it)*Vdelta(:,istep)
       Gkgtr(:,istep+1,it)  = Gkgtr(:,istep,it)*conjg(Udelta(:,istep))+Ikgtr0(:,it)*conjg(Vdelta(:,istep))
    enddo
    Gkgtr(:,istep+1,istep)=(Gkless(:,istep,istep)-xi)*conjg(Udelta(:,istep))+Ikgtr0(:,istep)*conjg(Vdelta(:,istep))
    Gkless(:,istep+1,istep+1)=Gkless(:,istep,istep)-xi*dt*Ikdiag(:)
    Gkgtr(:,istep+1,istep+1) = Gkless(:,istep+1,istep+1)-xi
    !forall(i=0:nstep,j=0:nstep,i>j)Gkless(:,i,j)=-conjg(Gkless(:,j,i))
    !forall(i=0:nstep,j=0:nstep,i<j)Gkgtr(:,i,j)=-conjg(Gkgtr(:,j,i))
    !$OMP PARALLEL PRIVATE(i,j)
    !$OMP DO
    do i=0,nstep
       do j=0,nstep
          if(i>j)then
             Gkless(:,i,j)=-conjg(Gkless(:,j,i))
          endif
          if(i<j)then
             Gkgtr(:,i,j)=-conjg(Gkgtr(:,j,i)) 
          endif
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    nk=0.d0
    nk=-2.d0*xi*sum(Gkless(:,istep+1,istep+1)*wt(:))
    if(ips==2)write(*,"(A,f12.9,f12.9)")'Ikdiag,nk:',maxval(abs(real(Ikdiag))),nk!maxval(abs(aimag(Ikdiag)))
    return
  end subroutine GFstep
  !******************************************************************
  !******************************************************************
  !******************************************************************









  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENTS : VER3.0_Apr2010. This version is based on the homogeneity
  ! with repects to k' of the matrices F_{k,k'}(t,t'), related to the 
  ! constant expression of Vhyb{k,k'}\== Vpd.
  ! A more general form has been tested in the past, using direct
  ! expression of the matrices in {k,k'}. See repo.
  !+-------------------------------------------------------------------+
  subroutine get_Ikcollision(Nt)
    integer,intent(in) :: Nt
    integer :: i,j,ik,jk,it,itp,itau,jtau,mk,hk,u
    real(4) :: t0,t1
    complex(8),dimension(Lk)   :: I1,I2
    complex(8),dimension(0:Nt) :: Vless,Vadv,Vret,Vgtr

    Ikless=zero; Ikgtr=zero
    if(Nt==0)return

    forall(i=0:Nt)
       Vless(i)=Sless(i,Nt)+S0less(i-Nt)
       Vadv(i)=conjg(Sret(Nt,i)+S0ret(Nt-i))
       Vret(i)=Sret(Nt,i)+S0ret(Nt-i)
       Vgtr(i)=Sgtr(Nt,i)+S0gtr(Nt-i)
    end forall

    !Get I_2^<(it,Nt) \forall it:0,...,Nt == t=0,...,T+\Delta
    !I1 = \int_0^it G^R_{k,h} * Sigma^<_{h,k'}
    !I2 = \int_0^Nt G^<_{k,h} * Sigma^Adv_{h,k'}
    !===============================================================
    do it=0,Nt
       I1=zero;I2=zero;
       do i=0,it
          I1 = I1 + GkretF(it,i)*Vless(i)
       enddo
       do i=0,Nt
          I2 = I2 + Gkless(:,it,i)*Vadv(i)
       enddo
       Ikless(:,it)=(I1+I2)*dt
    enddo

    !Get I_1^>(Nt,itp) itp:0,...,Nt == t' = 0,...,T+\Delta=Nt
    !I1 = \int_0^t=Nt S^Ret_{k,h} * G^>_{h,k'}
    !I2 = \int_0^t'=itp S^>_{k,h} * G^A_{h,k'}
    !===============================================================
    do itp=0,Nt
       I1=zero;I2=zero;
       do i=0,Nt
          I1 = I1 + Vret(i)*Gkgtr(:,i,itp)
       enddo
       do i=0,itp
          I2 = I2 + Vgtr(i)*conjg(GkretF(itp,i))
       enddo
       Ikgtr(:,itp)=(I1+I2)*dt
    enddo
    return
  contains
    function GkretF(i,j)      
      integer,intent(in)       :: i,j
      complex(8),dimension(Lk) :: GkretF
      GkretF(:) = heaviside(t(i)-t(j))*(Gkgtr(:,i,j)-Gkless(:,i,j))
    end function GkretF
  end subroutine get_Ikcollision
  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  : GET_OBSERVABLES
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine get_observables
    integer                                   :: i,j,it,ik,ix,iy,jx,jy,jk
    type(vect2D),dimension(0:nstep)           :: Jloc
    type(vect2D),dimension(0:nstep,0:Nx,0:Ny) :: Jkvec
    complex(8),dimension(-nstep:nstep)        :: gtless,gtgtr,gtret
    complex(8),dimension(2*nstep)             :: gfret,gfret_tmp
    real(8),dimension(0:nstep)                :: nt
    real(8),dimension(0:nstep,Lk)             :: nk,npi,sorted_nk,sorted_npi
    real(8),dimension(0:Nx,0:Ny,0:nstep)      :: nDens
    character(len=3)                          :: char
    !real(8),dimension(0:nstep,Lk)             :: sorted_epsipi

    !Get the current Jloc(t)=\sum_\ka J_\ka(t) = -e n_\ka(t)*v_\ka(t)
    !===============================================================
    Jloc=Vzero    
    nk=0.d0
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=0,nstep
          nk(i,ik)=-xi*Gkless(ik,i,i)
          Jkvec(i,ix,iy)=nk(i,ik)*velocity(kgrid(ix,iy) + t(i)*Ek)
          Jloc(i) = Jloc(i)  + wt(ik)*Jkvec(i,ix,iy) !nk(i,ik))*velocity( kgrid(ix,iy) + t(i)*Ek )
       enddo
    enddo

    !Sorting n(\ka,t)
    !===============================================================
    forall(i=0:nstep,j=1:Lk)sorted_nk(i,j)=nk(i,sorted_ik(j))
    call shift_kpoint( nk(0:nstep,1:Lk), npi(0:nstep,1:Lk))
    call shift_kpoint( sorted_nk(0:nstep,1:Lk), sorted_npi(0:nstep,1:Lk))
    !call shift_kpoint( transpose(sorted_epsikt(1:Lk,0:nstep)), sorted_epsipi(0:nstep,1:Lk))


    !Fermi surface
    !===============================================================
    forall(i=0:nstep,ik=1:Lk)nDens(ik2ix(ik),ik2iy(ik),i)=sorted_nk(i,ik)
    call plot_dislin3D("nVSk_tzero","$k_x$","$k_y$","$n(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0))
    call plot3Dintensity_movie("FSVSkVSt","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:nstep))
    forall(i=0:nstep,ik=1:Lk)nDens(ik2ix(ik),ik2iy(ik),i)=npi(i,ik)
    call plot3Dintensity_movie("FSVSpiVSt","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:nstep))

    !Plot currents:
    !===============================================================
    if(Efield/=0.d0)then
       call plot("JlocVStime",t(0:nstep),Jloc(0:nstep)%x+Jloc(0:nstep)%y,"",Xlabel="t",Ylabel="Jloc",wlp="wlp")
       if(Lk<2000)call plot_dislinVF("JfieldVSkVSt",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,t(0:nstep),Jkvec(0:nstep,:,:)%x,Jkvec(0:nstep,:,:)%y)
    endif


    !Plot the n_loc:
    !===============================================================
    forall(i=0:nstep)nt(i)=-xi*locGless(i,i)
    call plot("nVStime",t(0:nstep),2.d0*nt(0:nstep),Xlabel="t",Ylabel="n_loc",wlp="wlp")


    !Plot nkVSepsik --> Momenta distributions:
    call splot("nkVStimeVSepsk.ipt",t(0:nstep),sorted_epsik(1:Lk),sorted_nk(0:nstep,1:Lk))

    !Fermi distr. on the CoM reference (moving with Efield):
    !===============================================================
    call plot_dislin2Dmovie("nkVSepskVSt",sorted_epsik,sorted_nk(0:nstep,:),Xlabel="$\Huge\epsilon(k)$",Ylabel="$\Huge n_{k}(t)$",wlp="wlp")
    call plot_dislin2Dmovie("npiVSepskVSt",sorted_epsik,sorted_npi(0:nstep,:),Xlabel="$\Huge\epsilon(k)$",Ylabel="$\Huge n_{\pi}(t)$",wlp="wlp")    


    !Wigner transformed functions (2b completed):
    !===============================================================
    do i=0,nstep
       do j=0,nstep
          gtless(i-j)=locGless(i,j)
          gtgtr(i-j) =locGgtr(i,j)
       enddo
    enddo
    do i=-nstep,nstep
       gtret(i)=heaviside(t(i))*(gtgtr(i)-gtless(i))
    enddo
    call cfft_rt2rw(gtret,gfret,nstep) ;    gfret=gfret*dt ; call swap_fftrt2rw(gfret)
    call splot("locGless_t",t(0:nstep),aimag(gtless(0:nstep)),real(gtless(0:nstep)))
    call splot("locGgtr_t",t(0:nstep),aimag(gtgtr(0:nstep)),real(gtgtr(0:nstep)))
    call splot("locGret_t",t(0:nstep),aimag(gtret(0:nstep)),real(gtret(0:nstep)))
    call splot("locGret_realw",wrmini,aimag(gfret),real(gfret))
  end subroutine get_observables
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine shift_kpoint(arrayIN,arrayOUT)
    integer                 :: i,j,ik,ix,iy,jk,jx,jy
    real(8),dimension(0:,:) :: arrayIN,arrayOUT
    real(8),dimension(2)    :: pi_in
    integer,dimension(2)    :: pi_kcoord

    arrayOUT=0.d0
    do i=0,nstep
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)
          !find the Xcoord of shifted point:
          pi_in(1)=kgrid(ix,iy)%x + (-t(i))*Ek%x
          do j=1,100000
             if(pi_in(1) > pi) then
                pi_in(1)=pi_in(1) - pi2
             elseif(pi_in(1) < -pi) then
                pi_in(1)=pi_in(1) + pi2
             else
                exit
             endif
          enddo

          !find the Ycoord of shifted point:
          pi_in(2)=kgrid(ix,iy)%y + (-t(i))*Ek%y
          do j=1,100000
             if(pi_in(2) > pi) then
                pi_in(2)=pi_in(2) - pi2
             elseif(pi_in(2) < -pi) then
                pi_in(2)=pi_in(2) + pi2
             else
                exit
             endif
          enddo

          !FIND the kgrid point corresponding to Pi-point
          call find2Dmesh(kgrid(0:Nx,1)%x,kgrid(1,0:Ny)%y,pi_in,pi_kcoord)
          jx=pi_kcoord(1)-1 ; jy=pi_kcoord(2)-1
          if(jx < 0  .or. jx > Nx)print*,"error jx=",jx
          if(jy < 0  .or. jy > Ny)print*,"error jy=",jy
          jk=kindex(jx,jy)
          arrayOUT(i,ik)=arrayIN(i,jk)
       enddo
    enddo
  end subroutine shift_kpoint




  !#####################################################################
  !INITIAL CONDITIONS:
  !#####################################################################
  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine build_ic
    integer,parameter :: M=1024
    integer :: i,j,ik,im,itau,jtau,jk
    real(8) :: w,np,en,nloc,nk
    complex(8) :: iw
    complex(8),dimension(2*M) :: fg,Sm,eqSiw
    real(8),dimension(0:M) :: tmpGtau
    call dump("Building initial conditions:")
    nloc=0.d0
    icGkless=zero
    call getGmats(wr,eqSw,eqSiw,beta0)
    Sm=zero
    do i=1,2*M
       w=pi/beta0*dble(2*i-1)
       iw=cmplx(0.d0,w)
       do im=1,Lmu
          Sm(i)=Sm(i)+Vpd**2/(iw+xmu0-epsimu(im))/dble(Lmu)
       enddo
    enddo

    do ik=1,Lk
       call print_bar(ik,Lk)
       fg=zero
       do i=1,2*M
          w=pi/beta0*dble(2*i-1)
          iw=cmplx(0.d0,w)
          fg(i)=one/(iw+xmu0 - epsik(ik)  - eqSiw(i)- Sm(i))
       enddo
       call cfft_iw2it(fg,tmpGtau,beta0) !\tau\in[-beta,0
       nk=-tmpGtau(M) !from G(-tau) = -G(beta-tau)==> G(tau=0-)=-G(beta)
       !G_k^<(0,0)=i*G_k^M(\tau=0-)
       icGkless(ik)=xi*nk 
       call splot("nkVSepsk_init.ipt",epsik(ik),nk,append=TT)!write(200,*)epsik(ik),nk
       nloc=nloc+wt(ik)*fermi(epsik(ik),beta0)
    enddo
    call close_bar
    write(*,"(A,f9.6)")"nloc=",nloc
    return
  end subroutine build_ic
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine recover_ic()
    integer :: i,j,ik,jk,itau
    forall(ik=1:Lk)
       Gkless(ik,0,0)=icGkless(ik)
       Gkgtr(ik,0,0) =icGkless(ik)-xi
    end forall
    return
  end subroutine recover_ic
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine save_initial_condition()
    integer :: ik,jk
    integer :: itau
    forall(ik=1:Lk)icGkless(ik)=Gkless(ik,0,0)
  end subroutine save_initial_condition
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !#####################################################################
  !ALLOCATION
  !#####################################################################
  !+-------------------------------------------------------------------+
  !PROGRAM  : ALLOC_KFUNX
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine alloc_kfunx(char)
    character(len=1) :: char
    if(char=='a')then
       call dump("")
       call dump("Allocating:")
       allocate(Gkless(Lk,0:nstep,0:nstep),Gkgtr(Lk,0:nstep,0:nstep))
       allocate(Ikless(Lk,0:nstep),Ikgtr(Lk,0:nstep),Ikdiag(Lk))
       allocate(Ikless0(Lk,0:nstep),Ikgtr0(Lk,0:nstep))
       allocate(Udelta(Lk,0:nstep),Vdelta(Lk,0:nstep))
       call dump("done")
       call dump("")
    elseif(char=='d')then
       call dump("")
       call dump("Deallocating:")
       deallocate(Gkless,Gkgtr)      
       deallocate(Ikless,Ikgtr,Ikdiag)
       deallocate(Ikless0,Ikgtr0)
       deallocate(Udelta,Vdelta)
       call dump("done")
       call dump("")
    endif
    return
  end subroutine alloc_kfunx
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !#####################################################################
  !EVOLUTION OPERATORS U_\Delta, V_\Delta
  !#####################################################################
  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  elemental function Hbar(ik,istep)
    integer,intent(in) :: ik,istep  
    integer      :: i,j
    complex(8)   :: Hbar
    real(8)      :: tbar
    type(vect2D) :: kt
    !if(Efield/=0.d0)write(*,*)"Remember to take t-t` in Hbar routine!!!"
    !tbar=t(istep) + dt/2.d0
    !i=ik2ix(ik)
    !j=ik2iy(ik)
    !kt=kgrid(i,j) + tbar*Ek
    Hbar=epsikt(ik,istep)!epsk(kt)
  end function Hbar
  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  elemental function UdeltaF(ik,istep) 
    integer,intent(in)    :: ik,istep
    complex(8) :: UdeltaF
    complex(8) :: arg
    arg=Hbar(ik,istep)
    UdeltaF=exp(xi*arg*dt)
  end function UdeltaF
  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  elemental function VdeltaF(ik,istep)
    integer,intent(in)    :: ik,istep
    complex(8) :: VdeltaF
    complex(8) :: arg
    arg=Hbar(ik,istep)
    VdeltaF=exp(xi*arg*dt)
    if(abs(arg*dt) <= 1.d-5)then
       VdeltaF=xi*dt
    else
       VdeltaF=(VdeltaF-1.d0)/arg
    endif
  end function VdeltaF
  !******************************************************************
  !******************************************************************
  !******************************************************************
end module KADANOFBAYM
