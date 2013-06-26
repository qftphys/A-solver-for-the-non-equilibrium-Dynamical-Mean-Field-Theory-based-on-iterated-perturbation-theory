MODULE RESULTS
  USE NEQ_VARS_GLOBAL
  USE CONTOUR_GF
  USE ELECTRIC_FIELD
  USE BATH
  implicit none
  private

  !Add here the extra variables 
  integer                                :: i,j,ik,loop,narg,iarg,k,ia,ir,irel,iave
  complex(8),dimension(:,:),allocatable  :: locGret,Sret,G0ret
  logical                                :: file

  public :: plot_results

contains

  !+-------------------------------------------------------------------+
  subroutine plot_results(directory)
    character(len=*),optional :: directory
    character(len=32)         :: DIR

    DIR='./' ; if(present(directory))DIR=reg(directory)

    !call create_data_dir( reg(DIR) )
    call msg("Results are stored in "//reg(DIR) )   

    call evaluate_print_observables( reg(DIR) )

    return
  end subroutine plot_results
  !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  ! subroutine massive_allocation
  ! do i=1,nstep
  !    do j=1,nstep
  !       k = pack_index(i,j,Nstep)
  !       Sret(i,j)      = heaviside(t(i)-t(j))*(Sigma%gtr(i,j) - Sigma%less(i,j))
  !       locGret(i,j)   = heaviside(t(i)-t(j))*(locG%gtr(i,j) - locG%less(i,j))
  !       G0ret(i,j)     = heaviside(t(i)-t(j))*(G0%gtr(i,j)   - G0%less(i,j))
  !    enddo
  ! enddo
  !   ! allocate(trel(-nstep:nstep),tave(0:nstep))
  !   ! allocate(wgnGless(-nstep:nstep,0:nstep),&
  !   !      wgnGgtr(-nstep:nstep,0:nstep),     &
  !   !      wgnGret(-nstep:nstep,0:nstep))
  !   ! allocate(wgnSless(-nstep:nstep,0:nstep),&
  !   !      wgnSgtr(-nstep:nstep,0:nstep),     &
  !   !      wgnSret(-nstep:nstep,0:nstep))
  !   ! allocate(gfret_wgn(0:nstep,-nstep:nstep),gfless_wgn(0:nstep,-nstep:nstep),sfret_wgn(0:nstep,-nstep:nstep))
  !   ! allocate(nf_wgn(0:nstep,-nstep:nstep))
  ! end subroutine massive_allocation
  !+-------------------------------------------------------------------+

  ! !+-------------------------------------------------------------------+
  ! subroutine massive_deallocation
  !   deallocate(locGret,Sret,G0ret)
  !   deallocate(nk)
  !   ! call deallocate_gf(gf0)
  !   ! call deallocate_gf(gf)
  !   ! call deallocate_gf(sf)
  !   ! deallocate(trel,tave)
  !   ! deallocate(wgnGless,wgnGgtr,wgnSless,wgnSgtr,wgnGret,wgnSret)
  !   ! deallocate(gfret_wgn,gfless_wgn,sfret_wgn)
  !   ! deallocate(nf_wgn)
  !   if(fchi)deallocate(chi)
  ! end subroutine massive_deallocation
  ! !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  subroutine evaluate_print_observables(dir)
    character(len=*)                        :: dir
    integer                                 :: i,ik,ix,iy,it,is,step,iw,unit
    complex(8)                              :: I1,Ib
    type(vect2D)                            :: Ak,kt,Jk
    type(vect2D),dimension(nstep)           :: Jloc,Jheat,Xpos !local Current 
    type(vect2D),dimension(nstep,0:Nx,0:Ny) :: Jkvec                  !current vector field
    real(8),dimension(nstep,Lk)             :: npi                    !covariant occupation n(\pi=\ka+\Ekt) 
    real(8),dimension(nstep)                :: nt,Stot   !occupation(time)
    real(8),dimension(0:Nx,0:Ny,nstep)      :: nDens                  !occupation distribution on the k-grid
    real(8),dimension(nstep)                :: Ekin,Epot,Eb,Etot,doble!energies and double occ.
    real(8),allocatable,dimension(:)        :: sorted_epsik
    integer,allocatable,dimension(:)        :: sorted_ik
    real(8),dimension(nstep,Lk)             :: sorted_nk,sorted_npi   !sorted arrays
    real(8),dimension(nstep,Lk)             :: epi,sorted_epi   !sorted arrays
    integer,dimension(:),allocatable        :: reduced_ik             !reduced arrays
    real(8),dimension(:),allocatable        :: reduced_epsik
    real(8),dimension(:,:),allocatable      :: reduced_nk,reduced_npi,reduced_epi
    real(8),dimension(2,2,nstep,nstep)      :: oc,redoc
    !
    real(8)                                 :: ome(nstep)
    complex(8),allocatable                  :: ocw(:)
    real(8)                                 :: oct(-nstep:nstep)
    complex(8)                              :: swcond(2,2,0:nstep,nstep)


    call msg("Print Out Results (may take a while)")

    !SORTING:
    call msg("Sorting")
    allocate(sorted_epsik(Lk),sorted_ik(Lk))
    !
    sorted_epsik=epsik ; call sort_array(sorted_epsik,sorted_ik)
    !
    forall(i=1:nstep,ik=1:Lk)sorted_nk(i,ik) = nk(i,sorted_ik(ik))


    !COVARIANT transformation: \ka --> \pi = \ka + \Ek*t
    call msg("Covariant transformation")
    call shift_kpoint(nk(1:nstep,1:Lk), npi(1:nstep,1:Lk))
    call shift_kpoint(sorted_nk(1:nstep,1:Lk), sorted_npi(1:nstep,1:Lk))


    !REDUCTION of the k-grid:
    call msg("Reducing BZ")
    step=Lk/Lkreduced; if(step==0)step=1
    call square_lattice_reduxGrid_dimension(Lk,step,Lkreduced)
    allocate(reduced_ik(Lkreduced),&
         reduced_epsik(Lkreduced), &
         reduced_epi(1:nstep,Lkreduced),&
         reduced_nk(1:nstep,Lkreduced),reduced_npi(1:nstep,Lkreduced))
    call square_lattice_reduxGrid_index(Lk,step,reduced_ik)
    call square_lattice_reduxGrid_dispersion_array(sorted_epsik,reduced_ik,reduced_epsik)
    forall(ik=1:Lkreduced)reduced_nk(:,ik)  = sorted_nk(:,reduced_ik(ik))
    forall(ik=1:Lkreduced)reduced_npi(:,ik) = sorted_npi(:,reduced_ik(ik))


    !Get the CURRENT Jloc(t)=\sum_\ka J_\ka(t) = -e n_\ka(t)*v_\ka(t)
    call msg("Current Field")
    Jloc=Vzero ; Jheat=Vzero   ;Stot=0.d0
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=1,nstep
          Ak          = Afield(t(i),Ek)
          kt          = kgrid(ix,iy)-Ak
          epi(i,ik)   = square_lattice_dispersion(kt)
          Jk          = nk(i,ik)*square_lattice_velocity(kt)
          Jkvec(i,ix,iy)  = Jk
          Jloc(i)         = Jloc(i) +  wt(ik)*Jk !todo *2.d0 spin degeneracy
          Jheat(i)        = Jheat(i)+  wt(ik)*epi(i,ik)*Jk !todo *2.d0 spin degeneracy
          if(i>1)Stot(i)         = Stot(i) -  2.d0*wt(ik)*(nk(i,ik)*log(nk(i,ik)))!+(1.d0-nk(i,ik))*log(1.d0-nk(i,ik)))
       enddo
    enddo
    Stot(1)=0.d0

    forall(i=1:nstep,ik=1:Lk)sorted_epi(i,ik) = epi(i,sorted_ik(ik))
    forall(ik=1:Lkreduced)reduced_epi(:,ik) = sorted_epi(:,reduced_ik(ik))

    !OCCUPATION :
    call msg("Get occupation")
    do i=1,nstep
       nt(i)=-xi*pack_less_tri(i,i,locG)
    enddo

    !OCCUPATION density:
    forall(i=1:nstep,ik=1:Lk)nDens(ik2ix(ik),ik2iy(ik),i)=npi(i,ik)

    !ENERGY kinetic
    Ekin=zero
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=1,nstep
          Ak=Afield(t(i),Ek)
          Ekin(i) = Ekin(i) +  wt(ik)*square_lattice_dispersion(kgrid(ix,iy) - Ak)*nk(i,ik)
       enddo
    enddo


    !ENERGY potential && total 
    !Get Epot = <V>(it)= xi/2\lim_{t'-->t}\sum_k[xi\partial_t - h_{0,k}(t)] G^<_k(t,t') 
    ! by eq. of motion = xi/2\lim_{t'-->t}\sum_k{\delta^<(t,t')+\int_0^t S^R(t,z)*G^<_k(z,t')+\int_0^t' S^<(t,z)*G^A_k(z,t')}
    !                  = xi/2\lim_{t'-->t}{\delta^<(t,t')+\int_0^t S^R(t,z)*G^<_loc(z,t')+\int_0^t' S^<(t,z)*G^A_loc(z,t')}
    !                  = xi/2{1+\int_0^t [S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*G^A_loc(z,t)]}
    ! cnst disregarded = xi/2\int_0^t {S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*[G^R_loc(t,z)]^+}
    do it=1,nstep
       I1=zero; Ib=zero
       do i=1,it
          I1 = I1 + SretF(it,i)*pack_less_tri(i,it,locG) + pack_less(i,it,Nstep,Sigma)*conjg(GretF(it,i))
       enddo
       Epot(it)= -xi*I1*dt/2.d0
    enddo
    !Get Etot = Ekin + Epot
    Etot=Ekin+Epot

    do i=1,nstep
       Xpos(i)=0.d0
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)          
          do k=0,i
             Jk=Jkvec(k,ix,iy)
             Xpos(i)=Xpos(i)+Jk*wt(ik)*dt
          enddo
       enddo
    enddo

    !Double OCCUPATION:
    doble= 0.5d0*(2.d0*nt) - 0.25d0 ; if(U/=0)doble = Epot/U + 0.5d0*(2.d0*nt)- 0.25d0

    if(fchi)then
       !Get the optical conductivity \sigma(t,t') from the susceptibility \Chi:
       oc=0.d0
       do i=1,nstep
          do j=1,nstep
             do it=j,nstep
                oc(1,1,i,j)=oc(1,1,i,j)-chi(1,1,i,it)*dt
                oc(1,2,i,j)=oc(1,2,i,j)-chi(1,2,i,it)*dt
                oc(2,1,i,j)=oc(2,1,i,j)-chi(2,1,i,it)*dt
                oc(2,2,i,j)=oc(2,2,i,j)-chi(2,2,i,it)*dt
             enddo
          enddo
       enddo
       ! forall(i=1:nstep,j=1:nstep)oct(i-j) = oc(1,1,i,j)
       ! call splot("OC_t.neq",t(0:),oct(0:))
       ! !Partial FT to \sigma(t,w)
       ! ome = wr(nstep+1:2*nstep)
       ! swcond=0.d0
       ! do i=1,nstep
       !    do it=0,nstep
       !       do is=0,it
       !          swcond(1,1,it,i)=swcond(1,1,it,i)+oc(1,1,it,it-is)*exp(-xi*ome(i)*t(is))*dt
       !       enddo
       !    enddo
       ! enddo
       ! do it=0,nstep
       !    call splot("allOC_it_realw.neq",ome,swcond(1,1,it,:),append=.true.)
       ! end do
       ! allocate(ocw(nstep))
       ! ocw=zero
       ! do iw=1,nstep
       !    do it=0,nstep
       !       ocw(iw)=ocw(iw) + oct(it)*exp(xi*ome(iw)*t(it))*fmesh
       !    enddo
       ! enddo
       ! call splot("OC_realw.neq",ome,ocw)
    endif


    !====================================================================================
    !PRINT:
    !====================================================================================
    call msg("Print n(t)")
    unit=free_unit()
    open(unit,file=dir//"columns_info.neq")
    write(unit,"(A1,A19,12A20)")"#","1t","2Jx","3Jy","4n","5docc","6Etot","7Ekin","8Epot","9S","10Jhx","11Jhy","12Rx","13Ry"
    close(unit)
    open(unit,file=dir//"/observables.neq")
    do i=1,nstep
       write(unit,"(13F20.12)")t(i),Jloc(i)%x,Jloc(i)%y,2.d0*nt(i),doble(i),Etot(i),Ekin(i),Epot(i),Stot(i),Jheat(i)%x,Jheat(i)%y,Xpos(i)%x,Xpos(i)%y
    enddo
    close(unit)

    call msg("Print n(k,t)")
    call splot3d(dir//"/nVStimeVSepsk3D.neq",t,reduced_epsik,reduced_nk)
    do i=1,nstep
       call splot(dir//"/nVSepi.neq",reduced_epsik,reduced_npi(i,:),append=.true.)
    enddo

    !Fermi Surface plot:
    call msg("Print FS(k,t)")
    call splot3d(dir//"/3dFSVSpiVSt.neq",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,:))

    if(fchi)then
       call splot3d(dir//"/OC.neq",t,t,oc(1,1,:,:))
       !call splot3d(dir//"/redOC.neq",t,t,redoc(1,1,0:nstep,0:nstep))
    endif

    !Plot local functions:
    !===========================================================================
    if(plot3D)then
       call plot_keldysh_contour_tri_gf(t,G0,Nstep,dir//"/G0")
       call plot_keldysh_contour_tri_gf(t,locG,Nstep,dir//"/locG")
       call plot_keldysh_contour_gf(t,Sigma,dir//"/Sigma")
    endif

  end subroutine evaluate_print_observables
  !+-------------------------------------------------------------------+



  ! !+-------------------------------------------------------------------+
  ! subroutine get_plot_quasiWigner_functions(dir)
  !   character(len=*)  :: dir
  !   integer           :: i,j
  !   real(8),dimension(2*nstep)         :: phi
  !   complex(8),dimension(-nstep:nstep) :: gtkel
  !   complex(8),dimension(2*nstep)      :: gfkel
  !   forall(i=0:nstep,j=0:nstep)
  !      gf%less%t(i-j) = locG%less(i,j)
  !      gf%gtr%t(i-j)  = locG%gtr(i,j)
  !      gf%ret%t(i-j)  = locGret(i,j)
  !   end forall
  !   if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0
  !   call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep) ;    gf%ret%w=gf%ret%w*dt ; call swap_fftrt2rw(gf%ret%w)
  !   call splot(dir//"/locGless_t.neq",t(-nstep:nstep),gf%less%t,append=.true.)
  !   call splot(dir//"/locGgtr_t.neq",t(-nstep:nstep),gf%gtr%t,append=.true.)
  !   call splot(dir//"/locGret_t.neq",t(-nstep:nstep),gf%ret%t,append=.true.)
  !   call splot(dir//"/locGret_realw.neq",wr,gf%ret%w,append=.true.)
  !   call splot(dir//"/locDOS.neq",wr,-aimag(gf%ret%w)/pi)
  !   forall(i=0:nstep,j=0:nstep)
  !      gf0%less%t(i-j)= G0%less(i,j)
  !      gf0%gtr%t(i-j) = G0%gtr(i,j)
  !      gf0%ret%t(i-j) = G0ret(i,j)
  !      sf%less%t(i-j) = Sigma%less(i,j)
  !      sf%gtr%t(i-j)  = Sigma%gtr(i,j)
  !      sf%ret%t(i-j)  = Sret(i,j)
  !   end forall
  !   if(heaviside(0.d0)==1.d0)gf0%ret%t(0)=gf0%ret%t(0)/2.d0 !; gf0%ret%t(0)=-xi
  !   if(heaviside(0.d0)==1.d0)sf%ret%t(0)=sf%ret%t(0)/2.d0
  !   if(heaviside(0.d0)==1.d0)gf%ret%t(0)=gf%ret%t(0)/2.d0   !; gf%ret%t(0)=-xi
  !   call splot(dir//"/G0less_t.neq",t(-nstep:nstep),gf0%less%t)
  !   call splot(dir//"/G0gtr_t.neq",t(-nstep:nstep),gf0%gtr%t)
  !   call splot(dir//"/G0ret_t.neq",t(-nstep:nstep),gf0%ret%t)
  !   call splot(dir//"/Sless_t.neq",t(-nstep:nstep),sf%less%t)
  !   call splot(dir//"/Sgtr_t.neq",t(-nstep:nstep),sf%gtr%t)
  !   call splot(dir//"/Sret_t.neq",t(-nstep:nstep),sf%ret%t)
  !   !Obtain && plot Real frequency Functions:
  !   !===========================================================================
  !   call fftgf_rt2rw(gf0%ret%t,gf0%ret%w,nstep) ;    gf0%ret%w=gf0%ret%w*dt ; call swap_fftrt2rw(gf0%ret%w)
  !   call fftgf_rt2rw(gf%ret%t,gf%ret%w,nstep)   ;    gf%ret%w=gf%ret%w*dt   ; call swap_fftrt2rw(gf%ret%w)
  !   call fftgf_rt2rw(sf%ret%t,sf%ret%w,nstep)   ;    sf%ret%w=dt*sf%ret%w   ; call swap_fftrt2rw(sf%ret%w)
  !   call splot(dir//"/G0ret_realw.neq",wr,gf0%ret%w)
  !   call splot(dir//"/Sret_realw.neq",wr,sf%ret%w)
  !   call splot(dir//"/DOS.neq",wr,-aimag(gf%ret%w)/pi)
  !   forall(i=0:nstep,j=0:nstep)gtkel(i-j) = locG%less(i,j)+locG%gtr(i,j)
  !   call fftgf_rt2rw(gtkel,gfkel,nstep) ; gfkel=gfkel*dt ; call swap_fftrt2rw(gfkel)
  !   call splot(dir//"/locGkel_realw.neq",wr,gfkel)
  !   phi = xi*gfkel/aimag(gf%ret%w)/2.d0
  !   do i=1,2*nstep
  !      if(wr(i)>-4.d0)exit
  !   enddo
  !   call splot(dir//"/phi_realw.neq",wr(i:(2*nstep-i)),phi(i:(2*nstep-i)))
  !   return
  ! end subroutine get_plot_quasiWigner_functions
  ! !+-------------------------------------------------------------------+



  ! !+-------------------------------------------------------------------+
  ! subroutine plot_wigner_functions(dir)
  !   character(len=*)  :: dir
  !   complex(8)        :: delta
  !   call init_trel(t,trel,nstep)
  !   call init_tave(t,tave,nstep)
  !   !Perform the Wigner Rotation:
  !   call msg("Perform Wigner Rotation")
  !   wgnGless= wigner_transform(locG%less,nstep)
  !   wgnGret = wigner_transform(locGret,nstep)
  !   wgnSless= wigner_transform(Sigma%less,nstep)
  !   wgnSret = wigner_transform(Sret,nstep)
  !   call system("if [ ! -d WIGNER ]; then mkdir WIGNER; fi")
  !   ! call plot_3D("wgnGless3D","X","Y","Z",trel(-nstep:nstep)/dt,tave(1:nstep)/dt,wgnGless(-nstep:nstep,1:nstep))
  !   ! call plot_3D("wgnSless3D","X","Y","Z",trel(-nstep:nstep)/dt,tave(1:nstep)/dt,wgnSless(-nstep:nstep,1:nstep))
  !   call system("mv wgn*3D WIGNER/")
  !   delta=(one+xi)/dble(nstep)
  !   do ia=0,nstep
  !      call fftgf_rt2rw(wgnGret(:,ia),gfret_wgn(ia,:),nstep)  ;gfret_wgn(ia,:)=gfret_wgn(ia,:)*dt;call swap_fftrt2rw(gfret_wgn(ia,:))
  !      call fftgf_rt2rw(wgnGless(:,ia),gfless_wgn(ia,:),nstep);gfless_wgn(ia,:)=gfless_wgn(ia,:)*dt;call swap_fftrt2rw(gfless_wgn(ia,:))
  !      call fftgf_rt2rw(wgnSret(:,ia),sfret_wgn(ia,:),nstep)  ;sfret_wgn(ia,:)=sfret_wgn(ia,:)*dt;call swap_fftrt2rw(sfret_wgn(ia,:))
  !      call splot("WIGNER/wgnDOS.neq",wr,(-aimag(gfret_wgn(ia,:) - delta*dble(ia)*pi)/pi),append=.true.)
  !      call splot("WIGNER/wgnSigma_realw.neq",wr,(sfret_wgn(ia,:) + delta*dble(ia)),append=.true.)
  !      call splot("WIGNER/wgnGless_realw.neq",wr,(gfless_wgn(ia,:) + delta*dble(ia)),append=.true.)
  !      ! nf_wgn(ia,:) = -aimag(gfless_wgn(ia,:))!/aimag(gfret_wgn(ia,:))/pi2
  !      nf_wgn(ia,:) = -xi*gfless_wgn(ia,:)/aimag(gfret_wgn(ia,:))
  !      call splot("n_wgnVSepi.neq",wr(:),nf_wgn(ia,:),append=.true.)
  !   enddo
  !   call splot3d("WIGNER/wgndosVSrealwVStime.neq",tave(0:nstep),wr(-nstep:nstep),-aimag(gfret_wgn(0:nstep,-nstep:nstep))/pi)
  !   call splot3d("WIGNER/wgnnfVSrealwVStime.neq",tave(0:nstep),wr(-nstep:nstep),nf_wgn(0:nstep,-nstep:nstep))
  ! end subroutine plot_wigner_functions
  ! !+-------------------------------------------------------------------+



  !+-------------------------------------------------------------------+
  subroutine shift_kpoint(arrayIN,arrayOUT)
    integer                 :: i,j,ik,ix,iy,jk,jx,jy
    real(8),dimension(0:,:) :: arrayIN,arrayOUT
    real(8),dimension(2)    :: pi_in
    integer,dimension(2)    :: pi_kcoord
    type(vect2D)            :: Ak
    arrayOUT=0.d0
    do i=1,Nstep
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)
          !find the Xcoord of shifted point:
          Ak=Afield(t(i),Ek)
          pi_in(1)=kgrid(ix,iy)%x + Ak%x !+ (-t(i))*Ek%x
          do j=1,1000000
             if(pi_in(1) > pi) then
                pi_in(1)=pi_in(1) - pi2
             elseif(pi_in(1) < -pi) then
                pi_in(1)=pi_in(1) + pi2
             else
                exit
             endif
          enddo
          !find the Ycoord of shifted point:
          pi_in(2)=kgrid(ix,iy)%y + Ak%y !+ (-t(i))*Ek%y
          do j=1,1000000
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
  !+-------------------------------------------------------------------+


  ! !+-------------------------------------------------------------------+
  ! function wigner_transform(Gin,M)
  !   integer                       :: M,ir,ia,i,j
  !   complex(8),dimension(M,M) :: Gin
  !   complex(8),dimension(-M+1:M-1,M):: Gout,wigner_transform
  !   do ir=-nstep+1,nstep-1
  !      do ia=1,nstep
  !         forall(j=1:nstep,i=1:nstep, i-j==ir .AND. (i+j)/2==ia)
  !            Gout(ir,ia)=Gin(i,j)
  !         end forall
  !      enddo
  !   enddo
  !   wigner_transform=Gout
  ! end function wigner_transform
  ! !+-------------------------------------------------------------------+


  ! !+-------------------------------------------------------------------+
  ! subroutine init_trel(time,tr,M)
  !   integer                            :: i,j,ir
  !   integer,intent(in)                 :: M
  !   real(8),dimension(-M+1:M-1),intent(in) :: time
  !   real(8),dimension(-M+1:M-1),intent(out):: tr
  !   do ir=-nstep,nstep
  !      forall(j=1:nstep,i=1:nstep, i-j==ir)tr(ir)=time(i)-time(j)
  !   enddo
  !   return
  ! end subroutine init_trel
  ! !+-------------------------------------------------------------------+


  ! !+-------------------------------------------------------------------+
  ! subroutine init_tave(time,ta,M)
  !   integer                            :: i,j,ia
  !   integer,intent(in)                 :: M
  !   real(8),dimension(-M+1:M-1),intent(in) :: time
  !   real(8),dimension(M),intent(out) :: ta
  !   do ia=1,nstep
  !      forall(j=1:nstep,i=1:nstep, (i+j)/2==ia)ta(ia)=(time(i)+time(j))/2.d0
  !   enddo
  !   return
  ! end subroutine init_tave
  ! !+-------------------------------------------------------------------+


  !+-------------------------------------------------------------------+
  function GretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: GretF
    GretF = heaviside(t(i)-t(j))*(pack_gtr_tri(i,j,locG)-pack_less_tri(i,j,locG))
  end function GretF
  !-------------------------------------------------------!


  !-------------------------------------------------------!
  function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i)-t(j))*(pack_gtr(i,j,Nstep,Sigma)-pack_less(i,j,Nstep,Sigma))
  end function SretF
  !-------------------------------------------------------!


end MODULE RESULTS
