program getDATA
  USE VARS_GLOBAL
  USE ELECTRIC_FIELD
  USE FUNX_NEQ
  USE DLPLOT
  implicit none
  integer                                :: i,j,ik,loop,narg,iarg
  character(len=4)                       :: char
  complex(8),dimension(:,:),allocatable  :: locGret,Sret,impGret,G0ret
  logical                                :: file


  call read_input_init()
  DIR="RESULTS"
  FFILE="false.file"

  if(Lkreduced<200)Lkreduced=200
  call Build_2DSquareLattice(Nx,Ny,Lk)
  allocate(epsik(Lk),epsimu(Lmu),sorted_epsik(Lk),sorted_ik(Lk))
  call get_epsik(epsik,ts,0.d0,sorted_epsik,sorted_ik)
  call get_epsimu(emin,epsimu)

  !SET EXTERNAL ELECTRIC FIELD:
  Ek%x=1.d0;Ek%y=1.d0;Ek=(Efield*pi2/alat)*Ek !

  call system("if [ -e "//adjustl(trim(SRC))//".tgz ]; then tar xzvf "//adjustl(trim(SRC))//".tgz ;fi")
  call get_data()
  call system("if [ -e "//adjustl(trim(SRC))//".tgz ]; then rm -rf "//adjustl(trim(SRC))//" ;fi")

contains

  !+-------------------------------------------------------------------+
  subroutine get_data()

    !leggi la directory
    call system("if [ ! -d "//adjustl(trim(SRC))//" ]; then touch "//trim(adjustl(trim(FFILE)))//" ;fi")
    inquire(file=trim(adjustl(trim(FFILE))),EXIST=file)
    if(file)then
       call system("rm -f "//trim(adjustl(trim(FFILE))))
       print*,"----------------------------------------------------------------------"
       write(*,"(A)")"Directory "//trim(adjustl(trim(SRC)))//" does not exist: skipping"
       print*,"----------------------------------------------------------------------"
       print*,""
       return

    elseif(.not.file)then
       write(*,"(A)")"Results are in ",trim(adjustl(trim(DIR)))
       call system("if [ ! -d "//trim(adjustl(trim(DIR)))//" ]; then mkdir "//trim(adjustl(trim(DIR)))//"; fi")

       call massive_allocation()

       call get_Bath

       !Read the functions:
       call read_locgf(trim(adjustl(trim(SRC))))
       call read_sigma(trim(adjustl(trim(SRC))))
       call read_weissfld(trim(adjustl(trim(SRC))))
       call read_impgf(trim(adjustl(trim(SRC))))
       call read_nk(trim(adjustl(trim(SRC))))

       forall(i=0:nstep,j=0:nstep)
          locGret(i,j)   = heaviside(t(i)-t(j))*(locGgtr(i,j) - locGless(i,j))
          Sret(i,j)      = heaviside(t(i)-t(j))*(Sgtr(i,j)    - Sless(i,j))
          impGret(i,j)   = heaviside(t(i)-t(j))*(impGgtr(i,j) - impGless(i,j))
          G0ret(i,j)     = heaviside(t(i)-t(j))*(G0gtr(i,j)   - G0less(i,j))
       end forall

       call get_plot_quasiWigner_functions( trim(adjustl(trim(DIR))) )

       call evaluate_print_observables( trim(adjustl(trim(DIR))) )

       call massive_deallocation
    endif
    return
  end subroutine get_data
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine massive_allocation
    allocate(G0gtr(0:nstep,0:nstep),G0less(0:nstep,0:nstep))
    allocate(S0gtr(-nstep:nstep),S0less(-nstep:nstep))
    allocate(Sgtr(0:nstep,0:nstep),Sless(0:nstep,0:nstep))
    allocate(locGless(0:nstep,0:nstep),locGgtr(0:nstep,0:nstep))
    allocate(impGless(0:nstep,0:nstep),impGgtr(0:nstep,0:nstep))
    allocate(locGret(0:nstep,0:nstep))
    allocate(Sret(0:nstep,0:nstep))
    allocate(impGret(0:nstep,0:nstep))
    allocate(G0ret(0:nstep,0:nstep))
    allocate(nk(0:nstep,Lk))
    allocate(g0fret(2*nstep),g0fless(2*nstep),g0fgtr(2*nstep))
    allocate(gfret(2*nstep),gfless(2*nstep),gfgtr(2*nstep))
    allocate(sfret(2*nstep))
    allocate(g0tret(-nstep:nstep),g0tless(-nstep:nstep),g0tgtr(-nstep:nstep))
    allocate(gtret(-nstep:nstep),gtless(-nstep:nstep),gtgtr(-nstep:nstep))
    allocate(stret(-nstep:nstep),stless(-nstep:nstep),stgtr(-nstep:nstep))
  end subroutine massive_allocation
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine massive_deallocation
    deallocate(G0gtr,G0less)
    deallocate(S0gtr,S0less)
    deallocate(Sgtr,Sless)
    deallocate(locGless,locGgtr)
    deallocate(impGless,impGgtr)
    deallocate(locGret)
    deallocate(Sret)
    deallocate(impGret)
    deallocate(G0ret)
    deallocate(nk)
    deallocate(g0fret,g0fless,g0fgtr)
    deallocate(gfret,gfless,gfgtr)
    deallocate(sfret)
    deallocate(g0tret,g0tless,g0tgtr)
    deallocate(gtret,gtless,gtgtr)
    deallocate(stret,stless,stgtr)
  end subroutine massive_deallocation
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine read_locgf(dir)
    character(len=*) :: dir
    logical          :: file
    inquire(file=dir//"/locGless_data",EXIST=file)
    if(file)call read_data(locGless(0:nstep,0:nstep),dir//"/locGless_data",nstep)
    inquire(file=dir//"/locGgtr_data",EXIST=file)
    if(file)call read_data(locGgtr(0:nstep,0:nstep),dir//"/locGgtr_data",nstep)
  end subroutine read_locgf
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine read_sigma(dir)
    character(len=*) :: dir
    logical          :: file
    inquire(file=dir//"/Sless_data",EXIST=file)
    if(file)call read_data(Sless(0:nstep,0:nstep),dir//"/Sless_data",nstep)
    inquire(file=dir//"/Sgtr_data",EXIST=file)
    if(file)call read_data(Sgtr(0:nstep,0:nstep),dir//"/Sgtr_data",nstep)
  end subroutine read_sigma
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine read_weissfld(dir)
    character(len=*) :: dir
    logical          :: guess
    logical          :: file
    inquire(file=dir//"/guessG0less_data",EXIST=guess)
    if(guess)then
       call read_data(G0less(0:nstep,0:nstep),dir//"/guessG0less_data",nstep)
    else
       inquire(file=dir//"/G0less_data",EXIST=file)
       if(file)call read_data(G0less(0:nstep,0:nstep),dir//"/G0less_data",nstep)
    endif
    !
    inquire(file=dir//"/guessG0gtr_data",EXIST=guess)
    if(guess)then
       call read_data(G0gtr(0:nstep,0:nstep),dir//"/guessG0gtr_data",nstep)
    else
       inquire(file=dir//"/G0gtr_data",EXIST=file)
       if(file)call read_data(G0gtr(0:nstep,0:nstep),dir//"/G0gtr_data",nstep)
    endif
  end subroutine read_weissfld
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine read_impgf(dir)
    character(len=*) :: dir
    logical          :: file
    inquire(file=dir//"/impGless_data",EXIST=file)
    if(file)call read_data(impGless(0:nstep,0:nstep),dir//"/impGless_data",nstep)
    inquire(file=dir//"/impGgtr_data",EXIST=file)
    if(file)call read_data(impGgtr(0:nstep,0:nstep),dir//"/impGgtr_data",nstep)
  end subroutine read_impgf
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine read_nk(dir)
    character(len=*) :: dir
    logical          :: file
    inquire(file=dir//"/nk_data",EXIST=file)
    if(file)then
       open(20,file=dir//"/nk_data")
       do i=0,nstep
          do ik=1,Lk
             read(20,*)nk(i,ik)
          enddo
       enddo
       close(20)
    endif
  end subroutine read_nk
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine read_data(Gin,inFILE,M)
    character(len=*)              :: inFILE
    integer                       :: M,i,j
    complex(8),dimension(0:M,0:M) :: Gin
    open(20,file=trim(inFILE))
    do i=0,M
       do j=0,M
          read(20,*)Gin(i,j)
       enddo
    enddo
    close(20)
  end subroutine read_data
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  subroutine evaluate_print_observables(dir)
    character(len=*)                          :: dir
    integer                                   :: i,ik,ix,iy,it,step
    complex(8)                                :: I1,Ib
    real(8)                                   :: Wtot
    type(vect2D)                              :: Ak,kt,Jk
    type(vect2D),dimension(0:nstep)           :: Jloc,Jheat !local Current 
    type(vect2D),dimension(0:nstep,0:Nx,0:Ny) :: Jkvec,Tloc                  !current vector field
    real(8),dimension(0:nstep,Lk)             :: npi                    !covariant occupation n(\pi=\ka+\Ekt) 
    real(8),dimension(0:nstep)                :: nt,Jint,Stot   !occupation(time)
    real(8),dimension(0:Nx,0:Ny,0:nstep)      :: nDens                  !occupation distribution on the k-grid
    real(8),dimension(0:nstep)                :: Ekin,Epot,Eb,Etot,doble!energies and double occ.
    real(8),dimension(0:nstep,Lk)             :: sorted_nk,sorted_npi   !sorted arrays
    real(8),dimension(0:nstep,Lk)             :: epsikt,sorted_epsikt   !sorted arrays
    real(8),dimension(0:nstep,Lk)             :: sorted_epsipit         !sorted arrays
    integer,dimension(:),allocatable          :: reduced_ik             !reduced arrays
    real(8),dimension(:),allocatable          :: reduced_epsik
    real(8),dimension(:,:),allocatable        :: reduced_nk,reduced_npi,reduced_epsipit


    call dump("Print Out Results (may take a while):")
    !SORTING:
    write(*,"(A)")"Sorting:"
    forall(i=0:nstep,ik=1:Lk)sorted_nk(i,ik) = nk(i,sorted_ik(ik))


    !COVARIANT transformation: \ka --> \pi = \ka + \Ek*t
    write(*,"(A)")"Covariant transf.:"
    call shift_kpoint(nk(0:nstep,1:Lk), npi(0:nstep,1:Lk))
    call shift_kpoint(sorted_nk(0:nstep,1:Lk), sorted_npi(0:nstep,1:Lk))


    !REDUCTION of the k-grid:
    write(*,"(A)")"Reducing BZ:"
    step=Lk/Lkreduced; if(step==0)step=1
    call get_reduxGrid_dimension(Lk,step,Lkreduced)
    allocate(reduced_ik(Lkreduced),reduced_epsik(Lkreduced),reduced_epsipit(0:nstep,Lkreduced),&
         reduced_nk(0:nstep,Lkreduced),reduced_npi(0:nstep,Lkreduced))
    call get_reduxGrid_index(Lk,step,reduced_ik)
    call get_reduxGrid_epsik(sorted_epsik,reduced_ik,reduced_epsik)
    forall(ik=1:Lkreduced)reduced_nk(0:nstep,ik)  = sorted_nk(0:nstep,reduced_ik(ik))
    forall(ik=1:Lkreduced)reduced_npi(0:nstep,ik) = sorted_npi(0:nstep,reduced_ik(ik))

    !Get the CURRENT Jloc(t)=\sum_\ka J_\ka(t) = -e n_\ka(t)*v_\ka(t)
    write(*,"(A)")"Current Field:"
    Jloc=Vzero ; Jheat=Vzero   ;Stot=0.d0
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=0,nstep
          Ak          = Afield(t(i),Ek)
          kt          = kgrid(ix,iy)-Ak
          epsikt(i,ik)= epsk(kt)
          Jk          = nk(i,ik)*velocity(kt)
          Jkvec(i,ix,iy)  = Jk
          Jloc(i)         = Jloc(i) +  wt(ik)*Jk
          Jheat(i)        = Jheat(i)+  wt(ik)*epsikt(i,ik)*Jk
          Stot(i)         = Stot(i) -  wt(ik)*(nk(i,ik)*log(nk(i,ik))+(1.d0-nk(i,ik))*log(1.d0-nk(i,ik)))
       enddo
    enddo
    Stot(0)=0.d0
    forall(i=0:nstep,ik=1:Lk)sorted_epsikt(i,ik) = epsikt(i,sorted_ik(ik))
    call shift_kpoint(sorted_epsikt(0:nstep,1:Lk), sorted_epsipit(0:nstep,1:Lk))
    forall(ik=1:Lkreduced)reduced_epsipit(0:nstep,ik) = sorted_epsikt(0:nstep,reduced_ik(ik))

    !OCCUPATION :
    write(*,"(A)")"Get occupation"
    forall(i=0:nstep)nt(i)=-xi*locGless(i,i)

    !OCCUPATION density:
    forall(i=0:nstep,ik=1:Lk)nDens(ik2ix(ik),ik2iy(ik),i)=npi(i,ik)

    !ENERGY kinetic
    Ekin=zero
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=0,nstep
          Ak=Afield(t(i),Ek)
          Ekin(i) = Ekin(i) +  wt(ik)*epsk(kgrid(ix,iy) - Ak)*nk(i,ik)
       enddo
    enddo


    !ENERGY potential && total    if(Efield/=0.d0 .AND. plotVF)call plot_VF("vf_JfieldVSkVSt",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,Jkvec(0:nstep,:,:)%x,Jkvec(0:nstep,:,:)%y)
    !Get Epot = <V>(it)= xi/2 \lim_{t'-->t} \sum_k [xi\partial_t - h_{0,k}(t)] G^<_k(t,t') 
    ! by eq. of motion = xi/2 \lim_{t'-->t} \sum_k {\delta^<(t,t') + \int_0^t S^R(t,z)*G^<_k(z,t') + \int_0^t' S^<(t,z)*G^A_k(z,t')}
    !                  = xi/2 \lim_{t'-->t} {\delta^<(t,t') + \int_0^t S^R(t,z)*G^<_loc(z,t') + \int_0^t' S^<(t,z)*G^A_loc(z,t')}
    !                  = xi/2 {1 + \int_0^t [S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*G^A_loc(z,t)]}
    ! cnst disregarded = xi/2 \int_0^t {S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*[G^R_loc(t,z)]^+}
    do it=0,nstep
       I1=zero; Ib=zero
       do i=1,it
          I1 = I1 + SretF(it,i)*locGless(i,it) + Sless(it,i)*conjg(GretF(it,i))
          Ib = Ib + S0retF(it-i)*locGless(i,it) + S0less(it-i)*conjg(GretF(it,i))
       enddo
       Epot(it)= -xi*I1*dt/2.d0
       Eb(it)= -xi*Ib*dt/2.d0
    enddo
    !Get Etot = Ekin + Epot
    Etot=Ekin+Epot

    Jint=Jloc%x + Jloc%y
    Wtot=sum(Jint(0:))*Efield*dt

    !Double OCCUPATION:
    doble= 0.5d0*(2.d0*nt) - 0.25d0 ; if(U/=0)doble = Epot/U + 0.5d0*(2.d0*nt)- 0.25d0


    !PRINT:
    write(*,"(A)")"Print n(t):"
    call splot(dir//"/nVStime.ipt",t(0:nstep),2.d0*nt(0:nstep))
    !call plot("nVStime",t(0:nstep),2.d0*nt(0:nstep),Xlabel="t",Ylabel="<n>",wlp="wlp")


    write(*,"(A)")"Print J(t):"
    if(Efield/=0.d0)then
       call splot(dir//"/JlocVStime.ipt",t(0:nstep),Jloc(0:nstep)%x+Jloc(0:nstep)%y)
       call splot(dir//"/JheatVStime.ipt",t(0:nstep),Jheat(0:nstep)%x+Jheat(0:nstep)%y)
       !call plot("JlocVStime",t(0:nstep),Jloc(0:nstep)%x+Jloc(0:nstep)%y,"",Xlabel="t",Ylabel="<J>",wlp="wlp")
    endif

    write(*,"(A)")"Print Ex(t)"
    call splot(dir//"/EkinVStime.ipt",t(0:nstep),Ekin(0:nstep))
    call splot(dir//"/EpotVStime.ipt",t(0:nstep),Epot(0:nstep))
    call splot(dir//"/EhybVStime.ipt",t(0:nstep),Eb(0:nstep))
    call splot(dir//"/EtotVStime.ipt",t(0:nstep),Etot(0:nstep),Etot(0:nstep)+Eb(0:nstep))
    call splot(dir//"/WtotVSefield.ipt",Efield,Wtot)
    call splot(dir//"/StotVStime.ipt",t(0:nstep),Stot(0:nstep))
    !call plot("EkinVStime",t(0:nstep),Ekin(0:nstep),Xlabel="t",Ylabel="<K>",wlp="wlp")
    !call plot("EpotVStime",t(0:nstep),Epot(0:nstep),Xlabel="t",Ylabel="<U>",wlp="wlp")
    !call plot("EtotVStime",t(0:nstep),Etot(0:nstep),Xlabel="t",Ylabel="<E>",wlp="wlp")

    write(*,"(A)")"Print d(t)"
    call splot(dir//"/doccVStime.ipt",t(0:nstep),doble(0:nstep))
    !call plot("doccVStime",t(0:nstep),doble(0:nstep),Xlabel="t",Ylabel="<d>",wlp="wlp")


    !DISTRIBUTION:
    write(*,"(A)")"Print n(k,t):"
    call splot("nVStimeVSepsk3D.ipt",t(0:nstep),reduced_epsik,reduced_nk(0:nstep,:))
    !call splot("nVStimeVSepspit3D.ipt",t(0:nstep),reduced_epsik,reduced_npi(0:nstep,:),TT)
    ! do i=0,nstep
    !    call splot("nVSepsik.ipt",reduced_epsik,reduced_npi(i,:),append=TT)
    ! enddo

    !Fermi Surface plot:
    if(Efield/=0.d0 .or. Vpd/=0.0)then
       call plot_3D_movie("3dFSVSpiVSt","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:nstep))
    else
       call plot_3D("FSVSpi3D","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,nstep))
    endif

    !Current Vector Field:
    if(Efield/=0.d0 .AND. plotVF)call plot_VF("vf_JfieldVSkVSt",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,Jkvec(0:nstep,:,:)%x,Jkvec(0:nstep,:,:)%y)

    !Local functions:
    !===========================================================================
    if(plot3D)then
       if(Efield/=0.d0 .or. Vpd/=0.0)call plot_3D_surface_movie("3dFSVSpiVSt","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:nstep))
       call plot_movie("gif_nVSepsipiVSt",reduced_epsik,reduced_npi(0:nstep,:),Xlabel="$\Huge\epsilon(k)$",Ylabel="$\Huge n_{\pi}(t)$",wlp="wlp")
       if(loop==1)then
          call plot_3D("guessG0less3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
          call plot_3D("guessG0gtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0gtr(0:nstep,0:nstep))
       else
          call plot_3D("G0less3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
          call plot_3D("G0gtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0gtr(0:nstep,0:nstep))
       endif
       call plot_3D("locGless3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,locGless)
       call plot_3D("locGgtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,locGgtr)
       call plot_3D("impGless3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,impGless)
       call plot_3D("impGgtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,impGgtr)
       call plot_3D("Sless3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,Sless)
       call plot_3D("Sgtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,Sgtr)
       !call system("mv -vf *3D "//dir//"/ 2>/dec/null")
    endif

    call system("rm -rf "//dir//"/3d* 2>/dev/null")
    call system("rm -rf "//dir//"/*3D 2>/dev/null")
    call system("rm -rf "//dir//"/gif_* 2>/dev/null")
    call system("mv -vf PNG AGR gif_* vf_* 3d* *3D* "//dir//"/ 2>/dev/null")

    return
  end subroutine evaluate_print_observables
  !+-------------------------------------------------------------------+



  !+-------------------------------------------------------------------+
  subroutine get_plot_quasiWigner_functions(dir)
    character(len=*)  :: dir
    integer           :: i,j
    forall(i=0:nstep,j=0:nstep)
       gtless(i-j) = locGless(i,j)
       gtgtr(i-j)  = locGgtr(i,j)
       gtret(i-j)  = locGret(i,j)
    end forall
    if(heaviside(0.d0)==1.d0)gtret(0)=gtret(0)/2.d0
    call cfft_rt2rw(gtret,gfret,nstep) ;    gfret=gfret*dt ; call swap_fftrt2rw(gfret)
    call splot(dir//"/locGless_t.ipt",t(-nstep:nstep),gtless,TT)
    call splot(dir//"/locGgtr_t.ipt",t(-nstep:nstep),gtgtr,TT)
    call splot(dir//"/locGret_t.ipt",t(-nstep:nstep),gtret,TT)
    call splot(dir//"/locGret_realw.ipt",wr,gfret,TT)
    call splot(dir//"/locDOS.ipt",wr,-aimag(gfret)/pi,append=TT)

    forall(i=0:nstep,j=0:nstep)
       g0tless(i-j)= G0less(i,j)
       g0tgtr(i-j) = G0gtr(i,j)
       g0tret(i-j) = G0ret(i,j)
       gtless(i-j) = impGless(i,j)
       gtgtr(i-j)  = impGgtr(i,j)
       gtret(i-j)  = impGret(i,j)
       stless(i-j) = Sless(i,j)
       stgtr(i-j)  = Sgtr(i,j)
       stret(i-j)  = Sret(i,j)
    end forall
    if(heaviside(0.d0)==1.d0)g0tret(0)=g0tret(0)/2.d0 ; g0tret(0)=-xi
    if(heaviside(0.d0)==1.d0)stret(0)=stret(0)/2.d0
    if(heaviside(0.d0)==1.d0)gtret(0)=gtret(0)/2.d0   !; gtret(0)=-xi
    if(loop==1)then
       call splot(dir//"/guessG0less_t.ipt",t(-nstep:nstep),g0tless,append=TT)
       call splot(dir//"/guessG0gtr_t.ipt",t(-nstep:nstep),g0tgtr,append=TT)
       call splot(dir//"/guessG0ret_t.ipt",t(-nstep:nstep),g0tret,append=TT)
    else
       call splot(dir//"/G0less_t.ipt",t(-nstep:nstep),g0tless,append=TT)
       call splot(dir//"/G0gtr_t.ipt",t(-nstep:nstep),g0tgtr,append=TT)
       call splot(dir//"/G0ret_t.ipt",t(-nstep:nstep),g0tret,append=TT)
    endif
    call splot(dir//"/impGless_t.ipt",t(-nstep:nstep),gtless,append=TT)
    call splot(dir//"/impGgtr_t.ipt",t(-nstep:nstep),gtgtr,append=TT)
    call splot(dir//"/impGret_t.ipt",t(-nstep:nstep),gtret,append=TT)
    call splot(dir//"/Sless_t.ipt",t(-nstep:nstep),stless,append=TT)
    call splot(dir//"/Sgtr_t.ipt",t(-nstep:nstep),stgtr,append=TT)
    call splot(dir//"/Sret_t.ipt",t(-nstep:nstep),stret,append=TT)


    !Obtain && plot Real frequency Functions:
    !===========================================================================
    call cfft_rt2rw(g0tret,g0fret,nstep) ;    g0fret=g0fret*dt ; call swap_fftrt2rw(g0fret)
    call cfft_rt2rw(gtret,gfret,nstep)   ;    gfret=gfret*dt   ; call swap_fftrt2rw(gfret)
    call cfft_rt2rw(stret,sfret,nstep)   ;    sfret=dt*sfret   ; call swap_fftrt2rw(sfret)
    if(loop==1)then
       call splot(dir//"/guessG0ret_realw.ipt",wr,g0fret,append=TT)
    else
       call splot(dir//"/G0ret_realw.ipt",wr,g0fret,append=TT)
    endif
    call splot(dir//"/impGret_realw.ipt",wr,gfret,append=TT)
    call splot(dir//"/Sret_realw.ipt",wr,sfret,append=TT)
    call splot(dir//"/DOS.ipt",wr,-aimag(gfret)/pi,append=TT)
    return
  end subroutine get_plot_quasiWigner_functions
  !+-------------------------------------------------------------------+







  !+-------------------------------------------------------------------+
  subroutine shift_kpoint(arrayIN,arrayOUT)
    integer                 :: i,j,ik,ix,iy,jk,jx,jy
    real(8),dimension(0:,:) :: arrayIN,arrayOUT
    real(8),dimension(2)    :: pi_in
    integer,dimension(2)    :: pi_kcoord
    type(vect2D)            :: Ak
    arrayOUT=0.d0
    do i=0,nstep
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



  !+-------------------------------------------------------------------+
  function GretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: GretF
    GretF = heaviside(t(i)-t(j))*(locGgtr(i,j)-locGless(i,j))
  end function GretF
  !-------------------------------------------------------!

  !-------------------------------------------------------!
  function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i)-t(j))*(Sgtr(i,j)-Sless(i,j))
  end function SretF
  !-------------------------------------------------------!

  !-------------------------------------------------------!
  function S0retF(i)
    integer,intent(in) :: i
    complex(8)         :: S0retF
    S0retF = heaviside(t(i))*(S0gtr(i)-S0less(i))
  end function S0retF
  !-------------------------------------------------------!

end program getDATA
