!###############################################################
!     PROGRAM  : OBSERVABLES
!     TYPE     : Module
!     PURPOSE  : 
!     AUTHORS  : Adriano Amaricci
!     LAST UPDATE: 08/2010
!###############################################################
module OBSERVABLES
  use VARS_GLOBAL
  implicit none
  private
  !public :: print_observables_simple
  public :: evaluate_print_observables

contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : PRINT_OBSERVABLES_SIMPLE
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine print_observables_simple
    integer                          :: i,ik,ix,iy
    type(vect2D)                     :: Jk,Ak
    type(vect2D),dimension(0:nstep)  :: Jloc                   !local Current 
    real(8),dimension(0:nstep)       :: nt                     !occupation(time)
    forall(i=0:nstep)nt(i)=-xi*locGless(i,i)
    Jloc=Vzero    
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=0,nstep
          Ak= Afield(t(i),Ek)
          Jk= nk(i,ik)*velocity(kgrid(ix,iy) - Ak)
          Jloc(i) = Jloc(i) +  wt(ik)*Jk
       enddo
    enddo
    call splot("nVStime.ipt",t(0:nstep),2.d0*nt(0:nstep),append=TT)
    if(Efield/=0.d0)call splot("JlocVStime.ipt",t(0:nstep),Jloc(0:nstep)%x+Jloc(0:nstep)%y,append=TT)
    return
  end subroutine print_observables_simple




  !+-------------------------------------------------------------------+
  !PROGRAM  : EVALUATE_PRINT_OBSERVABLES
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine evaluate_print_observables(loop)
    integer                                   :: i,ik,ix,iy,it,loop,step
    complex(8)                                :: I1
    type(vect2D)                              :: Ak
    type(vect2D),dimension(0:nstep)           :: Jloc                   !local Current 
    type(vect2D),dimension(0:nstep,0:Nx,0:Ny) :: Jkvec                  !current vector field
    real(8),dimension(0:nstep,Lk)             :: npi                    !covariant occupation n(\pi=\ka+\Ekt) 
    real(8),dimension(0:nstep)                :: nt                     !occupation(time)
    real(8),dimension(0:Nx,0:Ny,0:nstep)      :: nDens                  !occupation distribution on the k-grid
    real(8),dimension(0:nstep)                :: Ekin,Epot,Etot,doble   !energies and double occ.
    !Sorted    
    real(8),dimension(0:nstep,Lk)             :: sorted_nk,sorted_npi
    !Reduced
    integer,dimension(:),allocatable          :: reduced_ik
    real(8),dimension(:),allocatable          :: reduced_epsik
    real(8),dimension(:,:),allocatable        :: reduced_nk,reduced_npi
    character(len=4)                          :: loop_char


    call dump("Print Out Results (may take a while):")

    !SORTING:
    print*,"sorting"
    forall(i=0:nstep,ik=1:Lk)sorted_nk(i,ik) = nk(i,sorted_ik(ik))


    !COVARIANT transformation: \ka --> \pi = \ka + \Ek*t
    print*,"covariant transf."
    call shift_kpoint(nk(0:nstep,1:Lk), npi(0:nstep,1:Lk))
    call shift_kpoint(sorted_nk(0:nstep,1:Lk), sorted_npi(0:nstep,1:Lk))


    !REDUCTION of the k-grid:
    print*,"reduction"
    step=Lk/Lkreduced; if(step==0)step=1
    call get_reduxGrid_dimension(Lk,step,Lkreduced)
    allocate(reduced_ik(Lkreduced),reduced_epsik(Lkreduced),reduced_nk(0:nstep,Lkreduced),reduced_npi(0:nstep,Lkreduced))
    call get_reduxGrid_index(Lk,step,reduced_ik)
    call get_reduxGrid_epsik(sorted_epsik,reduced_ik,reduced_epsik)
    forall(ik=1:Lkreduced)reduced_nk(0:nstep,ik)  = sorted_nk(0:nstep,reduced_ik(ik))
    forall(ik=1:Lkreduced)reduced_npi(0:nstep,ik)  = sorted_npi(0:nstep,reduced_ik(ik))

    !Get the CURRENT Jloc(t)=\sum_\ka J_\ka(t) = -e n_\ka(t)*v_\ka(t)
    print*,"get local current"
    Jloc=Vzero    
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=0,nstep
          Ak=Afield(t(i),Ek)
          Jkvec(i,ix,iy)  = nk(i,ik)*velocity(kgrid(ix,iy) + Ak)
          Jloc(i)         = Jloc(i) +  wt(ik)*Jkvec(i,ix,iy)
       enddo
    enddo


    !OCCUPATION :
    print*,"get occupation(t)"
    forall(i=0:nstep)nt(i)=-xi*locGless(i,i)


    !Obtain && plot REAL TIME Functions
    print*,"enter print GF"
    call get_plot_quasiWigner_functions()


    !PRINT:
    print*,"print n"
    call splot("nVStime.ipt",t(0:nstep),2.d0*nt(0:nstep),append=TT)
    call splot("nVStimeVSepsk3D.ipt",t(0:nstep),reduced_epsik,reduced_nk(0:nstep,:),TT)
    if(Efield/=0.d0)call splot("JlocVStime.ipt",t(0:nstep),Jloc(0:nstep)%x+Jloc(0:nstep)%y,append=TT)


    !More to print at the last loop:
    if(loop==nloop)then
       !OCCUPATION density:
       forall(i=0:nstep,ik=1:Lk)nDens(ik2ix(ik),ik2iy(ik),i)=npi(i,ik)


       !ENERGY kinetic
       Ekin=zero
       do ik=1,Lk
          ix=ik2ix(ik);iy=ik2iy(ik)
          do i=0,nstep
             Ekin(i) = Ekin(i) +  wt(ik)*epsk(kgrid(ix,iy) + t(i)*Ek)*nk(i,ik)
          enddo
       enddo
       !ENERGY potential && total
       !Get Epot = <V>(it)= xi/2 \lim_{t'-->t} \sum_k [xi\partial_t - h_{0,k}(t)] G^<_k(t,t') 
       ! by eq. of motion = xi/2 \lim_{t'-->t} \sum_k {\delta^<(t,t') + \int_0^t S^R(t,z)*G^<_k(z,t') + \int_0^t' S^<(t,z)*G^A_k(z,t')}
       !                  = xi/2 \lim_{t'-->t} {\delta^<(t,t') + \int_0^t S^R(t,z)*G^<_loc(z,t') + \int_0^t' S^<(t,z)*G^A_loc(z,t')}
       !                  = xi/2 {1 + \int_0^t [S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*G^A_loc(z,t)]}
       ! cnst disregarded = xi/2 \int_0^t {S^R(t,z)*G^<_loc(z,t) + S^<(t,z)*[G^R_loc(t,z)]^+}
       do it=0,nstep
          I1=zero
          do i=1,it
             I1 = I1 + SretF(it,i)*locGless(i,it) + Sless(it,i)*conjg(GretF(it,i))
          enddo
          Epot(it)= -xi*I1*dt/2.d0
       enddo
       !Get Etot = Ekin + Epot
       Etot=Ekin+Epot


       !DOUBLE OCCUPATION:
       doble= 0.5d0*(2.d0*nt) - 0.25d0 ; if(U/=0)doble = Epot/U + 0.5d0*(2.d0*nt)- 0.25d0


       call splot("EkinVStime.ipt",t(0:nstep),Ekin(0:nstep),append=TT)
       call splot("EpotVStime.ipt",t(0:nstep),Epot(0:nstep),append=TT)
       call splot("EtotVStime.ipt",t(0:nstep),Etot(0:nstep),append=TT)
       call splot("doccVStime.ipt",t(0:nstep),doble(0:nstep),append=TT)


       !DISTRIBUTION:
       write(loop_char,"(I4)")loop
       do i=0,nstep
          call splot("nVSepsipi_"//trim(adjustl(trim(loop_char)))//".ipt",reduced_epsik,reduced_npi(i,:),append=TT)
       enddo


       !Fermi Surface plot:
       if(Efield/=0.d0 .or. Vpd/=0.0)then
          call plot_3D_movie("FSVSpiVSt","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:nstep))
          call plot_3D_surface_movie("FSVSpiVSt","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,0:nstep))
       else
          call plot_3D("FSVSpi3D","$k_x$","$k_y$","$FS(k_x,k_y)$",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,nDens(0:Nx,0:Ny,nstep))
       endif


       !Current Vector Field:
       if(Efield/=0.d0)then
          call plot("JlocVStime",t(0:nstep),Jloc(0:nstep)%x+Jloc(0:nstep)%y,"",Xlabel="t",Ylabel="<J>",wlp="wlp")
          call plot_VF("JfieldVSkVSt",kgrid(0:Nx,0)%x,kgrid(0,0:Ny)%y,Jkvec(0:nstep,:,:)%x,Jkvec(0:nstep,:,:)%y)
       endif


       !Occupation && Energy && Double Occ.:
       call plot("nVStime",t(0:nstep),2.d0*nt(0:nstep),Xlabel="t",Ylabel="<n>",wlp="wlp")
       call plot("EkinVStime",t(0:nstep),Ekin(0:nstep),Xlabel="t",Ylabel="<K>",wlp="wlp")
       call plot("EpotVStime",t(0:nstep),Y1=Epot(0:nstep),Xlabel="t",Ylabel="<E>",wlp="wlp")
       call plot("EtotVStime",t(0:nstep),Y1=Etot(0:nstep),Xlabel="t",Ylabel="<U>",wlp="wlp")
       call plot("doccVStime",t(0:nstep),doble(0:nstep),Xlabel="t",Ylabel="<d>",wlp="wlp")
       call plot_movie("nVSepsipiVSt",reduced_epsik,reduced_npi(0:nstep,:),Xlabel="$\Huge\epsilon(k)$",Ylabel="$\Huge n_{\pi}(t)$",wlp="wlp")


       !Local functions:
       !===========================================================================
       if(plot3D)then
          call system("if [ ! -d 3Dplot ]; then mkdir 3Dplot; fi")
          call plot_3D("G0less3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
          call plot_3D("G0gtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,G0gtr(0:nstep,0:nstep))
          call plot_3D("locGless3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,locGless)
          call plot_3D("locGgtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,locGgtr)
          call plot_3D("impGless3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,impGless)
          call plot_3D("impGgtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,impGgtr)
          call plot_3D("Sless3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,Sless)
          call plot_3D("Sgtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,Sgtr)
          call system("mv -v *3D 3Dplot/")
       endif

    end if
    return
  end subroutine evaluate_print_observables
  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  :  
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine get_plot_quasiWigner_functions
    integer                                   :: i,j
    forall(i=0:nstep,j=0:nstep)
       stless(i-j) = Sless(i,j)
       stgtr(i-j)  = Sgtr(i,j)
       stret(i-j)  = heaviside(t(i-j))*(Sgtr(i,j)-Sless(i,j))
       g0tless(i-j)= G0less(i,j)
       g0tgtr(i-j) = G0gtr(i,j)
       g0tret(i-j) = heaviside(t(i-j))*(G0gtr(i,j)-G0less(i,j))
       gtless(i-j) = impGless(i,j)
       gtgtr(i-j)  = impGgtr(i,j)
       gtret(i-j)  = heaviside(t(i-j))*(impGgtr(i,j)-impGless(i,j))
    end forall
    if(heaviside(0.d0)==1.d0)stret(0)=stret(0)/2.d0
    if(heaviside(0.d0)==1.d0)g0tret(0)=g0tret(0)/2.d0 ; g0tret(0)=-xi
    if(heaviside(0.d0)==1.d0)gtret(0)=gtret(0)/2.d0   !; gtret(0)=-xi
    call splot("Sless_t.ipt",t(-nstep:nstep),stless,append=TT)
    call splot("Sgtr_t.ipt",t(-nstep:nstep),stgtr,append=TT)
    call splot("Sret_t.ipt",t(-nstep:nstep),stret,append=TT)
    call splot("G0less_t.ipt",t(-nstep:nstep),g0tless,append=TT)
    call splot("G0gtr_t.ipt",t(-nstep:nstep),g0tgtr,append=TT)
    call splot("G0ret_t.ipt",t(-nstep:nstep),g0tret,append=TT)
    call splot("impGless_t.ipt",t(-nstep:nstep),gtless,append=TT)
    call splot("impGgtr_t.ipt",t(-nstep:nstep),gtgtr,append=TT)
    call splot("impGret_t.ipt",t(-nstep:nstep),gtret,append=TT)


    !Obtain && plot Real frequency Functions:
    !===========================================================================
    call cfft_rt2rw(stret,sfret,nstep)   ;    sfret=dt*sfret   ; call swap_fftrt2rw(sfret)
    call cfft_rt2rw(g0tret,g0fret,nstep) ;    g0fret=g0fret*dt ; call swap_fftrt2rw(g0fret)
    call cfft_rt2rw(gtret,gfret,nstep)   ;    gfret=gfret*dt   ; call swap_fftrt2rw(gfret)
    call splot("Sret_realw.ipt",wr,sfret,append=TT)
    call splot("G0ret_realw.ipt",wr,g0fret,append=TT)
    call splot("impGret_realw.ipt",wr,gfret,append=TT)
    call splot("DOS.ipt",wr,-aimag(gfret)/pi,append=TT)

    return
  end subroutine get_plot_quasiWigner_functions
  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : SHIFT_KPOINT 
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
          pi_in(2)=kgrid(ix,iy)%y + (-t(i))*Ek%y
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
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  :  
  !TYPE     : FUNCTIONS
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function GretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: GretF
    GretF = heaviside(t(i)-t(j))*(locGgtr(i,j)-locGless(i,j))
  end function GretF
  !-------------------------------------------------------!
  function S0retF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: S0retF
    S0retF = heaviside(t(i)-t(j))*(S0gtr(i-j)-S0less(i-j))
  end function S0retF
  !-------------------------------------------------------!
  function SretF(i,j)      
    integer,intent(in) :: i,j
    complex(8)         :: SretF
    SretF = heaviside(t(i)-t(j))*(Sgtr(i,j)-Sless(i,j))
  end function SretF
  !******************************************************************
  !******************************************************************
  !******************************************************************

end module OBSERVABLES
