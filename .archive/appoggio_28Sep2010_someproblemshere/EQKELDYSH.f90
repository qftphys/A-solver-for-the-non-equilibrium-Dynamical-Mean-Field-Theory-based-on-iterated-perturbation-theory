MODULE EQKELDYSH
  !####################################################################
  !PROGRAM  : EQ_SOLUTION
  !TYPE     : module
  !PURPOSE  : Solve the model at EQUILIBRIUM using DMFT and Keldysh technique
  !AUTHORS  : Adriano Amaricci
  !INPUT    : 
  ! - neqloop = number of DMFT loop 
  !OUTPUT   : 
  ! - fg0   = real frequency Weiss field \cal{G}_0(\w) [BATH GF]
  ! - Sigma = Sigma_Adv on real axis 
  !####################################################################
  !LIBRARY:
  !===================
  USE TOOLS
  USE DLPLOT
  USE LATTICE
  USE PROGRESS
  USE FFTW
  USE GRIDS
  USE SPLINE
  !LOCAL:
  !===================
  USE VARS_GLOBAL
  implicit none

  !Greens functions
  !=========================================================
  private
  integer                          :: i,j,ik,im
  real(8)                          :: rfmesh
  real(8),dimension(:),allocatable :: eqwr,eqwm,eqt,ex
  !Frequency domain:
  complex(8),dimension(:),allocatable :: eqg0fret,eqg0fless,eqg0fgtr
  complex(8),dimension(:),allocatable :: eqgfret,eqgfless,eqgfgtr
  complex(8),dimension(:),allocatable :: eqsfret
  complex(8),dimension(:),allocatable :: eqsbfret
  !Time domain:
  complex(8),dimension(:),allocatable :: eqg0tret,eqg0tless,eqg0tgtr
  complex(8),dimension(:),allocatable :: eqgtret,eqgtless,eqgtgtr
  complex(8),dimension(:),allocatable :: eqstret,eqstless,eqstgtr
  !Matsubara
  complex(8),dimension(:),allocatable :: GMiw,G0Miw,SMiw,SbMiw
  !Tau
  real(8),allocatable,dimension(:)    :: Gt,G0t,St
  real(8),dimension(0:Ltau)           :: Gtau,G0tau,Stau

  public :: keldysheq

contains 
  subroutine keldysheq(neqloop)
    integer :: ieqloop,neqloop
    real(8) :: exx,wini
    !=========================================================


    write(*,*)"Get solution at EQUILIBRIUM "
    write(*,*)"----------------------------"

    !Init calculation
    !===============================================================
    call system("mkdir EQSOLUTION")   
    allocate(eqg0fret(2*L),eqg0fless(2*L),eqg0fgtr(2*L))
    allocate(eqgfret(2*L),eqgfless(2*L),eqgfgtr(2*L))
    allocate(eqsfret(2*L),eqsbfret(2*L))
    allocate(eqg0tret(-L:L),eqg0tless(-L:L),eqg0tgtr(-L:L))
    allocate(eqgtret(-L:L),eqgtless(-L:L),eqgtgtr(-L:L))
    allocate(eqstret(-L:L),eqstless(-L:L),eqstgtr(-L:L))
    allocate(GMiw(2*L),G0Miw(2*L),SMiw(2*L),SbMiw(2*L))
    allocate(Gt(0:L),G0t(0:L),St(0:L))


    rfmesh  = pi/dt/dble(L)
    wini    = -rfmesh*dble(L)
    allocate(eqwr(2*L),eqwm(2*L),eqt(-L:L))
    call init_wgrid(eqwr,wini,rfmesh)
    call init_wmgrid(eqwm,beta)
    call init_tgrid(eqt,dt,L)
    allocate(ex(-L:L))
    exx=-1.d0       
    do i=-L,L
       exx=-exx
       ex(i)=exx
    enddo
    write(*,"(A,f11.2)")'    U = ',U
    write(*,"(A,f11.2)")'    V = ',Vpd
    write(*,"(A,f11.2)")' beta = ',beta
    write(*,"(A,f11.6)")' Mesh = ',rfmesh

    !Set the contribution from the Bath:
    !===============================================================
    SbMiw=zero ;  eqsbfret  =zero
    !if(Vpd/=0.d0)call eq_build_bath
    call eq_guess_G0

    !Starts DMFT-loop
    eqg0fret   = zero; eqgfret    = zero ; eqsfret= zero
    do ieqloop=1,neqloop
       call print_bar(ieqloop,neqloop)

       !Impurity Solver: Sigma=U^2G0^3 && FFT t-->w
       !===============================================================
       call eq_Simpurity

       !Get G_loc^R=\sum_\kaG(\ka,w) &&  G_loc^{<,>}  && Guess/Update G0^R
       !===============================================================
       call eq_get_Gloc

       !Update WF G0 (SCC)
       !===============================================================
       call eq_update_G0

       call eq_dump_out(ieqloop,neqloop)
    enddo                     !here the dmft loops end up
    call close_bar

    !Store some arrays to provide correlated initial conditions:
    icG0tgtr  = ex*eqg0tgtr  !ex(-nstep:nstep)*eqg0tgtr(-nstep:nstep)
    icG0tless = ex*eqg0tless !ex(-nstep:nstep)*eqg0tless(-nstep:nstep)
    icG0w     = eqg0fret
    icSiw     = SMiw
    icGiw     = GMiw
    call cfft_iw2tau(SMiw,St,beta) ; call extract(St,Stau)
    icStau(0:Ltau) = Stau(0:Ltau)  ; forall(i=1:Ltau)icStau(-i)=-Stau(Ltau-i)
    icGtau(0:Ltau) = Gtau(0:Ltau)  ; forall(i=1:Ltau)icGtau(-i)=-Gtau(Ltau-i)

    deallocate(eqg0fret,eqg0fless,eqg0fgtr)
    deallocate(eqgfret,eqgfless,eqgfgtr)
    deallocate(eqsfret,eqsbfret)
    deallocate(eqg0tret,eqg0tless,eqg0tgtr)
    deallocate(eqgtret,eqgtless,eqgtgtr)
    deallocate(eqstret,eqstless,eqstgtr)
    deallocate(ex)
    deallocate(GMiw,G0Miw,SMiw,SbMiw)
    deallocate(Gt,G0t,St)
    deallocate(eqwr,eqwm,eqt)
  end subroutine keldysheq
  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  : GETBATH
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine eq_build_bath
    integer    :: im
    real(8)    :: w
    complex(8) :: iw

    eqsbfret  =zero
    do i=1,2*L
       w=eqwr(i);iw=cmplx(w,eps)
       do im=1,Lmu
          eqsbfret(i)=eqsbfret(i) + one/(iw-epsimu(im))/dble(Lmu)
       enddo
    enddo
    eqsbfret = Vpd**2*eqsbfret

    SbMiw=zero
    do i=1,2*L
       w=eqwm(i);iw=cmplx(0.d0,w)
       do im=1,Lmu
          SbMiw(i)=SbMiw(i)+Vpd**2/(iw-epsimu(im))/dble(Lmu)
       enddo
    enddo

  end subroutine eq_build_bath
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : GUESS_G0GL
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine eq_guess_G0()
    integer    :: ik
    real(8)    :: en,ngtr,nless
    complex(8) :: peso
    eqg0tless=zero ; eqg0tgtr  =zero
    do ik=1,Lk
       en=epsik(ik)
       ngtr=fermi(en,beta)-1.d0
       nless=fermi(en,beta)
       do i=-L,L
          peso=exp(-xi*en*eqt(i))
          eqg0tless(i)=eqg0tless(i) + xi*nless*peso*wt(ik)
          eqg0tgtr(i) =eqg0tgtr(i)  + xi*ngtr*peso*wt(ik)
       enddo
    enddo
    eqg0tless = ex*eqg0tless
    eqg0tgtr  = ex*eqg0tgtr
    do i=-L,L
       eqg0tret(i)=heaviside(eqt(i))*(eqg0tgtr(i) - eqg0tless(i))
    enddo
    call cfft_rt2rw(eqg0tret,eqg0fret,L) ;eqg0fret=dt*eqg0fret
  end subroutine eq_guess_G0
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PROGRAM  : eq_Simpurity
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine eq_Simpurity
    do i=-L,L
       eqstgtr(i)=(U**2)*(eqg0tgtr(i)**2)*eqg0tless(-i) 
       eqstless(i)=(U**2)*(eqg0tless(i)**2)*eqg0tgtr(-i)
       eqstret(i)=heaviside(eqt(i))*(eqstgtr(i) - eqstless(i))
    enddo
    if(heaviside(0.d0)==1.d0)eqstret(0)=eqstret(0)/2.d0
    call cfft_rt2rw(eqstret,eqsfret,L) ; eqsfret=dt*eqsfret
  end subroutine eq_Simpurity
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PROGRAM  : GET_GLOC
  !TYPE     : function
  !PURPOSE  : Evaluate the k-sum to get the local Green's functions (R,>,<)
  !+-------------------------------------------------------------------+
  subroutine eq_get_Gloc()
    complex(8) :: A,zetan
    real(8)    :: w
    eqgfret=zero
    do i=1,2*L
       w=eqwr(i)
       zetan=cmplx(w,eps)-eqsfret(i)-eqsbfret(i)
       do ik=1,Lk
          eqgfret(i)=eqgfret(i)+wt(ik)/(zetan-epsik(ik))
       enddo
       A=-aimag(eqgfret(i))/pi
       eqgfless(i)= pi2*xi*fermi(w,beta)*A
       eqgfgtr(i) = pi2*xi*(fermi(w,beta)-1.d0)*A
    enddo
    call cfft_rw2rt(eqgfless,eqgtless,L)  ; eqgtless=rfmesh/pi2*eqgtless
    call cfft_rw2rt(eqgfgtr,eqgtgtr,L)    ; eqgtgtr=rfmesh/pi2*eqgtgtr
    do i=-L,L
       eqgtret(i)=heaviside(eqt(i))*(eqgtgtr(i) - eqgtless(i))        
    enddo
    if(heaviside(0.d0)==1.d0)eqgtret(0)=eqgtret(0)/2.d0 ;eqgtret(0)=-xi
    return
  end subroutine eq_get_Gloc
  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  : GETBATH
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine eq_update_G0()
    complex(8) :: A,An
    real(8)    :: w
    eqg0fret=one/(one/eqgfret + eqsfret)
    do i=1,2*L
       w = eqwr(i)
       A = -aimag(eqg0fret(i))/pi
       An= A*fermi(w,beta)
       eqg0fless(i)= pi2*xi*An
       eqg0fgtr(i) = pi2*xi*(An-A)
    enddo
    call cfft_rw2rt(eqg0fless,eqg0tless,L)  ; eqg0tless=rfmesh/pi2*eqg0tless
    call cfft_rw2rt(eqg0fgtr, eqg0tgtr,L)   ; eqg0tgtr =rfmesh/pi2*eqg0tgtr
    do i=-L,L
       eqg0tret(i)=heaviside(eqt(i))*(eqg0tgtr(i) - eqg0tless(i))        
    enddo
    if(heaviside(0.d0)==1.d0)eqg0tret(0)=eqg0tret(0)/2.d0 ; eqg0tret(0)=-xi
  end subroutine eq_update_G0
  !******************************************************************
  !******************************************************************
  !******************************************************************








  !+----------------------------------------------------------------+
  !PROGRAM  : DUMP_OUT
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine eq_dump_out(loop,nloop)
    integer    :: i,ik,loop,nloop
    real(8)    :: w,Ekin,Epot,docc,nimp
    complex(8) :: iw

    call splot('EQSOLUTION/DOS.ipt',eqwr,-aimag(eqgfret)/pi,append=TT)
    call eq_write_realt_functions
    call eq_write_realw_functions

    if(loop==nloop)then
       call eq_write_imagw_functions
       call eq_write_imagt_functions

       !Get nimp:
       call cfft_iw2tau(GMiw,Gt,beta)
       nimp=-2.d0*(real(Gt(L)))

       !Energy && k-dispersed quantities
       Ekin=0.0
       do ik=1,Lk
          GMiw=zero
          do i=1,2*L
             w=eqwm(i)
             iw=cmplx(0.d0,w)
             GMiw(i)=one/(iw - epsik(ik) - SbMiw(i) - SMiw(i))
          enddo
          call cfft_iw2tau(GMiw,Gt,beta)
          icnk(ik)=-Gt(L) !from G(-tau) = -G(beta-tau)==> G(tau=0-)=-G(beta)
          Ekin=Ekin + wt(ik)*icnk(ik)*epsik(ik)
       enddo
       Epot=dot_product(conjg(SMiw(:)),GMiw(:))/beta/2.0
       docc=0.5*nimp  - 0.25d0
       if(u/=0.0)docc = Epot/U + 0.5*nimp - 0.25d0
       call splot("EQSOLUTION/nkVSepsik.ipt",epsik(1:Lk),icnk(1:Lk))
       call splot("EQSOLUTION/EtotVSu.ipt",u,Ekin+Epot)
       call splot("EQSOLUTION/EkinVSu.ipt",u,Ekin)
       call splot("EQSOLUTION/EpotVSu.ipt",u,Epot)
       call splot("EQSOLUTION/doccVSu.ipt",u,docc)
    endif

  contains
    subroutine eq_write_realt_functions()
      !Real time  Functions
      !===========================================================
      call splot('EQSOLUTION/Sret_t.ipt',eqt,ex*eqstret,append=TT)
      call splot('EQSOLUTION/Sgtr_t.ipt',eqt,ex*eqstgtr,append=TT)
      call splot('EQSOLUTION/Sless_t.ipt',eqt,ex*eqstless,append=TT)
      call splot('EQSOLUTION/G0ret_t.ipt',eqt,ex*eqg0tret,append=TT)
      call splot('EQSOLUTION/G0gtr_t.ipt',eqt,ex*eqg0tgtr,append=TT)
      call splot('EQSOLUTION/G0less_t.ipt',eqt,ex*eqg0tless,append=TT)
      call splot('EQSOLUTION/Gret_t.ipt',eqt,ex*eqgtret,append=TT)
      call splot('EQSOLUTION/Ggtr_t.ipt',eqt,ex*eqgtgtr,append=TT)
      call splot('EQSOLUTION/Gless_t.ipt',eqt,ex*eqgtless,append=TT)
    end subroutine eq_write_realt_functions
    !------------------------!
    !------------------------!
    !------------------------!
    subroutine eq_write_realw_functions()
      !REAL freq. functions
      !===========================================================
      call splot('EQSOLUTION/Sret_realw.ipt',eqwr,eqsfret,append=TT)
      call splot('EQSOLUTION/G0ret_realw.ipt',eqwr,eqg0fret,append=TT)
      call splot('EQSOLUTION/G0less_realw.ipt',eqwr,eqg0fless,append=TT)
      call splot('EQSOLUTION/G0gtr_realw.ipt',eqwr,eqg0fgtr,append=TT)
      call splot('EQSOLUTION/Gret_realw.ipt',eqwr,eqgfret,append=TT)
      call splot('EQSOLUTION/Gless_realw.ipt',eqwr,eqgfless,append=TT)
      call splot('EQSOLUTION/Ggtr_realw.ipt',eqwr,eqgfgtr,append=TT)
    end subroutine eq_write_realw_functions
    !------------------------!
    !------------------------!
    !------------------------!
    subroutine eq_write_imagw_functions()
      !Get Matsubara Green's function
      !===========================================================
      call getGmats(eqwr,eqgfret,GMiw,beta)
      call getGmats(eqwr,eqg0fret,G0Miw,beta)
      call getGmats(eqwr,eqsfret,SMiw,beta)
      call splot('EQSOLUTION/GM_iw.ipt',eqwm,GMiw)
      call splot('EQSOLUTION/G0M_iw.ipt',eqwm,G0Miw)
      call splot('EQSOLUTION/SigmaM_iw.ipt',eqwm,SMiw)
    end subroutine eq_write_imagw_functions
    !------------------------!
    !------------------------!
    !------------------------!
    subroutine eq_write_imagt_functions()
      !Imaginary Time functions
      !===========================================================
      call cfft_iw2tau(GMiw,Gt,beta)
      call cfft_iw2tau(G0Miw,G0t,beta)
      call extract(Gt,Gtau);call extract(G0t,G0tau)
      call splot('EQSOLUTION/GM_tau.ipt',tau,Gtau)
      call splot('EQSOLUTION/G0M_tau.ipt',tau,G0tau)
    end subroutine eq_write_imagt_functions
  end subroutine eq_dump_out
  !******************************************************************
  !******************************************************************
  !******************************************************************




end module EQKELDYSH


