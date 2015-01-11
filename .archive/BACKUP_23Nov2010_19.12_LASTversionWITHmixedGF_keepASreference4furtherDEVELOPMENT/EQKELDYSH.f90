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
  !LIB:
  USE TOOLS
  !LOCAL:
  USE VARS_GLOBAL
  implicit none

  !Greens functions
  !=========================================================
  private
  integer                          :: i,j,ik,im
  real(8),dimension(:),allocatable :: ex
  real(8),allocatable,dimension(:) :: eqnk,eqnm
  complex(8),dimension(:),allocatable :: eqg0fret,eqg0fless,eqg0fgtr
  complex(8),dimension(:),allocatable :: eqgfret,eqgfless,eqgfgtr
  complex(8),dimension(:),allocatable :: eqsfret
  complex(8),dimension(:),allocatable :: eqsbfret
  !Time domain:
  complex(8),dimension(:),allocatable :: eqg0tret,eqg0tless,eqg0tgtr
  complex(8),dimension(:),allocatable :: eqgtret,eqgtless,eqgtgtr
  complex(8),dimension(:),allocatable :: eqstret,eqstless,eqstgtr
  complex(8),dimension(:),allocatable :: eqsbtret
  !Matsubara
  complex(8),dimension(:),allocatable :: eqfg,eqfg0,eqs,eqsb
  !Tau
  real(8),allocatable,dimension(:)    :: eqfgtau,eqfg0tau,eqstau
  real(8),allocatable,dimension(:,:)  :: eqGktau

  public :: equilibium_ipt_keldysh

contains 
  !+-------------------------------------------------------------------+
  !PROGRAM  : KELDYSHEQ
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine equilibium_ipt_keldysh(neqloop)
    integer :: ieqloop,neqloop

    write(*,*)"Get solution at EQUILIBRIUM "
    Write(*,*)"----------------------------"

    !Init calculation
    call system("mkdir EQSOLUTION")   
    call main_alloc("a")
    call eq_build_bath()

    !=======DMFT loop=====================
    call eq_guess_G0
    do ieqloop=1,neqloop
       call print_bar(ieqloop,neqloop)
       if(ieqloop>1)call eq_update_G0
       call eq_Simpurity
       call eq_get_Gloc
       call eq_dump_out(ieqloop,neqloop)
    enddo
    call close_bar

    !Store some arrays to provide initial conditions:
    icG0tgtr  = ex*eqg0tgtr
    icG0tless = ex*eqg0tless
    icG0w     = eqg0fret
    icSiw     = eqs
    icGiw     = eqfg
    icnk      = eqnk
    call extract(eqstau,icStau(0:Ltau))  ; forall(i=1:Ltau)icStau(-i)=-icStau(Ltau-i)
    call extract(eqfgtau,icGtau(0:Ltau)) ; forall(i=1:Ltau)icGtau(-i)=-icGtau(Ltau-i)
#ifdef _mix
    do ik=1,Lk
       call extract(eqGktau(ik,0:L),icGktau(ik,0:Ltau))
    enddo
#endif
    call main_alloc("d")
    return
  end subroutine equilibium_ipt_keldysh
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PROGRAM  : MAIN_ALLOC
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine main_alloc(char)
    integer          :: i,im
    character(len=1) :: char
    real(8)          :: exx,w
    complex(8)       :: iw,gf
    if(char=="a")then
       allocate(eqg0fret(2*L),eqg0fless(2*L),eqg0fgtr(2*L))
       allocate(eqgfret(2*L),eqgfless(2*L),eqgfgtr(2*L))
       allocate(eqsfret(2*L))
       allocate(eqg0tret(-L:L),eqg0tless(-L:L),eqg0tgtr(-L:L))
       allocate(eqgtret(-L:L),eqgtless(-L:L),eqgtgtr(-L:L))
       allocate(eqstret(-L:L),eqstless(-L:L),eqstgtr(-L:L))
       allocate(eqsbfret(2*L),eqsbtret(-L:L))
       allocate(eqfg(2*L),eqfg0(2*L),eqs(2*L),eqsb(2*L))
       allocate(eqfgtau(0:L),eqfg0tau(0:L),eqstau(0:L))
       allocate(eqnk(Lk),eqnm(Lmu),eqGktau(Lk,0:L))
       allocate(ex(-L:L))
       exx=-1.d0       
       do i=-L,L
          exx=-exx
          ex(i)=exx
       enddo
    elseif(char=="d")then
       deallocate(eqg0fret,eqg0fless,eqg0fgtr)
       deallocate(eqgfret,eqgfless,eqgfgtr)
       deallocate(eqsfret)
       deallocate(eqg0tret,eqg0tless,eqg0tgtr)
       deallocate(eqgtret,eqgtless,eqgtgtr)
       deallocate(eqstret,eqstless,eqstgtr)
       deallocate(eqsbfret,eqsbtret)
       deallocate(eqfg,eqfg0,eqs,eqsb)
       deallocate(eqfgtau,eqfg0tau,eqstau)
       deallocate(eqnk,eqnm,eqGktau)
       deallocate(ex)
    endif
    return
  end subroutine main_alloc
  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  : GETBATH
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine eq_build_bath
    integer    :: im,i
    real(8)    :: w,en,nless,ngtr
    complex(8) :: iw,ggtr(-L:L),gless(-L:L),peso,gf
    complex(8) :: funcM(2*L)
    real(8)    :: funcT(0:L)

    ! !Real time => Real freq.
    ! eqsbtret  =zero
    ! gless=zero;ggtr=zero
    ! do im=1,Lmu
    !    en=epsimu(im)
    !    ngtr=fermi0(en,beta)-1.d0
    !    nless=fermi0(en,beta)
    !    do i=-L,L
    !       peso=densmu(im)*exp(-xi*en*eqt(i))*de
    !       gless(i)=gless(i) + xi*nless*peso
    !       ggtr(i) =ggtr(i)  + xi*ngtr*peso
    !    enddo
    ! enddo
    ! forall(i=-L:L)eqsbtret(i)=heaviside(eqt(i))*(ggtr(i)-gless(i))
    ! eqsbtret = Vpd**2*eqsbtret*ex
    ! call cfft_rt2rw(eqsbtret,eqsbfret,L) ; eqsbfret=dt*eqsbfret
    ! call splot("bathSret_t_1.ipt",eqt,eqsbtret)
    ! call splot("bathSret_realw_1.ipt",eqwr,eqsbfret)

    eqsbfret=zero
    do i=1,2*L
       w=eqwr(i);iw=cmplx(w,eps)
       do im=1,Lmu
          eqsbfret(i)=eqsbfret(i) + de*densmu(im)/(iw-epsimu(im))!/dble(Lmu)
       enddo
    enddo
    eqsbfret = eqsbfret*Vpd**2
    call cfft_rw2rt(eqsbfret,eqsbtret,L)  ; eqsbtret=rfmesh/pi2*eqsbtret
    ! call splot("bathSret_t_2.ipt",eqt,eqsbtret)
    ! call splot("bathSret_realw_2.ipt",eqwr,eqsbfret)

    do i=1,2*L
       w=eqwm(i);iw=cmplx(0.d0,w)
       do im=1,Lmu
          eqsb(i)=eqsb(i)+de*densmu(im)/(iw-epsimu(im))
       enddo
    enddo
    eqsb=eqsb*Vpd**2


    do im=1,Lmu
       funcM=zero
       do i=1,2*L
          w=eqwm(i);iw=cmplx(0.d0,w)
          funcM(i)=one/(iw-epsimu(im))
       enddo
       call cfft_iw2tau(funcM,funcT,beta)
       eqnm(im)=-funcT(L) !from G(-tau) = -G(beta-tau)==> G(tau=0-)=-G(beta)
    enddo

    return
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
    complex(8) :: peso,zetan
    !Real time part:
    eqg0tless=zero ; eqg0tgtr  =zero
    do ik=1,Lk
       en=epsik(ik)
       ngtr=fermi0(en,beta)-1.d0
       nless=fermi0(en,beta)
       do i=-L,L
          peso=exp(-xi*en*eqt(i))
          eqg0tless(i)=eqg0tless(i) + xi*nless*peso*wt(ik)
          eqg0tgtr(i) =eqg0tgtr(i)  + xi*ngtr*peso*wt(ik)
       enddo
    enddo
    eqg0tless = ex*eqg0tless
    eqg0tgtr  = ex*eqg0tgtr
    forall(i=-L:L)eqg0tret(i)=heaviside(eqt(i))*(eqg0tgtr(i) - eqg0tless(i))
    call cfft_rt2rw(eqg0tret,eqg0fret,L) ;eqg0fret=dt*eqg0fret
    !Imaginary time part: (guess only metallic part)
    eqfg0=zero
    do i=1,2*L
       zetan=xi*wm(i) + xmu
       do ik=1,Lk
          eqfg0(i)=eqfg0(i)+wt(ik)/(zetan-epsik(ik))
       enddo
    enddo
    call cfft_iw2tau(eqfg0,eqfg0tau,beta)
    return
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
    !Real time part:
    forall(i=-L:L)
       eqstgtr(i)=(U**2)*(eqg0tgtr(i)**2)*eqg0tless(-i) 
       eqstless(i)=(U**2)*(eqg0tless(i)**2)*eqg0tgtr(-i)
       eqstret(i)=heaviside(eqt(i))*(eqstgtr(i) - eqstless(i))
    end forall
    if(heaviside(0.d0)==1.d0)eqstret(0)=eqstret(0)/2.d0
    call cfft_rt2rw(eqstret,eqsfret,L) ; eqsfret=dt*eqsfret
    !Imaginary time part:
    forall(i=0:L)eqstau(i)=U**2*(eqfg0tau(i))**2*eqfg0tau(L-i)
    call cfft_tau2iw(eqstau,eqs,beta)
    return
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
    !real time part:
    eqgfret=zero
    do i=1,2*L
       w=eqwr(i)
       zetan=cmplx(w,eps)-eqsfret(i)-eqsbfret(i)
       do ik=1,Lk
          eqgfret(i)=eqgfret(i)+wt(ik)/(zetan-epsik(ik))
       enddo
       A=-aimag(eqgfret(i))/pi
       eqgfless(i)= pi2*xi*fermi0(w,beta)*A
       eqgfgtr(i) = pi2*xi*(fermi0(w,beta)-1.d0)*A
    enddo
    call cfft_rw2rt(eqgfless,eqgtless,L)  ; eqgtless=rfmesh/pi2*eqgtless
    call cfft_rw2rt(eqgfgtr,eqgtgtr,L)    ; eqgtgtr=rfmesh/pi2*eqgtgtr
    forall(i=-L:L)eqgtret(i)=heaviside(eqt(i))*(eqgtgtr(i) - eqgtless(i))        
    if(heaviside(0.d0)==1.d0)eqgtret(0)=eqgtret(0)/2.d0 ;eqgtret(0)=-xi
    !imaginary time part:
    eqfg=zero
    do i=1,2*L
       zetan=xi*wm(i) + xmu - eqs(i) - eqsb(i)
       do ik=1,Lk
          eqfg(i) = eqfg(i) + wt(ik)/(zetan - epsik(ik))
       enddo
    enddo
    call cfft_iw2tau(eqfg,eqfgtau,beta)
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
    !Real time part:
    eqg0fret=one/(one/eqgfret + eqsfret)
    do i=1,2*L
       w = eqwr(i)
       A = -aimag(eqg0fret(i))/pi
       An= A*fermi0(w,beta)
       eqg0fless(i)= pi2*xi*An
       eqg0fgtr(i) = pi2*xi*(An-A)
    enddo
    call cfft_rw2rt(eqg0fless,eqg0tless,L)  ; eqg0tless=rfmesh/pi2*eqg0tless
    call cfft_rw2rt(eqg0fgtr, eqg0tgtr,L)   ; eqg0tgtr =rfmesh/pi2*eqg0tgtr
    forall(i=-L:L)eqg0tret(i)=heaviside(eqt(i))*(eqg0tgtr(i) - eqg0tless(i))        
    if(heaviside(0.d0)==1.d0)eqg0tret(0)=eqg0tret(0)/2.d0 ; eqg0tret(0)=-xi
    !Imaginary time part:
    eqfg0 = one/(one/eqfg + eqs) ; call cfft_iw2tau(eqfg0,eqfg0tau,beta)
    return
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
    complex(8) :: funcM(2*L)
    real(8)    :: funcT(0:L)

    call splot('EQSOLUTION/DOS.ipt',eqwr,-aimag(eqgfret)/pi,append=TT)
    call eq_write_realt_functions
    call eq_write_realw_functions

    if(loop==nloop)then
       call eq_write_imagw_functions
       call eq_write_imagt_functions

       !Get nimp:
       call cfft_iw2tau(eqfg,eqfgtau,beta)
       nimp=-2.d0*(real(eqfgtau(L)))

       !Energy && k-dispersed quantities
       Ekin=0.0
       do ik=1,Lk
          funcM=zero
          do i=1,2*L
             w=eqwm(i)
             iw=cmplx(0.d0,w)
             funcM(i)=one/(iw - epsik(ik) - eqs(i) - eqsb(i))
          enddo
          call cfft_iw2tau(funcM,funcT,beta)
          eqnk(ik)=-funcT(L) !from G(-tau) = -G(beta-tau)==> G(tau=0-)=-G(beta)
          eqGktau(ik,0:L) = funcT(0:L)
          Ekin=Ekin + wt(ik)*eqnk(ik)*epsik(ik)
       enddo
       Epot=dot_product(conjg(eqs(:)),eqfg(:))/beta/2.0
       docc=0.5*nimp  - 0.25d0
       if(u/=0.0)docc = Epot/U + 0.5*nimp - 0.25d0
       call splot("EQSOLUTION/nkVSepsik.ipt",epsik(1:Lk),eqnk(1:Lk))
       call splot("EQSOLUTION/nmuVSepsimu.ipt",epsimu(1:Lmu),eqnm(1:Lmu))
       call splot("EQSOLUTION/EtotVSu.ipt",u,Ekin+Epot)
       call splot("EQSOLUTION/EkinVSu.ipt",u,Ekin)
       call splot("EQSOLUTION/EpotVSu.ipt",u,Epot)
       call splot("EQSOLUTION/doccVSu.ipt",u,docc)
    endif

  contains
    subroutine eq_write_realt_functions()
      !Real time  Functions
      !===========================================================
      call splot('EQSOLUTION/bathSret_t.ipt',eqt,ex*eqsbtret,append=TT)
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
      call splot('EQSOLUTION/bathSret_realq.ipt',eqwr,eqsbfret,append=TT)
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
      call splot('EQSOLUTION/GM_iw.ipt',eqwm,eqfg)
      call splot('EQSOLUTION/G0M_iw.ipt',eqwm,eqfg0)
      call splot('EQSOLUTION/SM_iw.ipt',eqwm,eqs)
      call splot('EQSOLUTION/bathSM_iw.ipt',eqwm,eqsb)
    end subroutine eq_write_imagw_functions
    !------------------------!
    !------------------------!
    !------------------------!
    subroutine eq_write_imagt_functions()
      !Imaginary Time functions
      !===========================================================
      call splot('EQSOLUTION/GM_tau.ipt',eqtau,eqfgtau)
      call splot('EQSOLUTION/G0M_tau.ipt',eqtau,eqfg0tau)
    end subroutine eq_write_imagt_functions
  end subroutine eq_dump_out
  !******************************************************************
  !******************************************************************
  !******************************************************************




end module EQKELDYSH


