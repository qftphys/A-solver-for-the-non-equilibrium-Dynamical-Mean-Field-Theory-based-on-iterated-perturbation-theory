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
  private

  integer                               :: i,j,ik
  real(8)                               :: w
  complex(8)                            :: zeta
  real(8),dimension(:),allocatable      :: exx

  type(keldysh_equilibrium_gf)          :: fgk0,fgk,sigmak

  type(matsubara_gf)                    :: fg,fg0,sigma

  real(8),allocatable,dimension(:,:)    :: eqGktau
  real(8),allocatable,dimension(:)      :: eqnk

  public :: equilibium_ipt_keldysh

contains 

  !--------------------------------------------------------!
  subroutine equilibium_ipt_keldysh(neqloop)
    integer :: ieqloop,neqloop
    write(*,*)"Get solution at EQUILIBRIUM "
    Write(*,*)"----------------------------"

    call system("mkdir EQSOLUTION")   
    call main_alloc("a")


    sigmak=zero
    do iloop=1,neqloop
       write(*,"(A,i5,A,i5)")"DMFT-loop",iloop,"/",neqloop
       call get_gloc()
       call update_g0

       call simpurity
       call print_out(iloop)
    enddo

    !Store some arrays to provide initial conditions:
    icG0w     = fgk%ret%w       !conjg(fg0)
    icnk      = eqnk

    call extract(sigma%tau,icStau(0:Ltau)); forall(i=1:Ltau)icStau(-i)=-icStau(Ltau-i)
    call extract(fg%tau,icGtau(0:Ltau)); forall(i=1:Ltau)icGtau(-i)=-icGtau(Ltau-i)
    call main_alloc("d")
    call dump("")
    return
  end subroutine equilibium_ipt_keldysh
  !--------------------------------------------------------!


  !--------------------------------------------------------!
  subroutine main_alloc(char)
    integer          :: i,im
    character(len=1) :: char
    real(8)          :: fxx,w
    complex(8)       :: iw,gf
    if(char=="a")then
       call allocate_gf(fg,L)
       call allocate_gf(fg0,L)
       call allocate_gf(sigma,L)
       call allocate_gf(fgk,L)
       call allocate_gf(fgk0,L)
       call allocate_gf(sigmak,L)
       allocate(eqnk(Lk),eqGktau(Lk,0:L))
       allocate(exx(-L:L))
       fxx=-1.d0       
       do i=-L,L
          fxx=-fxx
          exx(i)=fxx
       enddo
    elseif(char=="d")then
       deallocate(eqnk,eqGktau)
       deallocate(exx)
    endif
    return
  end subroutine main_alloc
  !--------------------------------------------------------!

  ! !--------------------------------------------------------!
  ! subroutine guess_g0
  !   sigmak=zero
  !   call get_gloc
  !   fgk0%ret%w=fgk%ret%w
  ! end subroutine guess_g0
  ! !--------------------------------------------------------!

  !--------------------------------------------------------!
  subroutine get_gloc
    real(8)    :: w
    complex(8) :: ifg
    fgk%ret%w=zero
    do i=-L,L
       w = eqwr(i);zeta=cmplx(w,eps)-sigmak%ret%w(i)
       ifg=zero
       do ik=1,Lk
          ifg=ifg+wt(ik)/(zeta-epsik(ik))
       enddo
       fgk%ret%w(i)=ifg
    enddo
  end subroutine get_gloc
  !--------------------------------------------------------!

  !--------------------------------------------------------!
  subroutine update_g0
    complex(8),dimension(-L:L)          :: gf
    gf = one/(one/fgk%ret%w + sigmak%ret%w)

    !Update Keldysh components of the Weiss Field:
    fgk0%less%w=less_component_w(gf,eqwr,beta)
    fgk0%gtr%w =gtr_component_w(gf,eqwr,beta)

    !FFT to real time:
    !G0^{<,>}(t) = FT[G0^{<,>}(w)]
    call cfft_rw2rt(fgk0%less%w,fgk0%less%t,L);   fgk0%less%t=xi*rfmesh*fgk0%less%t
    call cfft_rw2rt(fgk0%gtr%w ,fgk0%gtr%t ,L);   fgk0%gtr%t =xi*rfmesh*fgk0%gtr%t
    fgk0%ret%t=ret_component_t(fgk0%gtr%t,fgk0%less%t,t,L)

    return
  end subroutine update_g0
  !--------------------------------------------------------!


  !--------------------------------------------------------!
  subroutine simpurity
    do i=-L,L
       sigmak%less%t(i)=(U**2)*(fgk0%less%t(i)**2)*fgk0%gtr%t(-i)
       sigmak%gtr%t(i) =(U**2)*(fgk0%gtr%t(i)**2)*fgk0%less%t(-i)
    enddo
    sigmak%ret%t=ret_component_t(sigmak%gtr%t,sigmak%less%t,t,L)
    call cfft_rt2rw(sigmak%ret%t,sigmak%ret%w,L); sigmak%ret%w= dt*sigmak%ret%w
    call cfft_rt2rw(fgk0%ret%t,fgk0%ret%w,L)    ; fgk0%ret%w  = dt*fgk0%ret%w

    return
  end subroutine simpurity
  !--------------------------------------------------------!


  !--------------------------------------------------------!
  subroutine print_out(loop)
    integer    :: loop
    real(8)    :: w,Ekin,Epot,docc,nimp
    type(matsubara_gf) :: func

    call getGmats(eqwr,fgk%ret%w,fg%iw,beta)
    call cfft_iw2tau(fg%iw,fg%tau,beta)
    nimp=-2.d0*fg%tau(L)
    call splot("EQSOLUTION/nVSiloop.ipt",iloop,nimp,append=TT)

    if(loop==eqnloop)then
       call splot('EQSOLUTION/Sret_t.ipt',eqt,exx*sigmak%ret%t,append=TT)
       call splot('EQSOLUTION/Sless_t.ipt',eqt,exx*sigmak%less%t,append=TT)
       call splot('EQSOLUTION/Sgtr_t.ipt',eqt,exx*sigmak%gtr%t,append=TT)

       call splot('EQSOLUTION/G0ret_t.ipt',eqt,exx*fgk0%ret%t,append=TT)
       call splot('EQSOLUTION/G0gtr_t.ipt',eqt,exx*fgk0%gtr%t,append=TT)
       call splot('EQSOLUTION/G0less_t.ipt',eqt,exx*fgk0%less%t,append=TT)

       call splot('EQSOLUTION/Gret_realw.ipt',eqwr,fgk%ret%w,append=TT)
       call splot('EQSOLUTION/Sret_realw.ipt',eqwr,sigmak%ret%w,append=TT)
       call splot('EQSOLUTION/G0ret_realw.ipt',eqwr,fgk0%ret%w,append=TT)

       call splot('EQSOLUTION/DOS.ipt',eqwr,-aimag(fgk%ret%w)/pi,append=TT)

       ! fgk%ret%w = one/(one/fgk0%ret%w - sigmak%ret%w)
       ! call splot('EQSOLUTION/impGret_realw.ipt',eqwr,fgk%ret%w,append=TT)
       ! call splot('EQSOLUTION/impDOS.ipt',eqwr,aimag(fgk%ret%w)/pi,append=TT)

       call getGmats(eqwr,fgk%ret%w,fg%iw,beta)
       call getGmats(eqwr,fgk0%ret%w,fg0%iw,beta)
       call getGmats(eqwr,sigmak%ret%w,sigma%iw,beta)
       call splot('EQSOLUTION/GM_iw.ipt',eqwm,fg%iw)
       call splot('EQSOLUTION/G0M_iw.ipt',eqwm,fg0%iw)
       call splot('EQSOLUTION/SigmaM_iw.ipt',eqwm,sigma%iw)
       call splot('EQSOLUTION/GM_tau.ipt',eqtau,fg%tau)

       !Energy && k-dispersed quantities
       call allocate_gf(func,L)
       Ekin=0.0
       do ik=1,Lk
          func%iw=zero
          do i=1,L
             w=eqwm(i)
             func%iw(i)=one/(xi*w - epsik(ik) - sigma%iw(i))
          enddo
          call cfft_iw2tau(func%iw,func%tau,beta)
          eqnk(ik)=-func%tau(L)
          eqGktau(ik,0:L) = func%tau(0:L)
          Ekin=Ekin + wt(ik)*eqnk(ik)*epsik(ik)
       enddo
       Epot=dot_product(conjg(sigma%iw(:)),fg%iw(:))/beta/2.0
       docc=0.5*nimp  - 0.25d0
       if(u/=0.0)docc = Epot/U + 0.5*nimp - 0.25d0
       call splot("EQSOLUTION/nkVSepsik.ipt",epsik(1:Lk),eqnk(1:Lk))
       call splot("EQSOLUTION/EtotVSu.ipt",u,Ekin+Epot)
       call splot("EQSOLUTION/EkinVSu.ipt",u,Ekin)
       call splot("EQSOLUTION/EpotVSu.ipt",u,Epot)
       call splot("EQSOLUTION/doccVSu.ipt",u,docc)
       print*,nimp
    endif
    return
  end subroutine print_out
  !--------------------------------------------------------!

end module EQKELDYSH


