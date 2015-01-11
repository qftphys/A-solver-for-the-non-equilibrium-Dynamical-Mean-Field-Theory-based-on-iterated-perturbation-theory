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
  !LOCAL:
  !===================
  USE VARS_GLOBAL
  implicit none

  !Greens functions
  !=========================================================
  integer                   :: i,j,ik,im
  complex(8),dimension(2*L) :: g0fret,gfret,sfret,sbfret
  complex(8),dimension(2*L) :: g0fless,g0fgtr
  complex(8),dimension(-L:L):: g0tgtr,g0tless,g0tret

  complex(8),dimension(2*L) :: gb0fless,gb0fgtr,sbfless,sbfgtr
  complex(8),dimension(-L:L):: gb0tgtr,gb0tless

  complex(8),dimension(-L:L):: gtgtr,gtless,gtret
  complex(8),dimension(-L:L):: stgtr,stless,stret
  complex(8),dimension(-L:L):: sbtgtr,sbtless,sbtret
  complex(8),dimension(2*L) :: Adummy,Abdummy,Sdummy,Sm
  complex(8),dimension(-L:L):: At,Abt
  real(8),dimension(0:L)    :: gt,gt0

  public keldysheq

contains 
  subroutine keldysheq(sigma,fg0)
    complex(8),dimension(2*L) :: sigma,fg0
    integer :: ieqloop
    !=========================================================

    write(*,*)"Get solution at EQUILIBRIUM "
    write(*,*)"----------------------------"

    !Init calculation
    !===============================================================
    call system("mkdir EQSOLUTION")   

    write(*,"(A,f11.2)")'    U = ',U
    write(*,"(A,f11.2)")'    V = ',Vpd
    write(*,"(A,f11.2)")' beta = ',beta
    write(*,"(A,f11.6)")' Mesh = ',fmesh

    !Set the contribution from the Bath:
    !===============================================================
    Sm=zero
    if(Vpd/=0.d0)call build_bath

    !Starts DMFT-loop
    g0fret   = zero; gfret    = zero ; sfret= zero
    do ieqloop=1,neqloop
       call print_bar(ieqloop,neqloop)

       !Get G_loc^R=\sum_\kaG(\ka,w) &&  G_loc^{<,>}  && Guess/Update G0^R
       !===============================================================
       select case(ieqloop)
       case (1)
          call guess_G0gl
       case default
          call update_G0gl
       end select

       !Impurity Solver: Sigma=U^2G0^3 && FFT t-->w
       !===============================================================
       call Simpurity

       call dump_out(ieqloop)
    enddo                     !here the dmft loops end up
    call close_bar
    !+======================================================-
    sigma=sfret
    fg0=g0fret
  end subroutine keldysheq
  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : GETBATH
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine build_bath
    integer    :: im
    real(8)    :: w,en,ngtr,nless,ex,A,An
    complex(8) :: peso,iw
    real(8),dimension(-L:L)   :: exa
    sbfret  =zero
    do i=1,2*L
       w=wr(i);iw=cmplx(w,eps)
       do im=1,Lmu
          sbfret(i)=sbfret(i) + one/(iw-epsimu(im))/dble(Lmu)
       enddo
    enddo
    do i=1,n
       w = wr(i)
       A = -aimag(sbfret(i))/pi
       An= A*fermi(w,beta)
       sbfless(i)= pi2*xi*An
       sbfgtr(i) = pi2*xi*(An-A)
    enddo
    sbfless= Vpd**2*sbfless
    sbfgtr = Vpd**2*sbfgtr
    sbfret = Vpd**2*sbfret

    Sm=zero
    do i=1,2*L
       w=wm(i);iw=cmplx(0.d0,w)
       do im=1,Lmu
          Sm(i)=Sm(i)+Vpd**2/(iw-epsimu(im))/dble(Lmu)
       enddo
    enddo
    ex=-1.d0       
    do i=-L,L
       ex=-ex
       exa(i)=ex
    enddo

  end subroutine Build_bath
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : GUESS_G0GL
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine guess_G0gl()
    integer    :: ik
    real(8)    :: en,ngtr,nless
    complex(8) :: peso
    g0tless=zero
    g0tgtr  =zero
    do ik=1,Lk
       en=epsik(ik)
       ngtr=fermi(en,beta)-1.d0
       nless=fermi(en,beta)
       do i=-L,L
          peso=exp(-xi*en*t(i))
          g0tless(i)=g0tless(i)+ xi*nless*peso*wt(ik)
          g0tgtr(i) =g0tgtr(i) + xi*ngtr*peso*wt(ik)
       enddo
    enddo
  end subroutine guess_G0gl
  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PROGRAM  : GETBATH
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine update_G0gl()
    complex(8) :: A,An,Sf
    real(8)    :: ex,w
    complex(8),dimension(-L:L)    :: gammatR
    complex(8)                :: zetan
    gfret=zero
    do i=1,n        
       zetan=wr(i)+xi*eps-sfret(i)-sbfret(i)
       do ik=1,Lk
          gfret(i)=gfret(i)+wt(ik)/(zetan-epsik(ik))
       enddo
    enddo
    g0fret=one/(one/gfret + sfret)
    do i=1,n
       w = wr(i)
       A = -aimag(g0fret(i))/pi
       An= A*fermi(w,beta)
       g0fless(i)= pi2*xi*An
       g0fgtr(i) = pi2*xi*(An-A)
    enddo
    call cfft_rw2rt(g0fless,g0tless,L)  ; g0tless=fmesh/pi2*g0tless
    call cfft_rw2rt(g0fgtr, g0tgtr,L)   ; g0tgtr =fmesh/pi2*g0tgtr
    do i=-L,L
       g0tret(i)=heaviside(t(i))*(g0tgtr(i) - g0tless(i))        
    enddo
    if(heaviside(0.d0)==1.d0)g0tret(0)=g0tret(0)/2.d0
  end subroutine update_G0gl
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : GETBATH
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine Simpurity
    do i=-L,L
       stgtr(i)=(U**2)*(g0tgtr(i)**2)*g0tless(-i) 
       stless(i)=(U**2)*(g0tless(i)**2)*g0tgtr(-i)
       stret(i)=heaviside(t(i))*(stgtr(i) - stless(i))
    enddo
    if(heaviside(0.d0)==1.d0)stret(0)=stret(0)/2.d0
    call cfft_rt2rw(stret,sfret,L) ;sfret=dt*sfret
  end subroutine Simpurity
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+----------------------------------------------------------------+
  !PROGRAM  : DUMP_OUT
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine dump_out(loop)
    real(8) :: time
    integer :: i,loop
    real(8) :: ex,nk(Lk),w
    real(8),dimension(-L:L)   :: exa
    complex(8),dimension(2*L) :: fg
    real(8),dimension(0:L)    :: tmpGtau
    real(8),dimension(0:Ltau) :: GtauR,G0tauR
    complex(8) :: iw

    call getGmats(wr,gfret,Adummy,beta);call cfft_iw2it(Adummy,gt,beta)
    open(100,file="EQSOLUTION/nVSiloop.ipt",access="append")
    write(100,*)dble(loop),-2.d0*(real(gt(L)))
    close(100)
    if(loop==neqloop)then
       call get_Ggtrless
       !Adv - real time  Functions
       !===========================================================
       ex=-1.d0       
       do i=-L,L
          ex=-ex
          exa(i)=ex
       enddo
       call splot('EQSOLUTION/Sigmaret_t.ipt',t,exa*aimag(stret),exa*real(stret))
       call splot('EQSOLUTION/Sigmagtr_t.ipt',t,exa*aimag(stgtr),exa*real(stgtr))
       call splot('EQSOLUTION/Sigmaless_t.ipt',t,exa*aimag(stless),exa*real(stless))
       call splot('EQSOLUTION/G0ret_t.ipt',t,exa*aimag(g0tret),exa*real(g0tret))
       call splot('EQSOLUTION/G0gtr_t.ipt',t,exa*aimag(g0tgtr),exa*real(g0tgtr))
       call splot('EQSOLUTION/G0less_t.ipt',t,exa*aimag(g0tless),exa*real(g0tless))
       call splot('EQSOLUTION/Gret_t.ipt',t,exa*aimag(gtret),exa*real(gtret))
       call splot('EQSOLUTION/Ggtr_t.ipt',t,exa*aimag(gtgtr),exa*real(gtgtr))
       call splot('EQSOLUTION/Gless_t.ipt',t,exa*aimag(gtless),exa*real(gtless))

       !Adv - real freq. functions
       !===========================================================
       call splot('EQSOLUTION/Gret_realw.ipt',wr,aimag(gfret),real(gfret))
       call splot('EQSOLUTION/Sigmaret_realw.ipt',wr,aimag(sfret),real(sfret))
       call splot('EQSOLUTION/G0ret_realw.ipt',wr,aimag(g0fret),real(g0fret))
       call splot('EQSOLUTION/DOS.ipt',wr,-aimag(gfret)/pi)

       !Get Matsubara Green's function
       !===========================================================
       call getGmats(wr,gfret,Adummy,beta)
       call getGmats(wr,g0fret,Abdummy,beta)
       call getGmats(wr,sfret,Sdummy,beta)
       call splot('EQSOLUTION/GM_iw.ipt',wm,aimag(Adummy),real(Adummy))
       call splot('EQSOLUTION/G0M_iw.ipt',wm,aimag(Abdummy),real(Abdummy))
       call splot('EQSOLUTION/SigmaM_iw.ipt',wm,aimag(Sdummy),real(Sdummy))
       do ik=1,Lk
          fg=zero
          do i=1,2*L
             w=wm(i)
             iw=cmplx(0.d0,w)
             fg(i)=one/(iw - epsik(ik) - Sm(i) - Sdummy(i))
          enddo
          call cfft_iw2it(fg,tmpGtau,beta) !\tau\in[-beta,0
          nk(ik)=-tmpGtau(L) !from G(-tau) = -G(beta-tau)==> G(tau=0-)=-G(beta)
       enddo
       call splot("EQSOLUTION/nkVSepsik.ipt",epsik(1:Lk),nk(1:Lk))


       !Imaginary Time functions
       !===========================================================
       call cfft_iw2it(Adummy,gt,beta)
       call cfft_iw2it(Abdummy,gt0,beta)
       call extract(gt,GtauR)
       call extract(gt0,G0tauR)
       call splot('EQSOLUTION/GM_tau.ipt',tau,GtauR)
       call splot('EQSOLUTION/G0M_tau.ipt',tau,G0tauR)
    endif
  end subroutine dump_out
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PROGRAM  : GETBATH
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine get_G0gtrless()
    complex(8),dimension(-L:L)    :: gammatR
    real(8)                       :: ex
    Adummy=g0fret
    call cfft_rw2rt(Adummy,g0tret,L) ; g0tret=fmesh/pi2*g0tret
    forall(i=-L:L)gammatR(i)=1.d0+gtret(i)*stret(i)*dt
    do i=-L,L
       write(101,*)t(i),aimag(g0tret(i))
       g0tless(i)=(1.d0/gammatR(i))*gtless(i)*(1.d0/conjg(gammatR(-i))) - g0tret(i)*stless(i)*conjg(g0tret(-i))
       g0tgtr(i)=(1.d0/gammatR(i))*gtgtr(i)*(1.d0/conjg(gammatR(-i)))  - g0tret(i)*stgtr(i)*conjg(g0tret(-i))
       !g0tret(i)=heaviside(t(i))*(g0tgtr(i) - g0tless(i))
    enddo
    write(101,*)""
  end subroutine get_G0gtrless
  !******************************************************************
  !******************************************************************
  !******************************************************************
  !+-------------------------------------------------------------------+
  !PROGRAM  : GETGTRLESS
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine get_Ggtrless()
    complex(8),dimension(-L:L)    :: gammatR
    real(8) :: ex
    !Get G^{<,>}=\gamma^R*G0^{<,>}*\gamma^A + G^R*Sigma^{<,>}*G^A
    !\gamma^R=[1-G^R*Sigma^R]
    !\gamma^A=[1-Sigma^A*G^A]=[\gamma^R]^+
    Adummy=gfret
    call cfft_rw2rt(Adummy,gtret,L) ; gtret=fmesh/pi2*gtret !;gtret(0)=-xi
    forall(i=-L:L)gammatR(i)=1.d0+gtret(i)*stret(i)*dt
    do i=-L,L
       gtless(i)=gammatR(i)*g0tless(i)*conjg(gammatR(-i)) + gtret(i)*stless(i)*conjg(gtret(-i))
       gtgtr(i)=gammatR(i)*g0tgtr(i)*conjg(gammatR(-i)) + gtret(i)*stgtr(i)*conjg(gtret(-i))
       gtret(i)=heaviside(t(i))*(gtgtr(i) - gtless(i))
    enddo
  end subroutine get_Ggtrless
  !******************************************************************
  !******************************************************************
  !******************************************************************






end module EQKELDYSH


