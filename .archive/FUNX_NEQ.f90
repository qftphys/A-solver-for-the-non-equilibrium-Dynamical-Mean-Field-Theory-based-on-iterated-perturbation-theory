module FUNX_NEQ
  !###############################################################
  !     PROGRAM  : FUNCS_GLOBAL
  !     TYPE     : Module
  !     PURPOSE  : Constructs some functions used in other places. 
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !###############################################################
  use VARS_GLOBAL
  use FUNX_GLOBAL
  use FFTW
  use MKL95_BLAS
  use MKL95_PRECISION

  implicit none
  private
  public At, &
       Bt, &
       Xc, &
       getG0guess, &
       getG0free, &
       getG0analytic, &
       dyson4G0ret, &
       update_G0gtr_less, &
       update_G0lceil, &
       update_sigma
  save
contains
  !+----------------------------------------------------------------+
  !PROGRAM  : At
  !TYPE     : function
  !PURPOSE  : get A=int_t'^t cos(E*tau)*dtau
  ! = sin(Et')-sin(Et)=2*sin(E(t'-t)/2)*cos(E(t+t')/2)
  !+----------------------------------------------------------------+
  function At(i,j)
    integer :: i,j
    real(8) :: At
    real(8) :: alfa,beta,trel,tavg
    At=0.d0
    alfa=Efield*t(j) !E*t`
    beta=Efield*t(i) !E*t 
    At=sin(beta)-sin(alfa)
    At=At/Efield
  end function At



  !+------------------------------------------------------------------+
  !PROGRAM  : Bt
  !TYPE     : function
  !PURPOSE  : get B=int_t'^t sin(E*tau)*dtau
  ! = cos(Et)-cos(Et') = 2
  !+------------------------------------------------------------------+
  function Bt(i,j)
    integer :: i,j
    real(8) :: Bt
    real(8) :: alfa,beta,trel,tavg
    Bt=0.d0
    alfa=Efield*t(j) !E*t'
    beta=Efield*t(i) !E*t
    Bt=cos(alfa)-cos(beta)
    Bt=Bt/Efield
  end function Bt


  !+------------------------------------------------------------------+
  !PROGRAM  : Xc
  !TYPE     : function
  !PURPOSE  : build  the function X^c(e,t,t') equal to the term in front
  !of the exponential in the expressions for G_0(t,t') 
  !COMMENT  : The function is defined WITHOUT the imaginary unit xi
  !+------------------------------------------------------------------+
  function Xc(cc,x)
    real(8) :: Xc,x
    character(len=1) :: cc
    Xc=0.d0
    if(cc == '<')then
       Xc=fermi(x)
    elseif(cc == '>')then
       Xc=(1.d0-fermi(x))
    endif
  end function Xc


  !+-------------------------------------------------------------------+
  !PROGRAM  : GETG0GUESS
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial fields G0^{<,>,R,\lceil} used to obtain
  !the self-energies functions.
  !+-------------------------------------------------------------------+
  subroutine getG0guess
    integer :: i,j,k,itau,jtau
    real(8) :: w
    complex(8) :: peso,A
    real(8),dimension(2*L) :: rho,nfp,nfm

    rho=aimag(eqG0w)/pi

    do i=1,2*L
       w=wmin+dble(i-1)*fmesh
       nfp(i)=-Xc('>',w)
       nfm(i)=Xc('<',w)
    enddo
    G0gtr  =zero
    G0less =zero
    G0ret  =zero
    G0lceil=zero
    G0rceil=zero
    do i=0,L
       do j=0,L
          do k=1,2*L
             w=wmin+dble(k-1)*fmesh
             peso=exp(-xi*w*At(i,j))
             peso=peso*exp(-ts**2*Bt(i,j)**2/4.d0)
             peso=peso*exp(xi*xmu*(t(i)-t(j)))
             G0gtr(i,j)=G0gtr(i,j)+xi*fmesh*rho(k)*nfp(k)*peso
             G0less(i,j)=G0less(i,j)+xi*fmesh*rho(k)*nfm(k)*peso
          enddo
          !See below: this maybe NOT general
          !G0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j)-G0less(i,j))
       enddo
       !do jtau=0,Ltau
       !G0lceil(i,jtau)=-xi*G0ret(i,0)*eqG00tau(Ltau-jtau)
       !G0rceil(jtau,i)=-xi*eqG00tau(jtau)*conjg(G0ret(i,0))
       !enddo
    enddo


    !Improved version of Gret: 
    G0gtr=conjg(G0less)
    do i=0,L
       do j=0,L
          G0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j)-G0less(i,j))
       enddo
       G0ret(i,i)=-xi
       do jtau=0,Ltau
          G0lceil(i,jtau)=-xi*G0ret(i,0)*eqG00tau(Ltau-jtau)
          G0rceil(jtau,i)=-xi*eqG00tau(jtau)*conjg(G0ret(i,0))
       enddo
    enddo
    call plot_dislin("G0less_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),G0less(0:nstep,0:nstep))
    call plot_dislin("G0gtr_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),G0gtr(0:nstep,0:nstep))
    call plot_dislin("G0ret_t1t2","X","Y","Z",  &
         t(0:nstep),t(0:nstep),G0ret(0:nstep,0:nstep))
    call plot_dislin("G0lceil_t1t2","X","Y","Z",  &
         t(0:nstep),tau,G0lceil(0:nstep,0:Ltau))
    call plot_dislin("G0rceil_t1t2","X","Y","Z",  &
         tau,t(0:nstep),G0rceil(0:Ltau,0:nstep))
    return
  end subroutine getG0guess





  !+-------------------------------------------------------------------+
  !PROGRAM  : GETG0FREE
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial fields G0^{<,>,R,\lceil} used to obtain
  !the self-energies functions.
  !+-------------------------------------------------------------------+
  subroutine getG0free
    integer :: i,j,ik,itau,jtau,istep,it
    real(8) :: w,en,intE
    complex(8) :: peso,A,nt
    real(8) :: rho,nfp,nfm
    complex(8),dimension(0:nstep,0:nstep) :: dummy,Unity,dummy1
    G0gtr  =zero
    G0less =zero
    G0ret  =zero
    G0lceil=zero
    G0rceil=zero
    do i=0,nstep
       do j=0,nstep
          do ik=1,Lk
             en=epsik(ik)
             nfp=-Xc('>',en)
             nfm=Xc('<',en)
             intE=int_Epskt(ik,i,j)
             peso=exp(-xi*intE)
             peso=peso*exp(xi*xmu*(t(i)-t(j)))
             G0gtr(i,j)=G0gtr(i,j)+xi*nfp*peso/dble(Lk)
             G0less(i,j)=G0less(i,j)+xi*nfm*peso/dble(Lk)
          enddo
          !G0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j)-G0less(i,j))
       enddo
    enddo
    forall(i=0:L) G0ret(i,i)=-xi
    !Improved version of Gret: 
    G0gtr=conjg(G0less)
    do i=0,nstep
       do j=0,nstep
          G0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j)-G0less(i,j))
       enddo
       G0ret(i,i)=-xi
       do jtau=0,Ltau
          G0lceil(i,jtau)=-xi*G0ret(i,0)*eqG00tau(Ltau-jtau)
          G0rceil(jtau,i)=-xi*eqG00tau(jtau)*conjg(G0ret(i,0))
       enddo
    enddo
    do i=0,nstep
       nt=-xi*G0less(i,i)
       write(50,*)t(i),real(nt),aimag(nt)
    enddo
    call plot_dislin("G0less_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),G0less(0:nstep,0:nstep))
    call plot_dislin("G0gtr_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),G0gtr(0:nstep,0:nstep))
    call plot_dislin("G0ret_t1t2","X","Y","Z",  &
         t(0:nstep),t(0:nstep),G0ret(0:nstep,0:nstep))
    call plot_dislin("G0lceil_t1t2","X","Y","Z",  &
         t(0:nstep),tau,G0lceil(0:nstep,0:Ltau))
    call plot_dislin("G0rceil_t1t2","X","Y","Z",  &
         tau,t(0:nstep),G0rceil(0:Ltau,0:nstep))



    ppG0gtr  =zero
    ppG0less =zero
    ppG0ret  =zero
    ppG0lceil=zero
    ppG0rceil=zero
    do i=0,nstep
       do j=0,nstep
          do ik=1,Lk
             en=epsik(ik)
             nfp=-Xc('>',en)
             nfm=Xc('<',en)
             peso=exp(-xi*en*(t(i)-t(j)))
             peso=peso*exp(xi*xmu*(t(i)-t(j)))
             ppG0gtr(i,j)=ppG0gtr(i,j)+xi*nfp*peso/dble(Lk)
             ppG0less(i,j)=ppG0less(i,j)+xi*nfm*peso/dble(Lk)
          enddo
       enddo
    enddo

    !forall(i=0:L) ppG0ret(i,i)=-xi
    !Improved version of Gret: 
    ppG0gtr=conjg(ppG0less)
    do i=0,nstep
       do j=0,nstep
          ppG0ret(i,j)=heaviside(t(i)-t(j))*(ppG0gtr(i,j)-ppG0less(i,j))
       enddo
       ppG0ret(i,i)=-xi
       do jtau=0,Ltau
          ppG0lceil(i,jtau)=-xi*ppG0ret(i,0)*eqG00tau(Ltau-jtau)
          ppG0rceil(jtau,i)=-xi*eqG00tau(jtau)*conjg(ppG0ret(i,0))
       enddo
    enddo
    do i=0,nstep
       nt=-xi*ppG0less(i,i)
       write(51,*)t(i),real(nt),aimag(nt)
    enddo

    call plot_dislin("ppG0less_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),ppG0less(0:nstep,0:nstep))
    call plot_dislin("ppG0gtr_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),ppG0gtr(0:nstep,0:nstep))
    call plot_dislin("ppG0ret_t1t2","X","Y","Z",  &
         t(0:nstep),t(0:nstep),ppG0ret(0:nstep,0:nstep))
    call plot_dislin("ppG0lceil_t1t2","X","Y","Z",  &
         t(0:nstep),tau,ppG0lceil(0:nstep,0:Ltau))
    call plot_dislin("ppG0rceil_t1t2","X","Y","Z",  &
         tau,t(0:nstep),ppG0rceil(0:Ltau,0:nstep))

    return
  end subroutine getG0free



  !+-------------------------------------------------------------------+
  !PROGRAM  : GETG0
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial fields G0^{<,>,R,\lceil} used to obtain
  !the self-energies functions.
  !+-------------------------------------------------------------------+
  subroutine getG0analytic
    integer :: i,j,ik,itau,jtau
    real(8) :: w,en,intEp,intEm,gk
    complex(8) :: peso,A,nt
    real(8) :: rho,nfp,nfm
    complex(8),dimension(Lk,0:nstep) :: Dk
    print*,""
    print*,"Entro getG0"
    call getEplus; call getEminus
    G0gtr  =zero
    G0less =zero
    G0ret  =zero
    do i =0,nstep
       do j =0,nstep
          do ik=1,Lk
             en=epsik(ik)
             nfp=-Xc('>',en)
             nfm=Xc('<',en)
             gk   = gammak(ik,i)
             intEp= int_Eplus(ik,i,j)
             intEm= int_Eminus(ik,i,j)
             peso=gk**2*exp(-xi*intEp) + exp(-xi*intEm)
             peso=peso/(1.d0+gk**2)
             G0gtr(i,j)=G0gtr(i,j)+xi*nfp*peso/dble(Lk)
             G0less(i,j)=G0less(i,j)+xi*nfm*peso/dble(Lk)
          enddo
          G0ret(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j)-G0less(i,j))
       enddo
    enddo
    forall(i=0:L) G0ret(i,i)=-xi
    do i=0,nstep
       nt=-xi*G0less(i,i)
       write(52,*)t(i),real(nt),aimag(nt)
    enddo
    call plot_dislin("analyticG0less_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),G0less(0:nstep,0:nstep))
    call plot_dislin("analyticG0gtr_t1t2","X","Y","Z", &
         t(0:nstep),t(0:nstep),G0gtr(0:nstep,0:nstep))
    call plot_dislin("analyticG0ret_t1t2","X","Y","Z",  &
         t(0:nstep),t(0:nstep),G0ret(0:nstep,0:nstep))
    return
  end subroutine getG0analytic








  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine dyson4G0ret()
    integer :: i,j
    complex(8),dimension(0:nstep,0:nstep) :: dummy
    dummy=locGret 
    call InvMat(dummy(0:nstep,0:nstep),nstep+1)
    call plot_dislin("locGinvret_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),dummy(0:nstep,0:nstep)) 
    !Dyson Equation for Ret. component \GG_0^R=[G^R^-1 - Sigma^R]^-1
    G0ret=zero
    G0ret(0:nstep,0:nstep)= dummy - Sret(0:nstep,0:nstep)    
    call plot_dislin("locG0invret_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),G0ret(0:nstep,0:nstep))
    !Invert G0ret(0:L,0:L)
    call InvMat(G0ret(0:nstep,0:nstep),nstep+1)
    call plot_dislin("locG0ret_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),G0ret(0:nstep,0:nstep))
    return
  end subroutine dyson4G0ret



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine update_G0gtr_less()
    integer :: i,j
    complex(8),dimension(0:nstep,0:nstep) :: G0adv
    complex(8),dimension(0:nstep,0:nstep) :: Unity,OBJret,OBJadv
    complex(8),dimension(0:nstep,0:nstep) :: dummy_locGret,dummy_locGadv
    complex(8),dimension(0:nstep,0:nstep) :: dummy_locGless
    complex(8),dimension(0:nstep,0:Ltau)  :: dummy_locGlceil
    complex(8),dimension(0:Ltau,0:nstep)  :: dummy_locGrceil
    complex(8),dimension(0:nstep,0:nstep) :: Op1,Op2,Op3,Op4,Op5
    include 'updateG0less.f90'
    call plot_dislin("locG0less_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),G0less(0:nstep,0:nstep))
    G0gtr=conjg((G0less))       
    call plot_dislin("locG0gtr_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),G0gtr(0:nstep,0:nstep))
    !Get the ADV. components.
    do i=0,nstep
       do j=0,nstep
          G0adv(i,j)=heaviside(t(i)-t(j))*(G0gtr(i,j)-G0less(i,j))
       enddo
    enddo
    call plot_dislin("locG0ret1_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),G0adv(0:nstep,0:nstep))
    return
  end subroutine update_G0gtr_less


  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine update_G0lceil()
    integer :: i,jtau,ndim
    do i=0,nstep
       do jtau=0,Ltau
          G0lceil(i,jtau)=-xi*G0ret(i,0)*G0lceil(0,jtau)
       enddo
    enddo
    G0rceil=conjg(transpose(G0lceil))
    call plot_dislin("locG0tau_t1t2","X","Y","Z",&
         t(0:nstep),tau,G0lceil(0:nstep,0:Ltau))
    return
  end subroutine update_G0lceil


  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine update_sigma
    integer :: i,j,itau,jtau,ndim
    character(len=64) :: cmd
    !Get Sigma^>,<
    do i=0,L
       do j=0,L
          Sgtr(i,j)=Vpd**2*ppG0gtr(i,j)+&
               (U**2)*(G0gtr(i,j)**2)*G0less(j,i)
          Sless(i,j)=Vpd**2*ppG0less(i,j)+&
               (U**2)*(G0less(i,j)**2)*G0gtr(j,i)
       enddo
    enddo
!    Sgtr(L,L)=Vpd**2*ppG0gtr(L,L)+&
!         (U**2)*(G0gtr(L,L)**2)*G0less(L-1,L-1)
!    Sless(L,L)=Vpd**2*ppG0less(L,L)+&
!         (U**2)*(G0less(L,L)**2)*G0gtr(L-1,L-1)
    !Get Sigma^t  (stored in Sigma^ret)
    do i=0,L
       do j=0,L
          Sret(i,j)=heaviside(t(i)-t(j))*Sgtr(i,j) +   &
               heaviside(t(j)-t(i))*Sless(i,j)
          !Real Sret:
          !Sret(i,j)=heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
       enddo
    enddo
    !Sigma^R = Sigma^t - Sigma^<
    Sret=Sret-Sless
    Sadv=transpose(conjg(Sret))
    !Get Sigma^\r/lceil
    do i=0,L
       do jtau=0,Ltau
          Slceil(i,jtau)=Vpd**2*ppG0lceil(i,jtau)+&
               (U**2)*(G0lceil(i,jtau)**2)*G0rceil(jtau,i)
          Srceil(jtau,i)=Vpd**2*ppG0rceil(jtau,i)+&
               (U**2)*(G0rceil(jtau,i)**2)*G0lceil(i,jtau)
       enddo
    enddo
    !Srceil=conjg(transpose(Slceil))
    Smatsubara=0.d0
    do itau=0,Ltau
       !do jtau=0,itau
       !Smatsubara(itau,jtau)=eqStau(itau-jtau)
       !Smatsubara(jtau,itau)=Smatsubara(itau,jtau)
       !enddo
       Smatsubara(itau,itau)=eqStau(itau)       
    enddo
    call plot_dislin("Sless_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),Sless(0:nstep,0:nstep))
    call plot_dislin("Sgtr_t1t2","X","Y","Z",&
         t(0:nstep),t(0:nstep),Sgtr(0:nstep,0:nstep))
    call plot_dislin("Sret_t1t2","X","Y","Z",& 
         t(0:nstep),t(0:nstep),Sret(0:nstep,0:nstep))
    call plot_dislin("Sadv_t1t2","X","Y","Z",& 
         t(0:nstep),t(0:nstep),Sadv(0:nstep,0:nstep))
    call plot_dislin("Slceil_t1t2","X","Y","Z",&
         t(0:nstep),tau,Slceil(0:nstep,0:Ltau))
    call plot_dislin("Srceil_t1t2","X","Y","Z",&
         tau,t(0:nstep),Srceil(0:Ltau,0:nstep))
    return
  end subroutine update_sigma






  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : function
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function deltak(ik,it)
    integer :: ik,it
    real(8) :: deltak
    deltak = epsik(ik) - epsikt(ik,it)
  end function deltak

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : function
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function gammak(ik,it)
    integer :: ik,it
    real(8) :: gammak
    gammak=2.d0*Vpd/(deltak(ik,it)+sqrt((deltak(ik,it))**2 + 4*Vpd**2))
  end function gammak

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : function
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine getEplus
    integer :: ik,it
    if(.not.allocated(Eplus))allocate(Eplus(Lk,0:L))
    do ik=1,Lk
       do it=0,L
          Eplus(ik,it)=epsikt(ik,it) + epsik(ik) + sqrt(deltak(ik,it)**2 + 4*Vpd**2)
          Eplus(ik,it)=Eplus(ik,it)/2.d0
       enddo
    enddo
  end subroutine GetEplus

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : function
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine getEminus
    integer :: ik,it
    if(.not.allocated(Eminus))allocate(Eminus(Lk,0:L))
    do ik=1,Lk
       do it=0,L
          Eminus(ik,it)=epsikt(ik,it) + epsik(ik) - sqrt(deltak(ik,it)**2 + 4*Vpd**2)
          Eminus(ik,it)=Eminus(ik,it)/2.d0
       enddo
    enddo
  end subroutine GetEminus

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : function
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function int_Epskt(ik,it,jt)
    integer :: i,ik,it,jt
    real(8) :: ti,tj,sgn
    real(8) :: int_Epskt,dummy
    ti=t(it)
    tj=t(jt)
    sgn=1.d0
    if(tj > ti)sgn=-1.d0    
    dummy=0.d0
    do i=jt,it,sgn!jt,it-1,sgn
       dummy=dummy + sgn*epsikt(ik,i)*dt
    enddo
    int_Epskt=dummy
    if(it==jt)int_Epskt=0.d0
  end function int_Epskt


  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : function
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function int_Eplus(ik,it,jt)
    integer :: i,ik,it,jt
    real(8) :: ti,tj,sgn
    real(8) :: int_Eplus,dummy
    ti=t(it)
    tj=t(jt)
    sgn=1.d0 
    if(tj > ti)sgn=-1.d0
    dummy=0.d0
    do i=jt,it,sgn
       dummy=dummy + sgn*Eplus(ik,i)*dt
    enddo
    int_Eplus=dummy   
    if(it==jt)int_Eplus=0.d0
  end function int_Eplus

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : function
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function int_Eminus(ik,it,jt)
    integer :: i,ik,it,jt
    real(8) :: ti,tj,sgn
    real(8) :: int_Eminus,dummy
    ti=t(it)
    tj=t(jt)
    sgn=1.d0 
    if(tj > ti)sgn=-1.d0
    dummy=0.d0
    do i=jt,it,sgn
       dummy=dummy + sgn*Eminus(ik,i)*dt
    enddo
    int_Eminus=dummy   
    if(it==jt)int_Eminus=0.d0
  end function int_Eminus

end module FUNX_NEQ
