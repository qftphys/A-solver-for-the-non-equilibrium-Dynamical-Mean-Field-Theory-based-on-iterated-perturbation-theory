!###############################################################
!     PROGRAM  : BUILD THERMOSTATING BATHS
!     AUTHORS  : Adriano Amaricci
!###############################################################
module NEQ_THERMOSTAT
  USE NEQ_VARS_GLOBAL
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE ARRAYS
  USE FUNCTIONS
  USE IOTOOLS
  implicit none
  private
  real(8),allocatable,dimension(:) :: bath_dens,wfreq

  public                           :: get_thermostat_bath

contains
  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the Bath part of the system using exclusively time
  !dependent formalism. The bath is not interacting so it is 
  ! time translation invariant.
  !+-------------------------------------------------------------------+
  subroutine get_thermostat_bath(S0,params)
    type(kb_contour_gf)                   :: S0
    type(kb_contour_params)               :: params
    integer                               :: N,L
    integer                               :: iw,i,j,k,ints
    real(8)                               :: en,w,dw,wmax
    real(8)                               :: ngtr,nless,arg,Vhopping
    complex(8)                            :: peso,sless,sgtr    
    complex(8),dimension(:,:),allocatable :: S0gtr
    logical                               :: check

    if(.not.params%status)stop "thermostat/get_thermostat_bath: params not allocated. "
    if(.not.S0%status)stop "thermostat/get_thermostat_bath: S0 not allocated. "
    N = params%Ntime      !<== work with the MAX size of the contour
    L = params%Ntau

    write(*,"(A)")"Bath coupling is:"//txtfy(Vbath)
    Vhopping=sqrt(Vbath*2.d0*Wbath)
    write(*,"(A)")"Bath hopping, width are:"//reg(txtfy(Vhopping))//","//reg(txtfy(2.d0*Wbath))

    allocate(bath_dens(Lfreq),wfreq(Lfreq))
    wmax  = Wbath
    wfreq = linspace(-wmax,wmax,Lfreq,mesh=dw)
    select case(reg(bath_type))
    case("bethe")
       call get_bath_bethe_dos()
    case("gaussian")
       call get_bath_gaussian_dos()
    case ("flat")
       call get_bath_flat_dos()
    case("pgflat")
       call get_bath_pgflat_dos()
    case("gapflat")
       wmax  = Wbath+Wgap       !2.d0*Wbath+Wgap
       wfreq = linspace(-wmax,wmax,L,mesh=dw)
       call get_bath_gapflat_dos()
    case default
       stop "Bath type: not supported."
    end select
    call splot("DOSbath.plot",wfreq,bath_dens,append=.true.)


    !Bath self-energies:
    check = check_dimension_kb_contour(S0,N,L)
    allocate(S0gtr(N,N))
    S0   = zero
    S0gtr= zero
    if(Vbath==0.d0)return
    do iw=1,Lfreq
       en   = wfreq(iw)
       nless= fermi(en,beta)
       ngtr = fermi(en,beta)-1.d0 !it absorbs the minus sign of the greater functions
       do i=1,N
          !Less component:
          do j=1,i
             peso=exp(-xi*(params%t(i)-params%t(j))*en)
             S0%less(i,j)= S0%less(i,j) + xi*Vhopping**2*nless*peso*bath_dens(iw)*dw
             S0gtr(i,j) = S0gtr(i,j)    + xi*Vhopping**2*ngtr*peso*bath_dens(iw)*dw
          enddo
          !Lmix component:
          if(en>=0.d0)then
             do j=0,L
                if(beta*en>20.d0)then
                   peso=exp(-xi*en*params%t(i))*exp(-en*(beta-params%tau(j)))
                else
                   peso=exp(-xi*en*params%t(i))*exp(-en*(beta-params%tau(j)))/(1.d0+exp(-en*beta))
                endif
                S0%lmix(i,j) = S0%lmix(i,j) + xi*Vhopping**2*peso*bath_dens(iw)*dw
             enddo
          else
             do j=0,L
                if(beta*en<-20.d0)then
                   peso=exp(-xi*en*params%t(i))*exp(en*params%tau(j))
                else
                   peso=exp(-xi*en*params%t(i))*exp(-en*(beta-params%tau(j)))/(1.d0+exp(-en*beta))
                endif
                S0%lmix(i,j) = S0%lmix(i,j) + xi*Vhopping**2*peso*bath_dens(iw)*dw
             enddo
          endif
          !
       enddo
    enddo
    !Ret component:
    S0%ret = S0gtr-S0%less
    forall(i=1:params%Ntime,j=1:params%Ntime,i<j)S0%less(i,j)=-conjg(S0%less(j,i))
    call plot_kb_contour_gf("Sbath",S0,params)
  end subroutine get_thermostat_bath







  !+-----------------------------------------------------------------+
  !PURPOSE  : Build constant BATH 
  !+-----------------------------------------------------------------+
  subroutine get_bath_flat_dos()
    integer    :: i
    real(8)    :: w
    do i=1,Lfreq
       w=wfreq(i)
       bath_dens(i)= step(Wbath-abs(w))/(2.d0*Wbath)
    enddo
  end subroutine get_bath_flat_dos

  subroutine get_bath_pgflat_dos()
    integer    :: i
    real(8)    :: w,norm
    norm=(Walpha+1.d0)/(2.d0*Wbath**(Walpha+1.d0))
    do i=1,Lfreq
       w=wfreq(i)
       if(abs(w)>wbath)then
          bath_dens(i)=0.d0
       else
          bath_dens(i)=norm*abs(w)**Walpha
       end if
    end do
  end subroutine get_bath_pgflat_dos

  subroutine get_bath_gapflat_dos()
    integer    :: i
    real(8)    :: w,rho
    rho=1.d0/2.d0/(Wbath)!-Wgap)
    do i=1,Lfreq
       w=wfreq(i)
       if(abs(w)<Wbath+Wgap.AND.abs(w)>Wgap)then
          bath_dens(i)=rho
       else
          bath_dens(i)=0.d0         
       end if
    end do
  end subroutine get_bath_gapflat_dos

  subroutine get_bath_gaussian_dos()
    integer    :: i,ik
    real(8)    :: w,sig,alpha
    complex(8) :: gf,zeta
    bath_dens = exp(-0.5d0*(wfreq/Wbath)**2)/(sqrt(pi2)*Wbath) !standard Gaussian
    !bath_dens = exp(-((wfreq)/Wbath)**2)/(sqrt(pi)*Wbath) !Camille's choice
    !    !!w/ erf in frquency space: coded from Abramowitz-Stegun
    ! do i=-L,L
    !    !w=wfreq(i)
    !    !zeta=cmplx(w,eps,8)
    !    !sig=aimag(zeta)/abs(dimag(zeta))
    !    !gf=-sig*xi*sqrt(pi)*wfun(zeta/Wbath)/Wbath
    !    !bath_dens(i)=-aimag(gf)/pi
    ! enddo
  end subroutine get_bath_gaussian_dos

  subroutine get_bath_bethe_dos()
    integer    :: i,ik
    real(8)    :: w,sig,alpha
    complex(8) :: gf,zeta
    do i=1,Lfreq
       w=wfreq(i)
       zeta=cmplx(w,eps,8)
       gf=gfbether(w,zeta,wbath/2.d0)
       bath_dens(i)=-aimag(gf)/pi
    enddo
  end subroutine get_bath_bethe_dos



  !********************************************************************
  !********************************************************************
  !********************************************************************



end module NEQ_THERMOSTAT

