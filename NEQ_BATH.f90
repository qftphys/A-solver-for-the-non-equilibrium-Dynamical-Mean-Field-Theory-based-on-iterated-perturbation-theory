!###############################################################
!     PROGRAM  : BUILD THERMOSTATING BATHS
!     AUTHORS  : Adriano Amaricci
!###############################################################
module NEQ_BATH
  USE NEQ_VARS_GLOBAL
  implicit none
  private
  integer                          :: Lw
  real(8),allocatable,dimension(:) :: bath_dens,wfreq

  public                           :: get_thermostat_bath

contains
  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the Bath part of the system using exclusively time
  !dependent formalism. The bath is not interacting so it is 
  ! time translation invariant.
  !+-------------------------------------------------------------------+
  subroutine get_thermostat_bath()
    integer          :: iw,i,j
    real(8)          :: en,w,dw,wfin,wini
    complex(8)       :: peso
    real(8)          :: ngtr,nless,arg,lambda
    complex(8)       ::  S0gtr(Nstep,Nstep)
    ! complex(8)       :: s0less_w(Lw),s0less_t(-Lw/2:Lw/2)
    ! complex(8)       :: s0gtr_w(Lw),s0gtr_t(-Lw/2:Lw/2)
    call msg("Get Bath Type: "//reg(bath_type),id=0)
    call msg("Bath coupling is:"//txtfy(Vbath))
    call msg("Bath width is:"//txtfy(2.d0*Wbath))
    lambda=Vbath*2.d0*Wbath
    call msg("Bath coupling amplitude is:"//txtfy(lambda))
    Lw=L
    allocate(bath_dens(Lw),wfreq(Lw))
    wfin  = 2.d0*Wbath ; wini=-wfin
    wfreq = linspace(wini,wfin,Lw,mesh=dw)

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
       wfin  = 2.d0*Wbath+Wgap ; wini=-wfin
       wfreq = linspace(wini,wfin,Lw,mesh=dw)
       call get_bath_gapflat_dos()

    case default
       call abort("Bath type:"//reg(bath_type)//" not supported.")

    end select

    S0=zero
    S0gtr=zero
    if(Vbath/=0.d0)then
       do iw=1,Lw
          en   = wfreq(iw)
          nless= fermi(en,beta)
          ngtr = fermi(en,beta)-1.d0 !it absorbs the minus sign of the greater functions
          do i=1,nstep
             do j=1,nstep
                peso=exp(-xi*(time(i)-time(j))*en)
                S0%less(i,j)=S0%less(i,j)+ xi*lambda*nless*peso*bath_dens(iw)*dw
                S0gtr(i,j) =S0gtr(i,j)   + xi*lambda*ngtr*peso*bath_dens(iw)*dw
             enddo
          enddo
       enddo
       forall(i=1:Nstep,j=1:Nstep)S0%ret(i,j) = step(time(i)-time(j))*(S0gtr(i,j)-S0%less(i,j))
    endif

    if(mpiID==0)then
       call splot("DOSbath.lattice",wfreq,bath_dens,append=.true.)
       if(Vbath/=0.d0.AND.plot3D)call plot_keldysh_contour_gf(S0,time,reg(plot_dir)//"/S0")
    endif

  end subroutine get_thermostat_bath




  !+-----------------------------------------------------------------+
  !PURPOSE  : Build constant BATH 
  !+-----------------------------------------------------------------+
  subroutine get_bath_flat_dos()
    integer    :: i
    real(8)    :: w
    do i=1,Lw
       w=wfreq(i)
       bath_dens(i)= heaviside(Wbath-abs(w))/(2.d0*Wbath)
    enddo
  end subroutine get_bath_flat_dos

  subroutine get_bath_pgflat_dos()
    integer    :: i
    real(8)    :: w,norm
    norm=(Walpha+1.d0)/(2.d0*Wbath**(Walpha+1.d0))
    do i=1,Lw
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
    do i=1,Lw
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
    ! do i=-Lw,Lw
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
    do i=1,Lw
       w=wfreq(i)
       zeta=cmplx(w,eps,8)
       gf=gfbether(w,zeta,wbath/2.d0)
       bath_dens(i)=-aimag(gf)/pi
    enddo
  end subroutine get_bath_bethe_dos



  !********************************************************************
  !********************************************************************
  !********************************************************************



end module NEQ_BATH
