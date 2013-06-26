!###############################################################
!     PROGRAM  : BATH
!     AUTHORS  : Adriano Amaricci
!     LAST UPDATE: 07/2009
!###############################################################
module BATH
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
    integer          :: iw,i,j,k
    real(8)          :: en,w,dw,wfin,wini
    complex(8)       :: peso
    real(8)          :: ngtr,nless,arg,lambda
    call msg("Get Bath. Type: "//bold(reg(bath_type))//" dissipative bath")
    call msg("Using bath coupling is:"//txtfy(Vbath))
    lambda=Vbath*2.d0*Wbath
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
       call error("Bath type:"//reg(bath_type)//" not supported.")
    end select

    S0=zero
    if(Vbath/=0.d0)then
       do iw=1,Lw
          en   = wfreq(iw)
          nless= fermi(en,beta)
          ngtr = fermi(en,beta)-1.d0!it absorbs the minus sign of the greater functions
          do i=1,Nfit
             do j=1,Nfit
                peso=exp(-xi*(tfit(i)-tfit(j))*en)
                k = pack_index(i,j,Nfit)
                S0%less(k)=S0%less(k)+ xi*lambda*nless*peso*bath_dens(iw)*dw
                S0%gtr(k) =S0%gtr(k) + xi*lambda*ngtr*peso*bath_dens(iw)*dw
             enddo
          enddo
       enddo
    endif

    call splot("DOSbath.neq",wfreq,bath_dens,append=.true.)
    if(Vbath/=0.d0.AND.plot3D)call plot_keldysh_contour_gf(tfit,S0,"S0bath")
  end subroutine get_thermostat_bath


  !********************************************************************
  !********************************************************************
  !********************************************************************


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



end module BATH

