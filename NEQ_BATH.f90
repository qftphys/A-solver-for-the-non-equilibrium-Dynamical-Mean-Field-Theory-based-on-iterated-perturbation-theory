!###############################################################
!     PROGRAM  : BUILD THERMOSTATING BATHS
!     AUTHORS  : Adriano Amaricci
!###############################################################
module NEQ_BATH
  USE NEQ_VARS_GLOBAL
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
  subroutine get_thermostat_bath()
    integer          :: iw,i,j,ints
    real(8)          :: en,w,dw,wmax
    complex(8)       :: peso
    real(8)          :: ngtr,nless,arg,Vhopping
    complex(8)       ::  S0gtr(Nstep,Nstep)

    call msg("Get Bath Type: "//reg(bath_type),id=0)
    call msg("Bath coupling is       :"//txtfy(Vbath))
    Vhopping=sqrt(Vbath*2.d0*Wbath)
    call msg("Bath hopping, width are:"//reg(txtfy(Vhopping))//","//reg(txtfy(2.d0*Wbath)))

    allocate(bath_dens(L),wfreq(L))
    wmax  = Wbath
    wfreq = linspace(-wmax,wmax,L,mesh=dw)

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
       call abort("Bath type:"//reg(bath_type)//" not supported.")

    end select

    S0%less=zero
    S0%ret =zero
    S0gtr  =zero
    if(Vbath/=0.d0)then
       do i=1,nstep
          do j=1,i
             do iw=1,L
                en   = wfreq(iw)
                nless= fermi(en,beta)
                ngtr = fermi(en,beta)-1.d0 !it absorbs the minus sign of the greater functions
                peso=exp(-xi*(time(i)-time(j))*en)
                S0%less(i,j)=S0%less(i,j) + xi*Vhopping**2*nless*peso*bath_dens(iw)*dw
                S0gtr(i,j)  =S0gtr(i,j)   + xi*Vhopping**2*ngtr*peso*bath_dens(iw)*dw
             enddo
          enddo
       enddo

       S0%ret = S0gtr-S0%less
       forall(i=1:Nstep,j=1:Nstep,i<j)
          S0%less(i,j)=-conjg(S0%less(j,i))
          S0gtr(i,j)  =-conjg(S0gtr(j,i))
       end forall
       call splot("DOSbath.lattice",wfreq,bath_dens)
       if(plot3D)call plot_keldysh_contour_gf(S0,time,reg(plot_dir)//"/S0")

       !<<<DEBUG
       ints=200
       do j=1,Nstep,Nstep/10
          ints=ints+1
          rewind(ints)
          do i=1,Nstep
             write(ints,"(7F26.16)")time(i),dimag(S0%less(i,j)),dreal(S0%less(i,j)),&
                  dimag(S0gtr(i,j)),dreal(S0gtr(i,j)),&
                  dimag(S0%ret(i,j)),dreal(S0%ret(i,j))
          enddo
       enddo
       !>>>DEBUG
    endif



  end subroutine get_thermostat_bath




  !+-----------------------------------------------------------------+
  !PURPOSE  : Build constant BATH 
  !+-----------------------------------------------------------------+
  subroutine get_bath_flat_dos()
    integer    :: i
    real(8)    :: w
    do i=1,L
       w=wfreq(i)
       bath_dens(i)= step(Wbath-abs(w))/(2.d0*Wbath)
    enddo
  end subroutine get_bath_flat_dos

  subroutine get_bath_pgflat_dos()
    integer    :: i
    real(8)    :: w,norm
    norm=(Walpha+1.d0)/(2.d0*Wbath**(Walpha+1.d0))
    do i=1,L
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
    do i=1,L
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
    do i=1,L
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

