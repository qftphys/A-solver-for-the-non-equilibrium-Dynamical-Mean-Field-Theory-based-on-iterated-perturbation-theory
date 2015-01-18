!#####################################################################
!     PURPOSE  : Function for the Electric field vector potential
!     AUTHORS  : Adriano Amaricci
!#####################################################################
MODULE ELECTRIC_FIELD
  USE NEQ_CONTOUR
  USE NEQ_VARS_GLOBAL
  USE CONSTANTS
  USE FUNCTIONS
  USE DERIVATE
  USE IOTOOLS
  implicit none
  private
  public :: Afield
  public :: set_efield_vector

  !ELECTRIC FIELD
  !=========================================================  
  real(8),dimension(3) :: Ek         !Electric field vector potential and vector

contains

  !+------------------------------------------------------------+
  !PURPOSE: set the normalized electric field versors using given direction
  !+------------------------------------------------------------+
  subroutine set_efield_vector(time)
    real(8),dimension(:) :: time
    real(8)              :: modulo
    integer              :: i
    logical              :: check
    !Normalize the Electric Field components
    !Keep unaltered the Electric Field Strenght Efield=E0
    modulo=sqrt(Ex**2+Ey**2)
    if(modulo/=0.d0)then
       Ex=Ex/modulo
       Ey=Ey/modulo
    endif
    Ek = [Ex,Ey,0d0]
    print*,"|E|=E0="//trim(txtfy(Efield/modulo))
    check=.false.
    check=field_type=="dc".OR.&
         field_type=="ac".OR.&
         field_type=="acdc".OR.&
         field_type=="pulse".OR.&
         field_type=="ramp"
    if(.not.check)stop "ELECTRIC_FIELD/Afield: wrong field_type. set:dc,ac,acdc,pulse,ramp"
    call print_field(time)
  end subroutine set_efield_vector


  subroutine print_field(t)
    real(8),dimension(:) :: t
    integer              :: i
    real(8),dimension(3) :: A
    real(8),dimension(size(t)) :: Ax,Ay,Ex,Ey
    do i=1,size(t)
       A=Afield(t(i))
       Ax(i)=A(1)
       Ay(i)=A(2)
    enddo
    Ex = deriv(Ax,dt)
    Ey = deriv(Ay,dt)
    open(10,file="Avector_shape.ipt")
    open(11,file="Efield_shape.ipt")
    do i=1,size(t)
       write(10,*)t(i),Afield(t(i))
       write(11,*)t(i),Ex(i),Ey(i)
    enddo
    close(10)
    close(11)
    if(field_type=="ac")print*,"Root condition: "//trim(txtfy(bessel_j0(Efield/Omega0)))
  end subroutine print_field



  !+------------------------------------------------------------+
  !PURPOSE : 
  !+------------------------------------------------------------+
  function Afield(t)
    real(8),intent(in)      :: t
    real(8)                 :: ftime,tau0,tau1
    real(8),dimension(3)    :: Afield
    complex(8)              :: zp,zm
    select case(field_type)
    case ("dc")                !DC ELECTRIC FIELD:
       ftime=-(step(t-t0)*(t-t0 + (t1-t)*step(t-t1) - (t1-t0)*step(t0-t1)))
       Afield=Ek*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("ac")                  !AC ELECTRIC FIELD
       ftime=-sin(Omega0*(t-t0))/Omega0
       Afield=Ek*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("acdc")                !AC+DC ELECTRIC FIELD (super-bloch)
       !ftime=-(t+sin(Omega0*(t-t0))/Omega0)
       ftime =-sin(Omega0*(t-t0))/Omega0
       Afield=Ek*(Efield*ftime - E1*t)       !A(t) = E0*F(t)*(e_x + e_y)

    case("pulse")               !LIGHT PULSE (for Pump&Probe) 
       !Signal function:
       !cos(\Omega*(t-t0)-pi/2)Exp(-((t-t0)/tau0)**2)
       !sin(\Omega*(t-t0))Exp(-((t-t0)/tau0)**2)
       ! tau1=tau0/pi2
       ! zp=cmplx(t-t0,tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       ! zm=cmplx(t-t0,-tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       ! ftime =-real(sqrt(pi/2.d0)/2.d0*tau1*exp(-(tau1*w0)**2/2.d0)*(zerf(zm)+zerf(zp)),8)
       tau0 = Ncycles/Omega0
       zp = cmplx((t-t0)/tau0 , tau0*Omega0*pi/2.d0)
       zm = cmplx((t-t0)/tau0 ,-tau0*Omega0*pi/2.d0)
       ftime = -dimag(sqrt(pi)*tau0/4.d0*exp(-0.25d0*(tau0*Omega0*pi)**2)*(zerf(zp)-zerf(zm)))
       Afield=Ek*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("ramp")                !RAMP TO CONSTANT DC-FIELD:
       ftime=-(24.d0*pi*(t+(t-t0)*step(t-t0)+2.d0*(t1-t)*step(t-t0)*step(t-t1)-&
            2.d0*(t0-t1)*step(t-t0)*step(t0-t1))+                              &
            27.d0*t0*(step(t-t0)-1.d0)*Sin(pi*t/t0) - &
            t0*(step(t-t0)-1.d0)*Sin(3.d0*pi*t/t0))/48.d0/pi
       Afield=Ek*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

       !!add more here:
    end select
    !-----------------------------
  end function Afield




  !***************************************************************
  !***************************************************************
  !***************************************************************


end MODULE ELECTRIC_FIELD
