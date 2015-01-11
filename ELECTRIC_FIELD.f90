MODULE ELECTRIC_FIELD
  USE NEQ_CONTOUR
  USE NEQ_VARS_GLOBAL
  USE NEQ_INPUT_VARS
  USE CONSTANTS
  USE DERIVATE
  USE IOTOOLS
  USE FUNCTIONS
  implicit none
  private
  public :: Afield
  public :: set_efield_vector
  public :: print_field

  !ELECTRIC FIELD
  !=========================================================  
  real(8),dimension(3)                   :: Ak,Ek         !Electric field vector potential and vector
  public :: Ak
  public :: Ek

contains

  !+------------------------------------------------------------+
  !PURPOSE: set the normalized electric field versors using given direction
  !+------------------------------------------------------------+
  subroutine set_efield_vector(params)
    type(kb_contour_params) :: params
    real(8)          :: modulo
    integer          :: i
    logical          :: check
    !Normalize the Electric Field components
    !Keep unaltered the Electric Field Strenght Efield=E0
    modulo=sqrt(dot_product(Evect,Evect))
    if(modulo/=0.d0)Evect=Evect/modulo
    Ek = Evect
    write(*,*)"|E|=E0="//trim(txtfy(Efield/modulo))
    check=.false.
    check=field_type=="dc".OR.&
         field_type=="ac".OR.&
         field_type=="acdc".OR.&
         field_type=="pulse".OR.&
         field_type=="ramp"
    if(.not.check)stop "ELECTRIC_FIELD/set_efield_vector: wrong field_type. set:dc,ac,acdc,pulse,ramp"
    call print_field(params%t)
  end subroutine set_efield_vector

  subroutine print_field(t)
    real(8),dimension(:)       :: t
    integer                    :: i
    real(8),dimension(3)       :: A
    real(8),dimension(size(t)) :: Ax,Ay,Az,Ex,Ey,Ez
    do i=1,size(t)
       A=Afield(t(i))
       Ax(i)=A(1)
       Ay(i)=A(2)
       Az(i)=A(3)
    enddo
    Ex = deriv(Ax,dt)
    Ey = deriv(Ay,dt)
    Ez = deriv(Az,dt)
    open(10,file="Avector_shape.ipt")
    open(11,file="Efield_shape.ipt")
    do i=1,size(t)
       write(10,"(4F21.12)")t(i),Afield(t(i))
       write(11,"(4F21.12)")t(i),Ex(i),Ey(i),Ez(i)
    enddo
    close(10)
    close(11)
    if(field_type=="ac")write(*,*)"Root condition: "//trim(txtfy(bessel_j0(Efield/Omega0)))
  end subroutine print_field


  !+------------------------------------------------------------+
  !PURPOSE : 
  !+------------------------------------------------------------+
  function Afield(t)
    real(8),dimension(3) :: E
    real(8),intent(in)   :: t
    real(8)              :: ftime,tau0,tau1
    real(8),dimension(3) :: Afield
    complex(8)           :: zp,zm
    E = Ek
    select case(field_type)
    case ("dc")                !DC ELECTRIC FIELD:
       ftime=-(step(t-t0)*(t-t0 + (t1-t)*step(t-t1) - (t1-t0)*step(t0-t1)))
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("ac")                  !AC ELECTRIC FIELD
       ftime=-sin(Omega0*(t-t0))/Omega0
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("acdc")                !AC+DC ELECTRIC FIELD (super-bloch)
       !ftime=-(t+sin(Omega0*(t-t0))/Omega0)
       ftime =-sin(Omega0*(t-t0))/Omega0
       Afield=E*(Efield*ftime - E1*t)       !A(t) = E0*F(t)*(e_x + e_y)

    case("pulse")               !LIGHT PULSE (for Pump&Probe) 
       !Signal function:
       !cos(\Omega*(t-t0)-pi/2)Exp(-((t-t0)/tau0)**2)
       !sin(\Omega*(t-t0))Exp(-((t-t0)/tau0)**2)
       ! tau1=tau0/pi2
       ! zp=cmplx(t-t0,tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       ! zm=cmplx(t-t0,-tau1**2*w0,8)/(sqrt(2.d0)*tau1)
       ! ftime =-real(sqrt(pi/2.d0)/2.d0*tau1*exp(-(tau1*w0)**2/2.d0)*(zerf(zm)+zerf(zp)),8)
       tau0 = Ncycles/Omega0
       zp = cmplx((t-t0)/tau0 , tau0*Omega0*pi/2.d0,8)
       zm = cmplx((t-t0)/tau0 ,-tau0*Omega0*pi/2.d0,8)
       ftime = -sqrt(pi)*tau0/4.d0*exp(-0.25d0*(tau0*Omega0*pi)**2)*dimag(zerf(zp)-zerf(zm))
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

    case("ramp")                !RAMP TO CONSTANT DC-FIELD:
       ftime=-(24.d0*pi*(t+(t-t0)*step(t-t0)+2.d0*(t1-t)*step(t-t0)*step(t-t1)-&
            2.d0*(t0-t1)*step(t-t0)*step(t0-t1))+                              &
            27.d0*t0*(step(t-t0)-1.d0)*Sin(pi*t/t0) - &
            t0*(step(t-t0)-1.d0)*Sin(3.d0*pi*t/t0))/48.d0/pi
       Afield=E*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

       !!add more here:
    end select
    !-----------------------------
  end function Afield




end MODULE ELECTRIC_FIELD
