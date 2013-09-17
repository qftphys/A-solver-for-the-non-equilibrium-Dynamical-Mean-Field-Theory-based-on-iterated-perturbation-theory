!We use this code to solve 
!integral and integro-differential
!equations on 1-axis.
!This serves as an exercise for the 
!Kadanoff-Baym problem
MODULE RK24_VIDE
  !USE INTEGRATE
  implicit none
  private

  interface vie_rk2
     module procedure d_vie_rk2,c_vie_rk2
  end interface vie_rk2

  interface vide_rk2
     module procedure d_vide_rk2,c_vide_rk2
  end interface vide_rk2

  interface vie_rk4
     module procedure d_vie_rk4,c_vie_rk4
  end interface vie_rk4

  interface vide_rk4
     module procedure d_vide_rk4,c_vide_rk4
  end interface vide_rk4

  public :: vie_rk2
  public :: vie_rk4
  public :: vide_rk2,c_vide_rk2
  public :: vide_rk4

contains

  subroutine get_quadrature_weights(wt,nrk)
    real(8),dimension(:) :: wt
    integer,optional     :: nrk
    integer              :: nrk_
    integer              :: N
    nrk_=4;if(present(nrk))nrk_=nrk
    N=size(wt)
    if(nrk_==4)then
       select case(n)           !n>=3
       case (1)
          wt = 1.d0
       case (2)
          wt = 0.5d0
       case (3)                 !simpson's rule
          wt(1)=1.d0/3.d0
          wt(2)=4.d0/3.d0
          wt(3)=1.d0/3.d0
       case(4)                  !simpson's 3/8 rule
          wt(1)=3.d0/8.d0
          wt(2)=9.d0/8.d0
          wt(3)=9.d0/8.d0
          wt(4)=3.d0/8.d0
       case(5)                  !composite simpson's rule (E,O n)
          wt(1)=1.d0/3.d0
          wt(2)=4.d0/3.d0
          wt(3)=2.d0/3.d0
          wt(4)=4.d0/3.d0
          wt(5)=1.d0/3.d0
       case default             !"gregory's rule" or newton-cote n>=6

          ! !Simpson's 2
          ! if(mod(n-1,2)==0)then
          !    wt(1)=1.d0/3.d0
          !    wt(n)=1.d0/3.d0
          !    wt(2:n-1:2)=4.d0/3.d0
          !    wt(3:n-2:2)=2.d0/3.d0
          ! else
          !    wt(1)=1.d0/3.d0
          !    wt(2:n-4:2)=4.d0/3.d0
          !    wt(3:n-5:2)=2.d0/3.d0
          !    wt(n-3)=17.d0/24.d0
          !    wt(n-2)=9.d0/8.d0
          !    wt(n-1)=9.d0/8.d0
          !    wt(n)=3.d0/8.d0
          ! endif

          ! This has smaller oscillations, though still finite
          ! Gregory's
          wt(1)=3.d0/8.d0
          wt(2)=7.d0/6.d0
          wt(3)=23.d0/24.d0
          wt(4:n-3)=1.d0
          wt(n-2)=23.d0/24.d0
          wt(n-1)=7.d0/6.d0
          wt(n)=3.d0/8.d0

       end select

    elseif(nrk_==2)then
       wt(1) = 0.5d0
       wt(2:n-1)=1.d0
       wt(n) = 0.5d0
    else
       stop "error in +get_quadrature_weights: nrk != 2,4" 
    end if
  end subroutine get_quadrature_weights



  !#######################################################################
  !############################# Runge-Kutta 2 #########################
  !#######################################################################


  subroutine d_vie_rk2(h,f,q,ker)
    real(8), intent(in)                :: h !time-step
    real(8), dimension(:), intent(out) :: f !unknown function
    real(8), dimension(size(f))        :: wt
    integer                  :: n,j,nt
    real(8)                  :: a,b  
    interface
       function q(i)                        !termine noto
         integer, intent(in) :: i !time label
         real(8)             :: q !function value
       end function q
       function Ker(i,j)                      !Kernel
         integer, intent(in) :: i,j !time labels of the 
         real(8)             :: ker !kernel value
       end function Ker
    end interface
    Nt=size(f)
    f(1) = q(1)
    do n=2,Nt
       call get_quadrature_weights(wt(1:n),2)
       wt=wt*h
       b = q(n)
       do j=1,n-1
          b = b + ker(n,j)*f(j)*wt(j)
       end do
       a = 1.d0 - ker(n,n)*wt(n)
       f(n)=b/a
    enddo
  end subroutine d_vie_rk2
  !----------------------
  subroutine c_vie_rk2(h,f,q,ker)
    real(8), intent(in)                   :: h !time-step
    complex(8), dimension(:), intent(out) :: f !unknown function
    real(8), dimension(size(f))           :: wt
    integer                               :: n,j,nt
    complex(8)                            :: a,b  
    interface
       function q(i)                        !termine noto
         integer, intent(in)              :: i !time label
         complex(8)                       :: q !function value
       end function q
       function Ker(i,j)                      !Kernel
         integer, intent(in)              :: i,j !time labels of the 
         complex(8)                       :: ker !kernel value
       end function Ker
    end interface
    Nt=size(f)
    f(1) = q(1)
    do n=2,Nt
       call get_quadrature_weights(wt(1:n),2)
       wt=wt*h
       b = q(n)
       do j=1,n-1
          b = b + ker(n,j)*f(j)*wt(j)
       end do
       a = 1.d0 - ker(n,n)*wt(n)
       f(n)=b/a
    enddo
  end subroutine c_vie_rk2



  subroutine d_vide_rk2(h,f,f0,hk,q,ker)
    real(8), intent(in)                :: h !time-step
    real(8), intent(in)                :: f0 !initial-condition
    real(8), dimension(:), intent(inout) :: f !unknown function
    real(8), dimension(size(f))        :: wt,fexact,df
    integer                            :: n,j,nt,i
    real(8)                            :: a,b,c,dfold,time
    interface
       function q(i)                        !termine noto
         integer, intent(in)           :: i !time label
         real(8)                       :: q !function value
       end function q
       function hk(i)                        !linear coupling
         integer, intent(in)           :: i  !time label
         real(8)                       :: hk !function value
       end function hk
       function Ker(i,j)                      !Kernel
         integer, intent(in)           :: i,j !time labels of the 
         real(8)                       :: ker !kernel value
       end function Ker
    end interface
    Nt=size(f)
    ! !------
    ! !DEBUG:
    ! !Save input array as exact solution: used to bypass the block-by-block method temporarily.
    ! fexact=f
    ! do i=1,Nt
    !    time = dble(i-1)*h
    !    df(i) = -exp(-time)*cosh(time)+exp(-time)*sinh(time)
    !    !
    !    !df(i) = exp(-time)*(time-1.d0)
    ! enddo
    ! !-------
    f=0.d0
    f(1) = f0
    dfold= hk(1)*f(1) + q(1)
    ! time=(1.d0-1.d0)*h
    ! write(20,*)time,abs(f(1)-fexact(1)),abs(dfold-df(1))
    do n=2,Nt
       call get_quadrature_weights(wt(1:n),2);wt=wt*h
       B = f(n-1) + ( dfold + q(n) )*h/2.d0
       C = 0.d0
       do j=1,n-1
          C = C + ker(n,j)*f(j)*wt(j)
       enddo
       A = 1.d0 - hk(n)*h/2.d0 - ker(n,n)*h*h/4.d0
       B = B + C*h/2.d0
       f(n) = B/A
       dfold= hk(n)*f(n) + q(n) + C + ker(n,n)*f(n)*wt(n)
       ! time=(dble(n)-1.d0)*h
       ! write(20,*)time,abs(f(n)-fexact(n)),abs(dfold-df(n))
    enddo
  end subroutine d_vide_rk2

  subroutine c_vide_rk2(h,f,f0,hk,q,ker)
    real(8),intent(in)                     :: h  !time-step
    complex(8),intent(in)                  :: f0 !initial-condition
    complex(8),dimension(:), intent(inout) :: f  !unknown function
    real(8), dimension(size(f))             :: wt
    integer                                 :: n,j,nt
    complex(8)                              :: a,b,c,dfold
    interface
       function q(i)                        !termine noto
         integer, intent(in)                :: i !time label
         complex(8)                         :: q !function value
       end function q
       function hk(i)                        !linear coupling
         integer, intent(in)                :: i  !time label
         complex(8)                         :: hk !function value
       end function hk
       function Ker(i,j)                      !Kernel
         integer, intent(in)                :: i,j !time labels of the 
         complex(8)                         :: ker !kernel value
       end function Ker
    end interface
    Nt=size(f)
    f=cmplx(0.d0,0.d0,8)
    f(1) = f0
    dfold= hk(1)*f(1) + q(1)
    do n=2,Nt
       call get_quadrature_weights(wt(1:n),2)
       wt=wt*h
       B = f(n-1) + ( dfold + q(n) )*h/2.d0
       C = 0.d0
       do j=1,n-1
          C = C + ker(n,j)*f(j)*wt(j)
       enddo
       A = 1.d0 - hk(n)*h/2.d0 - ker(n,n)*h*h/4.d0
       B = B + C*h/2.d0
       f(n) = B/A
       dfold= hk(n)*f(n) + q(n) + C + ker(n,n)*f(n)*wt(n)
    enddo
  end subroutine c_vide_rk2





  !#######################################################################
  !############################# Runge-Kutta 4/5 #########################
  !#######################################################################





  subroutine d_vie_rk4(h,f,q,ker)
    real(8), intent(in)                :: h !time-step
    real(8), dimension(:), intent(out) :: f !unknown function
    real(8), dimension(size(f))        :: wt
    integer                            :: n,j,nt
    real(8)                            :: a,b,ker0,detA
    real(8),dimension(2,2) :: Alin,Ainv
    real(8),dimension(2)   :: blin
    interface
       function q(i)                        !termine noto
         integer, intent(in)           :: i !time label
         real(8)                       :: q !function value
       end function q
       function Ker(i,j)                      !Kernel
         integer, intent(in)           :: i,j !time labels of the 
         real(8)                       :: ker !kernel value
       end function Ker
    end interface
    Nt=size(f)
    !Solve n=1
    f(1) = q(1)
    !Solve n=2,3 together (block-by-block method) interpolate the kernel at the mid-point:
    ker0 = 3.d0/8.d0*ker(2,1) + 3.d0/4.d0*ker(2,2) - 1.d0/8.d0*ker(2,3)
    Alin(1,1) = h/2.d0*ker0 + h/6.d0*ker(2,2) - 1.d0
    Alin(1,2) =-h/12.d0*ker0
    Alin(2,1) = 4.d0*h/3.d0*ker(3,2)
    Alin(2,2) = h/3.d0*ker(3,3)-1.d0
    blin(1) = -q(2)-h/6.d0*ker(2,1)*f(1)-h/4.d0*ker0*f(1)
    blin(2) = -q(3)-h/3.d0*ker(3,1)*f(1)
    !
    detA=Alin(1,1)*Alin(2,2)-Alin(1,2)*Alin(2,1)
    Ainv(1,1)=Alin(2,2)/detA
    Ainv(2,2)=Alin(1,1)/detA
    Ainv(1,2)=-Alin(1,2)/detA
    Ainv(2,1)=-Alin(2,1)/detA
    f(2:3) = matmul(Ainv,blin)
    !Solve N>=4
    do n=4,Nt
       call get_quadrature_weights(wt(1:n),4)
       wt=wt*h
       b = q(n)
       do j=1,n-1
          b = b + ker(n,j)*f(j)*wt(j)
       end do
       a = 1.d0 - ker(n,n)*wt(j)
       f(n)=b/a
    enddo
  end subroutine d_vie_rk4
  !
  subroutine c_vie_rk4(h,f,q,ker)
    real(8), intent(in)                   :: h !time-step
    complex(8), dimension(:), intent(out) :: f !unknown function
    real(8), dimension(size(f))           :: wt
    integer                               :: n,j,nt
    complex(8)                            :: a,b,ker0,detA
    complex(8),dimension(2,2)             :: Alin,Ainv
    complex(8),dimension(2)               :: blin
    interface
       function q(i)                        !termine noto
         integer, intent(in)              :: i !time label
         complex(8)                       :: q !function value
       end function q
       function Ker(i,j)                      !Kernel
         integer, intent(in)              :: i,j !time labels of the 
         complex(8)                       :: ker !kernel value
       end function Ker
    end interface
    Nt=size(f)
    !Solve n=1
    f(1) = q(1)
    !Solve n=2,3 together (block-by-block method) interpolate the kernel at the mid-point:
    ker0 = 3.d0/8.d0*ker(2,1) + 3.d0/4.d0*ker(2,2) - 1.d0/8.d0*ker(2,3)
    Alin(1,1) = h/2.d0*ker0 + h/6.d0*ker(2,2) - 1.d0
    Alin(1,2) =-h/12.d0*ker0
    Alin(2,1) = 4.d0*h/3.d0*ker(3,2)
    Alin(2,2) = h/3.d0*ker(3,3)-1.d0
    blin(1) = -q(2)-h/6.d0*ker(2,1)*f(1)-h/4.d0*ker0*f(1)
    blin(2) = -q(3)-h/3.d0*ker(3,1)*f(1)
    !
    detA=Alin(1,1)*Alin(2,2)-Alin(1,2)*Alin(2,1)
    Ainv(1,1)=Alin(2,2)/detA
    Ainv(2,2)=Alin(1,1)/detA
    Ainv(1,2)=-Alin(1,2)/detA
    Ainv(2,1)=-Alin(2,1)/detA
    f(2:3) = matmul(Ainv,blin)
    !Solve N>=4
    do n=4,Nt
       call get_quadrature_weights(wt(1:n),4)
       wt=wt*h
       b = q(n)
       do j=1,n-1
          b = b + ker(n,j)*f(j)*wt(j)
       end do
       a = 1.d0 - ker(n,n)*wt(j)
       f(n)=b/a
    enddo
  end subroutine c_vie_rk4



  subroutine d_vide_rk4(h,f,f0,ek,q,ker)
    real(8), intent(in)                :: h !time-step
    real(8), intent(in)                :: f0 !initial-condition
    real(8), dimension(:), intent(out) :: f !unknown function
    real(8), dimension(size(f))        :: wt
    integer                            :: n,j,nt,it,i
    real(8)                            :: a,b,c,c1,c2,df1,df2,ker0,ker01,ker02,time,detA
    real(8),dimension(size(f))         :: df,fexact
    real(8),dimension(2,2) :: Alin,Ainv
    real(8),dimension(2)   :: blin
    interface
       function q(i)                        !termine noto
         integer, intent(in)           :: i !time label
         real(8)                       :: q !function value
       end function q
       function ek(i)                        !linear coupling
         integer, intent(in)           :: i  !time label
         real(8)                       :: ek !function value
       end function ek
       function Ker(i,j)                      !Kernel
         integer, intent(in)           :: i,j !time labels of the 
         real(8)                       :: ker !kernel value
       end function Ker
    end interface
    real(8) :: ran(3)
    Nt=size(f)
    !------
    !DEBUG:
    !Save input array as exact solution: used to bypass the block-by-block method temporarily.
    fexact=f
    do i=1,Nt
       time = dble(i-1)*h
       df(i) = -exp(-time)*cosh(time)+exp(-time)*sinh(time)
       !
       !df(i) = exp(-time)*(time-1.d0)
    enddo

    !-------
    f=0.d0
    !Solve n=1
    f(1) = f0
    df1  = ek(1)*f(1) + q(1)
    time=(1.d0-1.d0)*h
    write(40,*)time,abs(f(1)-fexact(1)),abs(df1-df(1))

    !Solve n=2&3 together (block-by-block method) interpolate the kernel at the mid-point:
    ker0 = 3.d0/8.d0*ker(2,1) + 3.d0/4.d0*ker(2,2) - 1.d0/8.d0*ker(2,3)
    ker01= 3.d0*ker0/2.d0 + ker(2,1)
    ker02= 3.d0*ker0      + ker(2,2)
    !
    Alin(1,1) =  2.d0/3.d0*ek(2)*h + h/9.d0*ker02*h - h/9.d0*ker(3,2)*h - 1.d0
    Alin(1,2) = -h/12.d0*ek(3)     - h/18.d0*ker0*h - h/12.d0*h/3.d0*ker(3,3)
    Blin(1)   = -(5.d0/12.d0*df1*h + 2.d0/3.d0*q(2)*h - h/12.d0*q(3) + &
         (1.d0 - h/12.d0*h/3.d0*ker(3,1) + h/12.d0*ker01*h)*f(1) )
    !
    Alin(2,1) =  4.d0/9.d0*h*h*ker(3,2) + 4.d0/3.0*h*ek(2) + 2.d0/9.d0*h*h*ker02
    Alin(2,2) = -h**2/9.d0*ker0 - 1.d0 + ek(3)*h/3.d0 + h/3.d0*h/3.d0*ker(3,3)
    Blin(2)   = -( h/3.d0*q(3) +  4.d0/3.d0*h*q(2) + h/3.d0*df1 + &
         (1.d0 + h*h/9.d0*ker(3,1) + 2.d0*h*h/9.d0*ker01)*f(1) )
    !
    detA=Alin(1,1)*Alin(2,2)-Alin(1,2)*Alin(2,1)
    Ainv(1,1)=Alin(2,2)/detA
    Ainv(2,2)=Alin(1,1)/detA
    Ainv(1,2)=-Alin(1,2)/detA
    Ainv(2,1)=-Alin(2,1)/detA
    f(2:3) = matmul(Ainv,blin)
    !
    ! f(2:3) = fexact(2:3)
    !
    !ONE MUST EVALUATE df1=f'(4-1=3) and df2=f'(4-2=2) to start the loop for n>3
    !df2=f'(n-2)=f'(2)
    df2=ek(2)*f(2) + q(2) + h/2.d0*(ker(2,1)*f(1)+ker(2,2)*f(2))
    time=(2.d0-1.d0)*h
    write(40,*)time,abs(f(2)-fexact(2)),abs(df2-df(2))
    write(*,*)time,abs(f(2)-fexact(2)),abs(df2-df(2))
    !df1=f'(n-1)=f'(3)
    df1=ek(3)*f(3) + q(3) + h/3.d0*(ker(3,1)*f(1) + 4.d0*ker(3,2)*f(2) + ker(3,3)*f(3))
    time=(3.d0-1.d0)*h
    write(40,*)time,abs(f(3)-fexact(3)),abs(df1-df(3))
    write(*,*)time,abs(f(3)-fexact(3)),abs(df1-df(3))



    ! !debug: USE 2ND ORDER TO GET THE FIRST TWO POINTS:
    ! !THIS GIVES A BETTER APPROXIMATION TO THE FIRST 
    ! !2 POINTS IN F: F(2:3)
    ! df1  = ek(1)*f(1) + q(1)
    ! print*,"debug: 2nd order for f(2:3)"
    ! n=2
    ! call get_quadrature_weights(wt(1:2),2);wt=wt*h
    ! A = 1.d0 - ek(2)*h/2.d0 - ker(2,2)*h*h/4.d0
    ! C = ker(2,1)*f(1)*wt(1)
    ! B = f(1) + ( df1 + q(2) )*h/2.d0 + C*h/2.d0
    ! f(2) = B/A
    ! df2= ek(2)*f(2) + q(2) + C + ker(2,2)*f(2)*wt(2)
    ! write(40,*)dble(n-1)*h,abs(f(n)-fexact(n)),abs(df1-df(n))
    ! write(*,*)dble(n-1)*h,abs(f(n)-fexact(n)),abs(df1-df(n))
    ! !
    ! n=3
    ! call get_quadrature_weights(wt(1:n),2);wt=wt*h
    ! A = 1.d0 - ek(n)*h/2.d0 - ker(n,n)*h*h/4.d0
    ! C = ker(n,1)*f(1)*wt(1)+ ker(n,2)*f(2)*wt(2)
    ! B = f(2) + ( df2 + q(n) )*h/2.d0 + C*h/2.d0
    ! f(n) = B/A
    ! df1= ek(n)*f(n) + q(n) + C + ker(n,n)*f(n)*wt(n)
    ! write(40,*)dble(n-1)*h,abs(f(n)-fexact(n)),abs(df1-df(n))
    ! write(*,*)dble(n-1)*h,abs(f(n)-fexact(n)),abs(df1-df(n))
    ! !

    do n=4,Nt
       call get_quadrature_weights(wt(1:n),4);wt=wt*h
       A = 1.d0 - h/3.d0*ek(n) - ker(n,n)*h/3.d0*wt(n)
       C = 0.d0
       do j=1,n-1
          C = C + ker(n,j)*f(j)*wt(j)
       enddo
       B = f(n-2) + ( df2 + 4.d0*df1 + q(n) )*h/3.d0 + C*h/3.d0
       f(n) = B/A
       df2  = df1
       df1  = ek(n)*f(n) + q(n) + C + ker(n,n)*f(n)*wt(n)

       time=dble(n-1)*h
       write(40,*)time,abs(f(n)-fexact(n)),abs(df1-df(n))
       write(*,*)time,abs(f(n)-fexact(n)),abs(df1-df(n))
    enddo
  end subroutine d_vide_rk4



  subroutine c_vide_rk4(h,f,f0,ek,q,ker)
    real(8), intent(in)                   :: h !time-step
    complex(8), intent(in)                :: f0 !initial-condition
    complex(8), dimension(:), intent(out) :: f !unknown function
    real(8), dimension(size(f))           :: wt
    integer                               :: n,j,nt,it,i
    complex(8)                            :: a,b,c,c1,c2,df1,df2,ker0,ker01,ker02,time,detA
    complex(8),dimension(size(f))         :: df,fexact
    complex(8),dimension(2,2)             :: Alin,Ainv
    complex(8),dimension(2)               :: blin
    interface
       function q(i)                        !termine noto
         integer, intent(in)              :: i !time label
         complex(8)                       :: q !function value
       end function q
       function ek(i)                        !linear coupling
         integer, intent(in)              :: i  !time label
         complex(8)                       :: ek !function value
       end function ek
       function Ker(i,j)                      !Kernel
         integer, intent(in)              :: i,j !time labels of the 
         complex(8)                       :: ker !kernel value
       end function Ker
    end interface
    Nt=size(f)
    ! !------
    ! !DEBUG:
    ! !Save input array as exact solution: used to bypass the block-by-block method temporarily.
    ! fexact=f
    ! do i=1,Nt
    !    time = dble(i-1)*h
    !    df(i) = -exp(-time)*cosh(time)+exp(-time)*sinh(time)
    !    !
    !    !df(i) = exp(-time)*(time-1.d0)
    ! enddo
    ! !-------
    f=0.d0
    !Solve n=1
    f(1) = f0
    df1  = ek(1)*f(1) + q(1)
    !Solve n=2&3 together (block-by-block method) interpolate the kernel at the mid-point:
    ker0 = 3.d0/8.d0*ker(2,1) + 3.d0/4.d0*ker(2,2) - 1.d0/8.d0*ker(2,3)
    ker01= 3.d0*ker0/2.d0 + ker(2,1)
    ker02= 3.d0*ker0      + ker(2,2)
    !
    Alin(1,1) =  2.d0/3.d0*ek(2)*h + h/9.d0*ker02*h - h/9.d0*ker(3,2)*h - 1.d0
    Alin(1,2) = -h/12.d0*ek(3)     - h/18.d0*ker0*h - h/12.d0*h/3.d0*ker(3,3)
    Blin(1)   = -(5.d0/12.d0*df1*h + 2.d0/3.d0*q(2)*h - h/12.d0*q(3) + &
         (1.d0 - h/12.d0*h/3.d0*ker(3,1) + h/12.d0*ker01*h)*f(1) )
    !
    Alin(2,1) =  4.d0/9.d0*h*h*ker(3,2) + 4.d0/3.0*h*ek(2) + 2.d0/9.d0*h*h*ker02
    Alin(2,2) = -h**2/9.d0*ker0 - 1.d0 + ek(3)*h/3.d0 + h/3.d0*h/3.d0*ker(3,3)
    Blin(2)   = -( h/3.d0*q(3) +  4.d0/3.d0*h*q(2) + h/3.d0*df1 + &
         (1.d0 + h*h/9.d0*ker(3,1) + 2.d0*h*h/9.d0*ker01)*f(1) )
    !
    detA=Alin(1,1)*Alin(2,2)-Alin(1,2)*Alin(2,1)
    Ainv(1,1)=Alin(2,2)/detA
    Ainv(2,2)=Alin(1,1)/detA
    Ainv(1,2)=-Alin(1,2)/detA
    Ainv(2,1)=-Alin(2,1)/detA
    f(2:3) = matmul(Ainv,blin)
    !ONE MUST EVALUATE df1=f'(4-1=3) and df2=f'(4-2=2) to start the loop for n>3
    !df2=f'(n-2)=f'(2)
    df2=ek(2)*f(2) + q(2) + h/2.d0*(ker(2,1)*f(1)+ker(2,2)*f(2))
    !df1=f'(n-1)=f'(3)
    df1=ek(3)*f(3) + q(3) + h/3.d0*(ker(3,1)*f(1) + 4.d0*ker(3,2)*f(2) + ker(3,3)*f(3))
    do n=4,Nt
       call get_quadrature_weights(wt(1:n),4)
       wt=wt*h
       A = 1.d0 - h/3.d0*ek(n) - ker(n,n)*h/3.d0*wt(n)
       B = f(n-2) + ( df2 + 4.d0*df1 + q(n) )*h/3.d0
       C = 0.d0
       do j=1,n-1
          C = C + ker(n,j)*f(j)*wt(j)
       enddo
       B = B + C*h/3.d0
       f(n) = B/A       
       df2  = df1
       df1  = ek(n)*f(n) + q(n) + C + ker(n,n)*f(n)*wt(n)
    enddo
  end subroutine c_vide_rk4




END MODULE RK24_VIDE

