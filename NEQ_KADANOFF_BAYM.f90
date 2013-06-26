MODULE NEQ_KADANOFF_BAYM
  USE NEQ_VARS_GLOBAL
  USE ELECTRIC_FIELD
  implicit none
  private


  real(8),dimension(4,4)   :: A
  real(8),dimension(4)     :: B,C

  !k-dependent GF:
  type(keldysh_contour_gf) :: Gk

  type(keldysh_contour_gf) :: tmpG1,tmpG2 

  type(keldysh_contour_gf) :: Icoll                          !+-> collision integral up to time tn 
  type(keldysh_contour_gf) :: delta_Icoll                    !+-> correction to the collision integral for times in [tn + t_rk]

  type(keldysh_contour_gf) :: Sfit

  public                   :: neq_get_localgf


contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the k-sum to get the local Green's functions (R,>,<)
  !+-------------------------------------------------------------------+
  subroutine neq_get_localgf()
    call kadanoff_baym_to_localgf()
    call print_out_Gloc()
  end subroutine neq_get_localgf


  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the local Green's functions G^<,>
  !along the real time axis, using a discretized verion of the Kadanoff-Baym 
  !equations.
  !+-------------------------------------------------------------------+
  subroutine kadanoff_baym_to_localgf()
    integer                  :: istep,i,j,ik,k

    call msg("Entering Kadanoff-Baym")

    call allocate_keldysh_contour_gf(Gk,Nstep*(Nstep+1)/2)
    call allocate_keldysh_contour_gf(Icoll,Nstep)
    call allocate_keldysh_contour_gf(tmpG1,Nstep)
    call allocate_keldysh_contour_gf(tmpG2,Nstep)
    call allocate_keldysh_contour_gf(delta_Icoll,Nstep)

    !Set to Zero GF and nk:
    locG   = zero
    nk     = 0.d0

    call RKcoeff(A,B,C)

    !=============START K-POINTS LOOP======================
    call start_timer
    do ik=1,Lk
       Gk=zero
       Icoll=zero
       tmpG1=zero
       tmpG2=zero
       delta_Icoll=zero
       !t-step loop
       do istep=1,nstep-1
          call step_dynamics(ik,istep)
       enddo

       !sum over k-point
       locG%less = locG%less + Gk%less*wt(ik)
       locG%gtr  = locG%gtr  + Gk%gtr*wt(ik)
       do istep=1,nstep
          k = pack_index_tri(istep,istep)
          nk(istep,ik)=-xi*Gk%less(k)
       enddo
       call eta(ik,Lk)
    enddo
    call stop_timer
    !=============END K-POINTS LOOP======================
    call deallocate_keldysh_contour_gf(tmpG1)
    call deallocate_keldysh_contour_gf(tmpG2)
    call deallocate_keldysh_contour_gf(delta_Icoll)
    call deallocate_keldysh_contour_gf(Icoll)
  end subroutine kadanoff_baym_to_localgf



  !+------------------------------------------------------------------+!
  SUBROUTINE step_dynamics(ik,itime)
    integer,intent(in)                  :: ik,itime
    integer                             :: ity,it,it_old
    integer                             :: irk,it_gf
    real(8)                             :: trk,tj,time
    !/********************\!    
    !+- 4-STEP PROCEDURE -+!
    !\********************/!    
    time = t(itime)
    !istep==0: provide initial conditions to the real-time part of the contour:
    if(itime==1)then
       it_gf = pack_index_tri(1,1)
       Gk%less(it_gf)= xi*fermi(epsik(ik),beta)
       Gk%gtr(it_gf) = xi*(fermi(epsik(ik),beta)-1.d0)
    endif

    call init_step(itime)

    do irk=1,4
       !vertical sweep t' at a fixed t=itime
       do ity=1,itime
          !+- INTERMEDIATE RK_VALUES -+!
          call get_intermediate_gf(ity,irk,itime,ik) !--> tmpG2 on intermediate/final point, tmpG1 from previous irk-step
          call get_Icoll(ity,irk,itime,ik)           !compute collision integrals lag-term (up to times t_n). It re-used in the next irk-step
          call get_delta_Icoll(ity,irk,itime,ik)     !correction to the collision integrals (= 0 if irk = 1).
          !+- UPDATE LESSER AND GREATER GF -+!       !at t_rk = t(itime) + C(irk)*dt
          !this routine *ADD* the irk-th contribution to the step-sum:
          !y_{n+1} = y_n + \sum_{\irk=1:4}B_\irk F(tmpG2_\irk) where tmpG2 is the actual tmp GF.
          call update_gf(ity,irk,itime,ik)
       end do

       !diagonal terms:
       call get_intermediate_diag_gf(irk,itime,ik)   !as before but just the diagonal propagation
       !this routine evaluate the lag-term of the diagonal collision integral
       !as this performs a loop 1:itime, it might be better to include this 
       !evaluation *inside* the main propagation loop above, updating the 
       !integral step-by-step in the loop.
       !compare with hey_green in slab_KB code.
       call get_diag_Icoll(irk,itime,ik)              !compute diagonal collision integral lag-term
       call get_delta_diag_Icoll(irk,itime,ik)  !correction to the diagonal collision integrals ( = 0 if irk = 1)
       !+- UPDATE DIAGONAL LESSER AND GREATER GREEN FUNCTION -+!
       call update_diag_gf(irk,itime,ik)
       !
       !COPY BACK TMPG1 <-- TMPG2
       do ity=1,itime+1
          tmpG1%less(ity) = tmpG2%less(ity)
          tmpG1%gtr(ity)  = tmpG2%gtr(ity)
       end do
    end do
  end subroutine step_dynamics




  !+- allocate and initialize alla quantities needed in the step routine -+!
  subroutine init_step(itime)
    integer :: ity,it_gf,itime
    do ity = 1,itime
       it_gf = pack_index_tri(itime,ity)
       tmpG1%less(ity) = Gk%less(it_gf)
       tmpG1%gtr(ity)  = Gk%gtr(it_gf)
       tmpG2%less(ity) = Gk%less(it_gf)
       tmpG2%gtr(ity)  = Gk%gtr(it_gf)
    end do
    ity = itime+1
    it_gf = pack_index_tri(itime,itime)
    tmpG1%less(ity) = Gk%less(it_gf)
    tmpG1%gtr(ity)  = Gk%gtr(it_gf)
    tmpG2%less(ity) = Gk%less(it_gf)
    tmpG2%gtr(ity)  = Gk%gtr(it_gf)
    delta_Icoll%less(:) = zero
    delta_Icoll%gtr(:)  = zero
  end subroutine init_step


  !+-------------------------------------------+!
  !+- COMPUTE THE INTERMEDIATE GREEN FUNCTION -+!
  !+-------------------------------------------+!
  SUBROUTINE get_intermediate_gf(ity,irk,itime,ik)
    integer,intent(in) :: ity,irk,itime,ik
    integer            :: it_gf,j
    real(8)            :: tj,ti,time
    it_gf = pack_index_tri(itime,ity)
    time  = t(itime)
    if(irk.eq.1) then
       tmpG1%less(ity) = Gk%less(it_gf)
       tmpG1%gtr(ity)  = Gk%gtr(it_gf)
       !
       tmpG2%less(ity) = tmpG1%less(ity)
       tmpG2%gtr(ity)  = tmpG1%gtr(ity)
    else
       j  = irk - 1              !this is because A_ij in RK coefficients reduces to \delta_{i}{i+1}
       ti = t(ity)
       tj = time + C(j)*dt
       tmpG2%less(ity) = Gk%less(it_gf) + A(irk,j)*dt*stepping_gf( & !add ik & itime-indices for the Hamiltonian
            ik,tj,tmpG1%less(ity),Icoll%less(ity),delta_Icoll%less(ity))
       !
       tmpG2%gtr(ity)  = Gk%gtr(it_gf)  + A(irk,j)*dt*stepping_gf( &
            ik,tj,tmpG1%gtr(ity),Icoll%gtr(ity),delta_Icoll%gtr(ity))
    end if
  end subroutine get_intermediate_gf




  !+--------------------------------------+!
  !+- COMPUTE TIME CORRELATION INTEGRALS -+!
  !+--------------------------------------+!
  subroutine get_Icoll(ity,irk,itime,ik)
    integer,intent(in) :: ity,irk,itime,ik
    integer            :: it,it_gf
    real(8)            :: wtrapz
    real(8)            :: ti,trk,time
    logical            :: swap_tt
    time= t(itime)
    trk = time + C(irk)*dt

    if(irk==3)return

    if(irk==1.AND.itime>2.AND.ity/=itime)then
       swap_tt = .false.
       wtrapz  = 0.5d0*dt
       it      = itime-1
       ti      = t(it)
       it_gf   = pack_index_tri(it,ity)
       Icoll%less(ity) = Icoll%less(ity) + timeCorr_ft_less(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz
       Icoll%gtr(ity)  = Icoll%gtr(ity)  + timeCorr_ft_gtr(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz
       it      = itime
       ti      = t(it)
       it_gf   = pack_index_tri(it,ity)
       Icoll%less(ity) = Icoll%less(ity) + timeCorr_ft_less(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz
       Icoll%gtr(ity)  = Icoll%gtr(ity)  + timeCorr_ft_gtr(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz

    else

       Icoll%less(ity)=zero
       Icoll%gtr(ity)=zero
       !+- integrate along the y-direction -+!
       if(ity.gt.1) then
          swap_tt = .true.
          do it = 1,ity  
             if(it.eq.1.or.it.eq.ity) then
                wtrapz = 0.5d0*dt
             else
                wtrapz = dt
             end if
             ti    = t(it)
             it_gf = pack_index_tri(ity,it)
             !
             Icoll%less(ity) = Icoll%less(ity) + &
                  timeCorr_ft_less(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz + &
                  timeCorr_ft_prime_less(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz
             !
             Icoll%gtr(ity)  = Icoll%gtr(ity)  + &
                  timeCorr_ft_gtr(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz + &
                  timeCorr_ft_prime_gtr(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz
          end do
       end if
       !
       !+- integrate along the x-direction -+!
       if(ity.lt.itime) then
          swap_tt = .false.
          do it = ity,itime 
             if(it.eq.ity.or.it.eq.itime) then
                wtrapz = 0.5d0*dt
             else
                wtrapz = dt
             end if
             ti    = t(it)
             it_gf = pack_index_tri(it,ity)
             !
             Icoll%less(ity) = Icoll%less(ity) + &
                  timeCorr_ft_less(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz
             !
             Icoll%gtr(ity)  = Icoll%gtr(ity)  + &
                  timeCorr_ft_gtr(trk,ti,Gk%less(it_gf),Gk%gtr(it_gf),swap_tt)*wtrapz
          end do
       end if

    end if

  end subroutine get_Icoll



  subroutine get_delta_icoll(ity,irk,itime,ik)
    integer,intent(in) :: ity,irk,itime,ik
    integer            :: j
    logical            :: swap_tt
    real(8)            :: time,tj,trk
    if(irk==1)return
    j   = irk - 1
    time= t(itime)
    tj  = time + C(j)*dt
    trk = time + C(irk)*dt
    swap_tt = .false.
    delta_Icoll%less(ity) = A(irk,j)*timeCorr_ft_less(trk,tj,tmpG1%less(ity),tmpG1%gtr(ity),swap_tt)*dt
    delta_Icoll%gtr(ity)  = A(irk,j)*timeCorr_ft_gtr(trk,tj,tmpG1%less(ity),tmpG1%gtr(ity),swap_tt)*dt
  end subroutine get_delta_icoll




  subroutine update_gf(ity,irk,itime,ik)
    integer,intent(in) :: ity,irk,itime,ik
    integer            :: it_gf,it_old
    real(8) :: trk,time
    it_gf  = pack_index_tri(itime+1,ity) 
    it_old = pack_index_tri(itime,ity)
    time= t(itime)
    trk = time + C(irk)*dt
    if(irk.eq.1) then
       Gk%less(it_gf) = Gk%less(it_old)
       Gk%gtr(it_gf)  = Gk%gtr(it_old)
    end if
    Gk%less(it_gf) = Gk%less(it_gf) + B(irk)*dt*    & 
         stepping_gf(ik,trk,tmpG2%less(ity),Icoll%less(ity),delta_Icoll%less(ity))
    Gk%gtr(it_gf)  = Gk%gtr(it_gf) + B(irk)*dt*      & 
         stepping_gf(ik,trk,tmpG2%gtr(ity),Icoll%gtr(ity),delta_Icoll%gtr(ity)) 
  end subroutine update_gf






  !##################################################################
  ! DIAG ROUTINES
  !##################################################################
  subroutine get_intermediate_diag_gf(irk,itime,ik)
    integer,intent(in) :: irk,itime,ik
    integer            :: it_old,j
    real(8) :: tj,time
    it_old = pack_index_tri(itime,itime) !this the OLD tip of the square. we are propagating to itime+1,itime+1

    if(irk.eq.1) then
       tmpG1%less(itime+1) = Gk%less(it_old)
       tmpG1%gtr(itime+1)  = Gk%gtr(it_old)
       tmpG2%less(itime+1) = tmpG1%less(itime+1)
       tmpG2%gtr(itime+1)  = tmpG1%gtr(itime+1)
    else
       j = irk - 1
       time=t(itime)
       tj = time + C(j)*dt
       tmpG2%less(itime+1) = Gk%less(it_old) + A(irk,j)*dt*&  
            diag_stepping_gf(ik,tj,tmpG1%less(itime+1),Icoll%less(itime+1),delta_Icoll%less(itime+1))
       !
       tmpG2%gtr(itime+1) = Gk%gtr(it_old) + A(irk,j)*dt*&
            diag_stepping_gf(ik,tj,tmpG1%gtr(itime+1),Icoll%gtr(itime+1),delta_Icoll%gtr(itime+1))
    end if
  end subroutine get_intermediate_diag_gf


  !+--------------------------------------+!
  !+- COMPUTE TIME CORRELATION INTEGRALS -+!
  !+--------------------------------------+!
  subroutine get_diag_Icoll(irk,itime,ik)
    integer,intent(in) :: irk,itime,ik
    integer            :: it,it_gf
    real(8)            :: wtrapz
    real(8)            :: ti,trk,time
    logical            :: swap_tt
    time= t(itime)
    trk = time + C(irk)*dt
    swap_tt=.true.
    Icoll%less(itime+1)=zero
    Icoll%gtr(itime+1) =zero
    !+- COLLISION INTEGRALS ALONG THE PROPAGATION FRONT -+!
    !
    if(itime==1)return
    do it=1,itime
       wtrapz=dt
       if(it==1.OR.it==itime)wtrapz = 0.5d0*dt
       ti = t(it)
       Icoll%less(itime+1) = Icoll%less(itime+1) + &
            timeCorr_ft_less(trk,ti,tmpG2%less(it),tmpG2%gtr(it),swap_tt)*wtrapz + &
            timeCorr_ft_prime_less(trk,ti,tmpG2%less(it),tmpG2%gtr(it),swap_tt)*wtrapz
       Icoll%gtr(itime+1) = Icoll%gtr(itime+1)+ &
            timeCorr_ft_gtr(trk,ti,tmpG2%less(it),tmpG2%gtr(it),swap_tt)*wtrapz + &
            timeCorr_ft_prime_gtr(trk,ti,tmpG2%less(it),tmpG2%gtr(it),swap_tt)*wtrapz
    enddo
  end subroutine get_diag_icoll


  subroutine get_delta_diag_icoll(irk,itime,ik)
    integer,intent(in) :: irk,itime,ik
    integer            :: j,ity
    logical            :: swap_tt
    real(8)            :: tj,time,trk
    if(irk.eq.1)return
    swap_tt = .true.
    j    = irk - 1
    time = t(itime)
    tj   = time + C(j)*dt
    trk  = time + C(irk)*dt
    ity  = itime+ 1
    delta_Icoll%less(itime+1) = A(irk,j)*( timeCorr_ft_less(trk,tj,tmpG1%less(itime+1),tmpG1%gtr(itime+1),swap_tt) + &
         timeCorr_ft_prime_less(trk,tj,tmpG1%less(itime+1),tmpG1%gtr(itime+1),swap_tt) )*dt         
    delta_Icoll%gtr(itime+1)  = A(irk,j)*( timeCorr_ft_gtr(trk,tj,tmpG1%less(itime+1),tmpG1%gtr(itime+1),swap_tt) + &
         timeCorr_ft_prime_gtr(trk,tj,tmpG1%less(itime+1),tmpG1%gtr(itime+1),swap_tt) )*dt
  end subroutine get_delta_diag_icoll


  subroutine update_diag_gf(irk,itime,ik)
    integer,intent(in) :: irk,itime,ik
    integer            :: it_gf,it_old
    real(8) :: time,trk
    it_gf  = pack_index_tri(itime+1,itime+1)
    it_old = pack_index_tri(itime,itime)
    time=t(itime)
    trk=time+C(irk)*dt
    if(irk.eq.1) then
       Gk%less(it_gf) = Gk%less(it_old)
       Gk%gtr(it_gf)  = Gk%gtr(it_old)
    end if
    Gk%less(it_gf) = Gk%less(it_gf) + B(irk)*dt*&
         diag_stepping_gf(ik,trk,tmpG2%less(itime+1),Icoll%less(itime+1),delta_Icoll%less(itime+1))
    Gk%gtr(it_gf) = Gk%gtr(it_gf) + B(irk)*dt*& 
         diag_stepping_gf(ik,trk,tmpG2%gtr(itime+1),Icoll%gtr(itime+1),delta_Icoll%gtr(itime+1)) 
  end subroutine update_diag_gf




  subroutine RKcoeff(A,B,C)
    real(8),dimension(4,4),intent(inout) :: A
    real(8),dimension(4),intent(inout)     :: B,C
    A=0.d0
    B=0.d0
    C=0.d0
    A(2,1) = 0.5d0
    A(3,2) = 0.5d0
    A(4,3) = 1.d0
    C(2) = 0.5d0
    C(3) = 0.5d0
    C(4) = 1.d0
    B(1) = 1.d0/6.d0
    B(2) = 1.d0/3.d0
    B(3) = 1.d0/3.d0
    B(4) = 1.d0/6.d0
  end subroutine rkcoeff





  ! ##################################################################
  ! ##################################################################
  ! ##################################################################
  ! ##################################################################
  ! ##################################################################





  !+-------------------------------+!
  !+- TIME CORRELATION INTEGRANDS -+!
  !+-------------------------------+!
  FUNCTION timeCorr_ft_less(t,tau,green_less,green_gtr,swap_tt)  result(ft_less)
    real(8)    :: t,tau
    complex(8) :: green_less
    complex(8) :: green_gtr
    complex(8) :: ft_less
    logical    :: swap_tt
    integer    :: it,jt,it_sigma
    complex(8) :: sigma_less,sigma_gtr,sigma_ret

    it = int(t/dt*2+1.d-5) + 1
    jt = int(tau/dt*2+1.d-5) + 1

    it_sigma = pack_index(it,jt,Nfit)

    Sigma_less = S0%less(it_sigma) !+Sigma_dmft_less(it_sigma)
    Sigma_gtr  = S0%gtr(it_sigma)

    Sigma_ret = Sigma_gtr - Sigma_less
    ft_less   = Sigma_ret*green_less
    if ( swap_tt ) ft_less = Sigma_ret*swap_time_gf(green_less)

  end function timecorr_ft_less



  FUNCTION timeCorr_ft_gtr(t,tau,green_less,green_gtr,swap_tt)  result(ft_gtr)
    real(8)    :: t,tau
    complex(8) :: green_less
    complex(8) :: green_gtr
    complex(8) :: ft_gtr
    logical    :: swap_tt
    complex(8) :: sigma_less,sigma_gtr,sigma_ret
    integer    :: it,jt,it_sigma

    it = int(t/dt*2+1.d-5) + 1
    jt = int(tau/dt*2+1.d-5) + 1

    it_sigma = pack_index(it,jt,Nfit)

    Sigma_less = S0%less(it_sigma) !+Sigma_dmft_less(it_sigma)
    Sigma_gtr  = S0%gtr(it_sigma)
    Sigma_ret = Sigma_gtr - Sigma_less
    ft_gtr   = Sigma_ret*green_gtr
    if ( swap_tt ) ft_gtr = Sigma_ret*swap_time_gf(green_gtr)
  end function timeCorr_ft_gtr







  function timeCorr_ft_prime_less(t,tau,green_less,green_gtr,swap_tt) result(ft_prime_less)
    real(8)                   :: t,tau
    complex(8)                :: green_less
    complex(8)                :: green_gtr
    complex(8)                :: ft_prime_less
    logical                   :: swap_tt
    complex(8)                :: Sigma_less,Sigma_gtr
    integer                   :: it,jt,it_sigma
    it = int(t/dt*2+1.d-5) + 1
    jt = int(tau/dt*2+1.d-5) + 1

    it_sigma = pack_index(it,jt,Nfit)

    Sigma_less    = S0%less(it_sigma) !+Sigma_dmft_less(it_sigma)
    Sigma_gtr     = S0%gtr(it_sigma)
    ft_prime_less = -Sigma_less*(green_gtr - green_less)
    if(swap_tt)ft_prime_less = -Sigma_less*(swap_time_gf(green_gtr)-swap_time_gf(green_less))
  end function timeCorr_ft_prime_less


  FUNCTION timeCorr_ft_prime_gtr(t,tau,green_less,green_gtr,swap_tt) result(ft_prime_gtr)
    real(8)                   :: t,tau
    complex(8)                :: green_less
    complex(8)                :: green_gtr
    complex(8)                :: ft_prime_gtr
    logical                   :: swap_tt
    integer                   :: it,jt,it_sigma
    complex(8)                :: Sigma_less,Sigma_gtr
    it = int(t/dt*2+1.d-5) + 1
    jt = int(tau/dt*2+1.d-5) + 1

    it_sigma = pack_index(it,jt,Nfit)

    Sigma_less    = S0%less(it_sigma) !+Sigma_dmft_less(it_sigma)
    Sigma_gtr     = S0%gtr(it_sigma)

    ft_prime_gtr = -Sigma_gtr*(green_gtr - green_less)
    if(swap_tt)ft_prime_gtr = -Sigma_gtr*(swap_time_gf(green_gtr)-swap_time_gf(green_less))

  end function timeCorr_ft_prime_gtr




  function stepping_gf(ik,time,green,Icoll,delta_Icoll) result(step_gf)
    integer    :: ik,itime
    complex(8) :: green
    complex(8) :: Icoll,delta_Icoll
    complex(8) :: step_gf
    integer :: i,j
    real(8) :: time,ekt
    type(vect2D) :: kt,Ak
    i       = ik2ix(ik)
    j       = ik2iy(ik)
    Ak      = Afield(time,Ek)
    kt      = kgrid(i,j) - Ak
    ekt     = square_lattice_dispersion(kt)
    step_gf = -xi*(Icoll + delta_Icoll) - xi*ekt*green
  end function stepping_gf



  function diag_stepping_gf(ik,time,green,Icoll,delta_Icoll) result(step_gf)
    integer    :: ik,itime
    complex(8) :: green
    complex(8) :: Icoll,delta_Icoll
    complex(8) :: step_gf
    real(8) :: time
    !here you don't have anything right now 'cause in the one-band simple model 
    !kinetic part cancels out. In more general models it DOES NOT!!
    step_gf = -xi*2.d0*dreal(Icoll+delta_Icoll)
  end function diag_stepping_gf


  function swap_time_gf(gf) result(swap_gf)
    complex(8),intent(in) :: gf
    complex(8)            :: swap_gf
    swap_gf = -conjg(gf)
  end function swap_time_gf



  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : print out useful information
  !+-------------------------------------------------------------------+
  subroutine print_out_Gloc()
    integer                       :: i,j,ix,iy,ik
    type(vect2D)                  :: Jk,Ak
    type(vect2D),dimension(nstep) :: Jloc                   !local Current 
    real(8),dimension(nstep)      :: nt,modJloc             !occupation(time)
    ! call write_keldysh_contour_gf(locG,"locG")
    call store_data("nk.neq",nk)
    do i=1,nstep
       nt(i)=-xi*pack_less_tri(i,i,locG)
    enddo
    Jloc=Vzero    
    do ik=1,Lk
       ix=ik2ix(ik);iy=ik2iy(ik)
       do i=1,nstep
          Ak= Afield(t(i),Ek)
          Jk= nk(i,ik)*square_lattice_velocity(kgrid(ix,iy) - Ak)
          Jloc(i) = Jloc(i) +  wt(ik)*Jk
       enddo
    enddo
    call splot("nVStime.neq",t,2.d0*nt,append=.true.)
    if(Efield/=0.d0)call splot("JlocVStime.neq",t,Jloc(:)%x,Jloc(:)%y,append=.true.)
    ! if(fchi)then
    !    call store_data("Chi_11.neq",chi(1,1,:,:))
    !    call store_data("Chi_12.neq",chi(1,2,:,:))
    !    call store_data("Chi_21.neq",chi(2,1,:,:))
    !    call store_data("Chi_22.neq",chi(2,2,:,:))
    ! endif
  end subroutine print_out_Gloc



  !******************************************************************
  !******************************************************************
  !******************************************************************



END MODULE NEQ_KADANOFF_BAYM
