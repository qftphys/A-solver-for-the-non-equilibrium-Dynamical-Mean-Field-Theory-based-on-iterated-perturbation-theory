!########################################################################
!PROGRAM  : KADANOFFBAYM
!TYPE     : module
!PURPOSE  : Evolve the Green's functions G^>,< on the real time axis
!given the initial condition at t=t'=0.
!AUTHORS  : Adriano Amaricci
!########################################################################
module KADANOFBAYM
  !LOCAL:
  USE VARS_GLOBAL
  USE FUNX_NEQ
  implicit none
  private
  public    :: get_Gloc_KadanoffBaym 
  save


contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : KADANOF_BAYM2GLOC
  !TYPE     : Subroutine
  !PURPOSE  : Build the local Green's functions G^<,>
  !along the real time axis, using a discretized verion of the Kadanoff-Baym 
  !equations. The evolution is performed by a two-pass procedure, following:
  ! H.S. K\"ohler, N.H.Kwong, Hashim A. Yousif
  !"A Fortran code for solving the Kadanoff-Baym equations for a homogeneous 
  !fermion system" Computer Physics Communications,123(1999),123-142
  !+-------------------------------------------------------------------+
  subroutine get_Gloc_KadanoffBaym(loop)
    integer   :: istep,i,j,ik,loop
    integer   :: STARTik,ENDik
    call dump("",lines=1)
    call dump("Entering Kadanoff-Baym")
    call dump("----------------------")

    !First loop ever, build the Initial Condition
    if(loop==1)call build_ic

    !Set to Zero loc GF:
    locGless=zero; locGgtr=zero
#ifdef _mix
    locGlceil=zero;locGrceil=zero
#endif

    !Set the evolution operators
    call buildUV

    !Start the time-stepping procedure: solve KB equations
    call dump("Solving KB equations:")

    !=============START K-POINTS LOOP======================
    do ik=1,Lk
       call print_bar(ik,Lk)
       Gkless=zero;   Gkgtr=zero
#ifdef _mix
       Gklceil=zero;  Gkrceil=zero
#endif
       !Recover Initial Condition for Gk^{<,>}
       call recover_ic(ik)
       !======T-STEP LOOP=====================
       do istep=0,nstep-1
          call GFstep(1,ik,istep) !1st-pass
          call GFstep(2,ik,istep) !2nd-pass
       enddo
       !call save_ic(ik)

       locGless(0:nstep,0:nstep) =locGless(0:nstep,0:nstep) + Gkless(0:nstep,0:nstep)*wt(ik)
       locGgtr(0:nstep,0:nstep)  =locGgtr(0:nstep,0:nstep)  + Gkgtr(0:nstep,0:nstep)*wt(ik)
#ifdef _mix
       locGlceil=locGlceil+ Gklceil(:,:)*wt(ik)
       locGrceil=locGrceil+ Gkrceil(:,:)*wt(ik)
#endif

       forall(istep=0:nstep)nk(istep,ik)=-xi*Gkless(istep,istep)
    enddo
    call close_bar
    !=============END K-POINTS LOOP======================
    call dump("")    

    !Gloc^>(t',t)= - Gloc^>(t,t')^T
    forall(i=0:nstep,j=0:nstep,i>j)locGless(i,j)=-conjg(locGless(j,i))
    forall(i=0:nstep,j=0:nstep,i<j)locGgtr(i,j) =-conjg(locGgtr(j,i))

    call print_out_Gloc()
    call dump("",lines=2)
    return
  end subroutine get_Gloc_KadanoffBaym
  !******************************************************************
  !******************************************************************
  !******************************************************************








  !+-------------------------------------------------------------------+
  !PROGRAM  : GFSTEP
  !TYPE     : Subroutine
  !PURPOSE  : Execute the 2pass procedure to solve KB equations:
  !+-------------------------------------------------------------------+
  subroutine GFstep(ips,ik,istep)
    integer :: ips,istep,ik
    integer :: i,j,itau
    integer :: it
    if(ips == 1)then
       !First Pass: get collision integrals up to t=T=istep
       Ikless0 = zero;Ikgtr0  = zero;Ikdiag  = zero
#ifdef _mix
       Iklceil0= zero
#endif
       call get_Ikcollision(istep)
       Ikless0  = Ikless
       Ikgtr0   = Ikgtr
       Ikdiag   = Ikgtr(istep)-Ikless(istep)
#ifdef _mix
       Iklceil0 = Iklceil
#endif
    elseif(ips==2)then
       !Second Pass: get collision integrals up to t=T+\Delta=istep+1
       call get_Ikcollision(istep+1)
       Ikless   = (Ikless  + Ikless0)/2.d0
       Ikgtr    = (Ikgtr   + Ikgtr0)/2.d0
       Ikdiag   = (Ikgtr(istep+1)-Ikless(istep+1) + Ikdiag)/2.d0
#ifdef _mix
       Iklceil = (Iklceil + Iklceil0)/2.d0
#endif
    endif

    !Evolve the solution of KB equations for all the k-points:
    forall(it=0:istep)
       Gkless(it,istep+1) = Gkless(it,istep)*conjg(Udelta(ik,istep))+Ikless(it)*conjg(Vdelta(ik,istep))
       Gkgtr(istep+1,it)  = Gkgtr(istep,it)*Udelta(ik,istep)+Ikgtr(it)*Vdelta(ik,istep)
    end forall
    Gkgtr(istep+1,istep)=(Gkless(istep,istep)-xi)*Udelta(ik,istep)+Ikgtr(istep)*Vdelta(ik,istep)
    Gkless(istep+1,istep+1)=Gkless(istep,istep)-xi*dt*Ikdiag
    Gkgtr(istep+1,istep+1) = Gkless(istep+1,istep+1)-xi

#ifdef _mix
    forall(itau=0:Ltau)Gklceil(istep+1,itau)=Gklceil(istep,itau)*Udelta(ik,istep)+Iklceil(itau)*Vdelta(ik,istep)
    forall(itau=0:Ltau)Gkrceil(itau,1:istep+1)=conjg(Gklceil(1:istep+1,Ltau-itau))
#endif

    !===================================================================
    !This is apparently the bottle-neck of this algorithm:
    !$OMP PARALLEL PRIVATE(i,j)
    !$OMP DO
    do i=0,istep+1!nstep
       do j=0,istep+1!nstep
          if(i>j)Gkless(i,j)=-conjg(Gkless(j,i))
          if(i<j)Gkgtr(i,j)=-conjg(Gkgtr(j,i)) 
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    return
  end subroutine GFstep
  !******************************************************************
  !******************************************************************
  !******************************************************************









  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENTS : VER3.0_Apr2010. This version is based on the homogeneity
  ! with repects to k' of the matrices F_{k,k'}(t,t'), related to the 
  ! constant expression of Vhyb{k,k'}\== Vpd.
  ! A more general form has been tested in the past, using direct
  ! expression of the matrices in {k,k'}. See repo.
  !+-------------------------------------------------------------------+
  subroutine get_Ikcollision(Nt)
    integer,intent(in)          :: Nt
    integer                     :: i,itau,it,itp
    complex(8)                  :: I1,I2,Ib
    complex(8),dimension(0:Nt)  :: Vadv,Vret
    complex(8),dimension(0:Nt)  :: Fadv,Fret
#ifdef _mix
    complex(8),dimension(0:Ltau):: Ftau
#endif

    Ikless=zero; Ikgtr=zero
#ifdef _mix
    Iklceil=zero
#endif
    if(Nt==0)return

    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=0,Nt
       Vret(i) = SretF(Nt,i)        !+S0retF(Nt-i)
       Vadv(i) = conjg(SretF(Nt,i)) !+ conjg(S0retF(Nt-i))!conjg(Vret(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !Get I^<(t=it,t'=Nt) it:0,...,Nt == t=0,...,T+\Delta
    !I1 = \int_0^{t=it}  G^R*Sigma^<
    !I2 = \int_0^{t'=Nt} G^<*Sigma^A
    !Ib = \int_0^\beta   G^\rceil*S^\rceil
    !===============================================================x
    do it=0,Nt
       I1=zero ; I2=zero
       forall(i=0:Nt)Fret(i) = GkretF(it,i)
       I1 = integrate(Fret(0:it)*Sless(0:it,Nt),it,dt)      
       I2 = integrate(Gkless(it,0:Nt)*Vadv(0:Nt),Nt,dt)
       Ikless(it)=I1 + I2

#ifdef _mix
       Ib=zero
       Ib=integrate(Gklceil(it,0:Ltau)*Srceil(0:Ltau,Nt),Ltau,dtau)
       Ikless(it)=Ikless(it) + Ib/xi
#endif
    enddo

    !Get I^>(t=Nt,t'=itp) itp:0,...,Nt == t' = 0,...,T+\Delta=Nt
    !I1 = \int_0^{t=Nt}   S^R*G^>
    !I2 = \int_0^{t'=itp} S^>*G^A
    !Ib = \int_0^\beta    S^\lceil*G^\rceil
    !===============================================================
    do itp=0,Nt
       I1=zero ; I2=zero
       forall(i=0:Nt)Fadv(i) = conjg(GkretF(itp,i))
       I1=integrate(Vret(0:Nt)*Gkgtr(0:Nt,itp),Nt,dt)
       I2=integrate(Sgtr(Nt,0:itp)*Fadv(0:itp),itp,dt)
       Ikgtr(itp)=I1 + I2

#ifdef _mix
       Ib=zero
       Ib=integrate(Slceil(Nt,0:Ltau)*Gkrceil(0:Ltau,itp),Ltau,dtau)
       Ikgtr(itp)=Ikgtr(itp) + Ib/xi
#endif 
    enddo


#ifdef _mix
    !Get I^\lceil(Nt,itau) itau:0,...,Ltau == tau=0,...,beta
    !I1 = \int_0^Nt    S^Ret*G^\rceil
    !Ib = \int_0^\beta S^\lceil* G^\Mats
    !===============================================================
    do itau=0,Ltau
       I1=zero ; Ib=zero
       forall(i=0:Ltau)Ftau(i)  = Gktau(i-itau)
       I1=integrate(Vret(0:Nt)*Gklceil(0:Nt,itau),Nt,dt)
       Ib=integrate(Slceil(Nt,0:Ltau)*Ftau(0:Ltau),Ltau,dtau)
       Iklceil(itau)=I1 + Ib/xi
    enddo
#endif

  contains
    pure function GkretF(i,j)      
      integer,intent(in) :: i,j
      complex(8)         :: GkretF
      GkretF = heaviside(t(i-j))*(Gkgtr(i,j)-Gkless(i,j))
    end function GkretF
    !-------------------------------------------------------!

    !-------------------------------------------------------!
    pure function S0retF(i)      
      integer,intent(in) :: i
      complex(8)         :: S0retF
      S0retF = heaviside(t(i))*(S0gtr(i)-S0less(i))
    end function S0retF
    !-------------------------------------------------------!

    !-------------------------------------------------------!
    pure function SretF(i,j)      
      integer,intent(in) :: i,j
      complex(8)         :: SretF
      SretF = heaviside(t(i-j))*(Sgtr(i,j)-Sless(i,j))
    end function SretF
    !-------------------------------------------------------!

    !-------------------------------------------------------!
    function integrate(F,n,h)
      integer    :: i,n
      complex(8) :: F(0:n),integrate,sum
      real(8)    :: h
      integrate = integrate_simple(F(0:n),n+1,h)
      !Simpson:
      ! if(n<3)then
      !    integrate = integrate_simple(F(0:n),n+1,h)
      ! else
      !    call simpsn(F(0:n),n+1,h,integrate)
      ! endif
    end function integrate
    !-------------------------------------------------------!

    !-------------------------------------------------------!
    function integrate_simple(F,n,h)
      integer    :: i,n
      complex(8) :: F(n),integrate_simple,sum
      real(8)    :: h
      sum=0.d0
      do i=1,n
         sum=sum+F(i)
      enddo
      integrate_simple=sum*h
    end function integrate_simple
    !-------------------------------------------------------!

    !-------------------------------------------------------!
    function integrate_trapzd(F,n,h)
      integer    :: i,n
      complex(8) :: F(0:n),integrate_trapzd,sum
      real(8)    :: h
      sum=0.d0
      sum=F(0)/2.d0
      do i=1,n-1
         sum=sum+F(i)
      enddo
      sum=sum+F(n)/2.d0
      integrate_trapzd=sum*h
    end function integrate_trapzd
    !-------------------------------------------------------!

    !-------------------------------------------------------!
    function integrate_simpson(F,n,h)
      integer    :: i,n
      complex(8) :: F(0:n),integrate_simpson,sum
      real(8)    :: h
      if(n<3)then
         sum=0.d0
         do i=0,n
            sum=sum+F(i)
         enddo
         integrate_simpson=sum*h
      else
         sum=F(0)+F(n)
         do i=1,n-1,2                !S=S+4*F(2j-1) forall j=1,n/2-1
            sum=sum+4.d0*F(i)
         enddo
         do i=2,n-2,2                !S=S+2*F(2j) forall j=1,n/2
            sum=sum+2.d0*F(i)
         enddo
         integrate_simpson=sum*h/3.d0
      endif
    end function integrate_simpson
    !-------------------------------------------------------!

    !-------------------------------------------------------!
    subroutine simpsn(y,ntab,h,result)
      integer :: ntab,i,n
      complex(8) :: result
      complex(8) :: y(ntab)
      complex(8) :: sum1
      real(8) :: del(3)
      real(8) :: f
      real(8) :: g(3)
      real(8) :: h
      real(8) :: pii(3)
      result = 0.d0
      if ( ntab <= 2 ) then
         write(*,'(A)') ""
         write(*,'(A)') "SIMPSN - Fatal error!"
         write(*,'(A,I8)') "NTAB < 2, NTAb = ", ntab
         stop
      end if
      if (mod(ntab,2)==0) then
         n= ntab - 1
      else
         n= ntab
      end if
      result= y(1) + y(n) + 4.d0*y(n-1)
      do i=2,n-2,2
         result= result + 4.d0*y(i) + 2.d0*y(i+1)
      end do
      result= h*result/3.d0
      if (mod(ntab,2) == 1) then
         return
      end if
      f      = h**3
      del(1) = h
      del(2) = -2.d0*h
      del(3) = h
      g(1)   = h
      g(2)   = 0.d0
      g(3)   = -h
      pii(1) = 0.d0
      pii(2) = -h**2
      pii(3) = 0.d0
      n = n-1
      sum1 = 0.d0
      do i = 1, 3
         sum1 = sum1 + y(n-1+i)*del(i)*(f/3.d0 - g(i)*0.5d0*h**2 + pii(i)*h)
      end do
      result= result + 0.5d0*sum1/h**3
      return
    end subroutine simpsn
    !-------------------------------------------------------!
  end subroutine get_Ikcollision
  !******************************************************************
  !******************************************************************
  !******************************************************************








  !+-------------------------------------------------------------------+
  !PROGRAM  : GFSTEP
  !TYPE     : Subroutine
  !PURPOSE  : Execute the 2pass procedure to solve KB equations:
  !+-------------------------------------------------------------------+
  subroutine print_out_Gloc()
    integer :: i,j
    real(8),dimension(0:nstep) :: nt
    forall(i=0:nstep,j=0:nstep)
       gtless(i-j) = locGless(i,j)
       gtgtr(i-j)  = locGgtr(i,j)
       gtret(i-j)  = heaviside(t(i-j))*(locGgtr(i,j)-locGless(i,j))
    end forall
    if(heaviside(0.d0)==1.d0)gtret(0)=gtret(0)/2.d0
    call cfft_rt2rw(gtret,gfret,nstep) ;    gfret=gfret*dt ; call swap_fftrt2rw(gfret)
    call splot("locGless_t.ipt",t,gtless,TT)
    call splot("locGgtr_t.ipt",t,gtgtr,TT)
    call splot("locGret_t.ipt",t,gtret,TT)
    call splot("locGret_realw.ipt",wr,gfret,TT)
    call splot("locDOS.ipt",wr,-aimag(gfret)/pi,append=TT)

#ifdef _mix
    call plot_3D("locGlceil3D","$\Delta t$","$\tau$","Z",t(0:nstep)/dt,tau,locGlceil(0:nstep,0:Ltau))
    call plot_3D("locGrceil3D","$\tau$","$\Delta t$","Z",tau,t(0:nstep)/dt,locGrceil(0:Ltau,0:nstep))
    call plot_3D("locGless3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,locGless)
    call plot_3D("locGgtr3D","X/$\Delta t$","Y/$\Delta t$","Z",t(0:nstep)/dt,t(0:nstep)/dt,locGgtr)
    forall(i=0:nstep)nt(i)=-xi*locGless(i,i)
    call splot("nVStime.ipt",t(0:nstep),2.d0*nt(0:nstep),append=TT)
#endif
    return
  end subroutine print_out_Gloc
  !******************************************************************
  !******************************************************************
  !******************************************************************








  !#####################################################################
  !INITIAL CONDITIONS:
  !#####################################################################
  !+-------------------------------------------------------------------+
  !PROGRAM  : BUILD_IC
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial condition for the solution of KB equations
  ! G_k^<(0,0)        = xi*G_k^M(tau=0-)
  ! G_k^>(0,0)        = xi*(G_k^<(0,0) - 1.0)
  ! G_k^\lceil(0,tau) = xi*G_k^M(tau<0) = -xi*G_k^M(beta-tau>0)
  ! G_k^\rceil(0,tau) = xi*G_k^M(tau>0)  
  !+-------------------------------------------------------------------+
  subroutine build_ic
    integer :: ik,i
    real(8) :: en,mu,invtemp
#ifdef _mix
    real(8)                   :: w
    complex(8)                :: zeta
    complex(8),dimension(2*L) :: Giw
    real(8),dimension(0:L)    :: Gt
    real(8),dimension(0:Ltau) :: Gtau
#endif

    call dump("Building initial conditions:")
    call system("if [ ! -d InitialConditions ]; then mkdir InitialConditions; fi")
    mu=xmu
    invtemp=beta
    if(iquench)mu=xmu0
    if(iquench)invtemp=beta0
    if(irdeq)then
       icGkless = xi*icnk
    else
       do ik=1,Lk
          en           = epsik(ik)-mu
          Icgkless(ik) = xi*fermi0(en,invtemp)
       enddo
    endif
    call splot("InitialConditions/icGklessVSepsik.ipt",epsik(1:Lk),aimag(icGkless(1:Lk)))

#ifdef _mix
    do ik=1,Lk
       call splot("InitialConditions/icnVSk.ipt",epsik(ik),-icGktau(ik,Ltau),append=TT)
       call splot("InitialConditions/icGklceilVSepsikVStau.ipt",tau(0:Ltau),icGktau(ik,0:Ltau),append=TT)
    enddo
#endif
    return
  end subroutine build_ic
  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine recover_ic(ik)
    integer :: ik,itau
    Gkless(0,0)=icGkless(ik)
    Gkgtr(0,0) =icGkless(ik)-xi
#ifdef _mix
    forall(itau=0:Ltau)
       Gklceil(0,itau) =-xi*icGktau(ik,Ltau-itau) !Gk^\lceil(0,tau) = xi*Gk^M(tau<0) = -xi*GM(beta-tau>0)
       Gkrceil(itau,0) = xi*icGktau(ik,itau)      !Gk^\rceil(0,tau) = xi*Gk^M(tau>0) 
       Gktau(itau)     = icGktau(ik,itau)
    end forall
    forall(itau=1:Ltau)Gktau(-itau)=-icGktau(ik,Ltau-itau)
#endif
    return
  end subroutine recover_ic
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine save_ic(ik)
    integer :: ik
    icGkless(ik)=Gkless(0,0)
  end subroutine save_ic
  !******************************************************************
  !******************************************************************
  !******************************************************************






  !#####################################################################
  !EVOLUTION OPERATORS U_\Delta, V_\Delta
  !#####################################################################
  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine buildUV
    integer :: ik,i
    forall(ik=1:Lk,i=0:nstep) 
       Udelta(ik,i)=UdeltaF(ik,i)
       Vdelta(ik,i)=VdeltaF(ik,i)
    end forall
  end subroutine buildUV
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  elemental function Hbar(ik,istep)
    integer,intent(in) :: ik,istep  
    integer      :: i,j
    complex(8)   :: Hbar
    real(8)      :: tbar
    type(vect2D) :: kt,Ak
    tbar=t(istep) + dt/2.d0
    i=ik2ix(ik)
    j=ik2iy(ik)
    Ak=Afield(tbar,Ek)
    kt=kgrid(i,j) - Ak !tbar*Ek
    Hbar=epsk(kt)
  end function Hbar
  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  elemental function UdeltaF(ik,istep) 
    integer,intent(in)    :: ik,istep
    complex(8) :: UdeltaF
    complex(8) :: arg
    arg=Hbar(ik,istep)
    UdeltaF=exp(-xi*arg*dt)
  end function UdeltaF
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  elemental function VdeltaF(ik,istep)
    integer,intent(in)    :: ik,istep
    complex(8) :: VdeltaF
    complex(8) :: arg
    arg=Hbar(ik,istep)
    VdeltaF=exp(-xi*arg*dt)
    if(abs(arg*dt) <= 1.d-5)then
       VdeltaF=-xi*dt
    else
       VdeltaF=(VdeltaF-1.d0)/arg
    endif
  end function VdeltaF
  !******************************************************************
  !******************************************************************
  !******************************************************************



end module KADANOFBAYM
