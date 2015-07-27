program neqDMFT
  USE NEQ_DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                 :: i,j,ik,itime,iloop,ix,iy,Lk
  logical                                 :: converged
  real(8)                                 :: wband
  character(len=16)                       :: finput
  type(kb_contour_gf),dimension(2,2)      :: Gloc
  type(kb_contour_gf),dimension(2,2)      :: Gwf
  type(kb_contour_sigma),dimension(2,2)   :: Sigma
  !
  type(kb_contour_gf)                     :: G0aux(2)
  type(kb_contour_gf),allocatable         :: Gk(:,:)
  type(kb_contour_gf),allocatable         :: Gk_aux(:,:)
  type(kb_contour_dgf),allocatable        :: dGk_aux(:,:)
  type(kb_contour_dgf),allocatable        :: dGk_aux_old(:,:)
  !
  type(kb_contour_gf)                     :: Kerk,Ker1,Ker2,Kfoo(7)
  !
  complex(8),dimension(:,:),allocatable   :: Ham
  real(8),dimension(:),allocatable        :: Delta,nt
  !RESULTS:
  real(8),dimension(:),allocatable        :: epsik,wt
  real(8),dimension(:,:),allocatable      :: nk
  real(8)                                 :: D,de,intwt

  !READ THE INPUT FILE (in vars_global):
  call parse_cmd_variable(finput,"FINPUT",default='inputNEQ.conf')
  call parse_input_variable(wband,"wband",finput,default=2d0,comment="half-bandwidth W=2t")
  call parse_input_variable(Lk,"Lk",finput,default=100,comment="Number of energy levels for Bethe DOS")
  call read_input_init(trim(finput))



  !BUILD TIME GRIDS AND NEQ-PARAMETERS:
  call allocate_kb_contour_params(cc_params,Ntime,Ntau,Niw)
  call setup_kb_contour_params(cc_params,dt,beta)


  !BUILD THE LATTICE STRUCTURE (in lib/square_lattice):
  allocate(epsik(Lk),wt(Lk))
  epsik = linspace(-wband,wband,Lk,mesh=de)
  do i=1,Lk
     wt(i) = dens_bethe(epsik(i),wband)
  enddo
  call splot("DOSbethe.ipt",epsik,wt)
  wt=wt/sum(wt)



  !ALLOCATE ALL THE FUNCTIONS INVOLVED IN THE CALCULATION:
  call allocate_kb_contour_sigma(Sigma,cc_params) !Self-Energy function
  call allocate_kb_contour_gf(Gloc,cc_params)     !Local Green's function
  call allocate_kb_contour_gf(Gwf,cc_params)      !Local Weiss-Field function
  call allocate_kb_contour_gf(G0aux,cc_params)
  allocate(Gk(2,Lk))
  allocate(Gk_aux(2,Lk))
  allocate(dGk_aux(2,Lk))
  allocate(dGk_aux_old(2,Lk))
  do ik=1,Lk
     call allocate_kb_contour_gf(Gk(:,ik),cc_params)
     call allocate_kb_contour_gf(Gk_aux(:,ik),cc_params)
     call allocate_kb_contour_dgf(dGk_aux(:,ik),cc_params)
     call allocate_kb_contour_dgf(dGk_aux_old(:,ik),cc_params)
  end do
  !
  call allocate_kb_contour_gf(Kerk,cc_params)
  call allocate_kb_contour_gf(Ker1,cc_params)
  call allocate_kb_contour_gf(Ker2,cc_params)
  do i=1,7
     call allocate_kb_contour_gf(Kfoo(i),cc_params)
  enddo
  !
  allocate(ham(cc_params%Ntime,Lk))
  allocate(nk(cc_params%Ntime,Lk))
  allocate(Delta(cc_params%Ntime))
  allocate(nt(cc_params%Ntime))


  !READ OR GUESS THE INITIAL WEISS FIELD G0 (t=t'=0)
  cc_params%Nt=1
  Gloc=zero
  call neq_continue_equilibirum(Gwf,Gk,Gk_aux,dGk_aux,Gloc,Sigma,epsik,wt,cc_params)
  Nt(1)    = measure_dens(Gloc(1,1),cc_params)
  Delta(1) = measure_delta(Gloc(1,2),cc_params)
  print*,Nt(1),Delta(1)
  do ik=1,Lk
     nk(1,ik)=dimag(Gk(1,ik)%less(1,1))
  enddo
  do i=1,cc_params%Ntime
     do ik=1,Lk
        ham(i,ik)=epsik(ik)
     enddo
  enddo
  call measure_observables(Gloc,Sigma,Gk,Ham,Wt,cc_params)
  !<DEBUG comment
  call plot_kb_contour_gf("Sigma",Sigma(1,1),cc_params)
  call plot_kb_contour_gf("Self",Sigma(1,2),cc_params)
  call plot_kb_contour_gf("Gloc",Gloc(1,1),cc_params)
  call plot_kb_contour_gf("Floc",Gloc(1,2),cc_params)
  call plot_kb_contour_gf("G0",Gwf(1,1),cc_params)
  call plot_kb_contour_gf("F0",Gwf(1,2),cc_params)
  !>DEBUG


  !START THE TIME_STEP LOOP  1<t<=Nt
  !AT EACH TIME_STEP PERFORM A FULL DMFT CALCULATION:
  do itime=2,cc_params%Ntime
     print*,""
     print*,"time step=",itime
     cc_params%Nt=itime
     !
     !
     !prepare the weiss-field at this actual time_step for DMFT:
     do i=1,2
        do j=1,2
           call extrapolate_kb_contour_gf(Gwf(i,j),cc_params)
           call extrapolate_kb_contour_gf(Gloc(i,j),cc_params)
        enddo
     enddo
     do i=1,2
        do ik=1,Lk
           dGk_aux_old(i,ik) = dGk_aux(i,ik)
        enddo
     enddo
     !
     !
     !
     !DMFT LOOP:
     iloop=0;converged=.false.
     do while(.not.converged.AND.iloop<nloop)
        iloop=iloop+1
        write(*,"(A,I4,A1)",advance='no')"dmft loop=",iloop," "


        !IMPURITY SOLVER: IPT.
        !call neq_solve_ipt(Gwf,Gloc,Sigma,cc_params)
        !<DEBUG We use the 1st order approximation for debug: Sigma%reg=0,Sigma%hf!=0
        nt(itime)   =measure_dens(Gloc(1,1),cc_params)
        delta(itime)=measure_delta(Gloc(1,2),cc_params)
        print*,nt(itime),delta(itime)
        Sigma(1,1)%hf(itime)=zero
        Sigma(1,2)%hf(itime)=-delta(itime)
        Sigma(2,1)%hf(itime)=-delta(itime)
        Sigma(2,2)%hf(itime)=zero
        do i=1,2
           do j=1,2
              Sigma%reg(i,j)=zero
           enddo
        enddo
        !>DEBUG


        !
        !PERFORM THE SELF_CONSISTENCY: get local GF + update Weiss-Field
        !solve KBE for the auxiliary functions:
        do ik=1,Lk
           ! solve a VIDE for g11k & f21k:
           !  - id/dt g11k = delta_cc + hk_HF_11*g11k + K*g11k; K=Sigma%reg(1,1)
           !  - id/dt f21k = delta_cc + hk_HF_22*f21k + K*f21k; K=Sigma%reg(2,2)
           call vide_kb_contour_gf( Ham(:,ik),Sigma(1,1),Gk_aux(1,ik),dGk_aux_old(1,ik),dGk_aux(1,ik),cc_params)
           call vide_kb_contour_gf(-Ham(:,ik),Sigma(2,2),Gk_aux(2,ik),dGk_aux_old(2,ik),dGk_aux(2,ik),cc_params)
           ! solve a VIE for Gk(1,1) and a convolution to get G(2,1;k)
           !build Ker = gk*(H^1{12}δ + S(1,2))*fbark*(H^{21}δ + S(2,1))
           call convolute_kb_contour_gf(      &
                G=[Gk_aux(1,ik),Gk_aux(2,ik)],&
                S=[Sigma(1,2),Sigma(2,1)],    &
                mask=[0,1,0,1],P=Kerk,cc_params)
           ! G(1,1;k) = gk + Ker*G(1,1;k) 
           call vie_kb_contour_gf(gk_aux(1,ik),Kerk,Gk(1,ik),cc_params)
           ! G(2,1;k) = fbark*Sigma(2,1)*G(1,1;k)
           call convolute_kb_contour_gf(  &
                G=[Gk_aux(2,ik),Gk(1,ik)],&
                S=[Sigma(2,1)],           &
                mask=[0,1,0],P=Gk(2,ik),params=cc_params)
        enddo

        call sum_kb_contour_gf(Gk(1,:),wt(:),Gloc(1,1),cc_params)
        call sum_kb_contour_gf(Gk(2,:),wt(:),Gloc(2,1),cc_params)
        call get_bar(Gloc(2,2),Gloc(1,1),cc_params)
        call get_bar(Gloc(1,2),Gloc(2,1),cc_params)
        !
        !
        !
        ! SOLVE FOR THE WEISS FIELDS: G0(1,1) & G0(2,1):
        ! solve the VIE for auxiliary weiss fields g0,fbar0
        ! g0    = 11 - K*g0    ; K=G(1,1)*Sigma(1,1) + G(1,2)*Sigma(2,1)
        ! fbar0 = 11 - K*fbar0 ; K=G(2,2)*Sigma(2,2) + G(2,1)*Sigma(1,2)
        call convolute_kb_contour_gf(Gloc(1,1),Sigma(1,1)%reg,Kfoo(1),cc_params)
        call convolute_kb_contour_gf(Gloc(1,2),Sigma(2,1)%reg,Kfoo(2),cc_params)
        call sum_kb_contour_gf(Kfoo(1),1d0,Kfoo(2),1d0,Kerk,cc_params)
        call vie_kb_contour_gf(Kerk,g0aux(1),cc_params)
        !
        call convolute_kb_contour_gf(Gloc(2,2),Sigma(2,2)%reg,Kfoo(1),cc_params)
        call convolute_kb_contour_gf(Gloc(2,1),Sigma(1,2)%reg,Kfoo(2),cc_params)
        call sum_kb_contour_gf(Ker0(1),1d0,Ker0(2),1d0,Kerk,cc_params)
        call vie_kb_contour_gf(Kerk,g0aux(2),cc_params)
        !
        ! Solve VIE for G0(1,1): G0(1,1) = K1 + K2*G0(1,1)
        ! K1 = g0*G - g0*G*Sreg*fbar0*Fbar - g0*F*Sigmabar*fbar0*Fbar
        ! K1 = Kfoo(1) + Kfoo(2) + Kfoo(3)
        call convolute_kb_contour_gf( g0aux(1),Gloc(1,1),Kfoo(1),cc_params)
        call convolute_kb_contour_gf( [Kfoo(1),Sigma(1,2)%reg,G0aux(2),Gloc(2,1)]           , Kfoo(2), cc_params)
        call convolute_kb_contour_gf( [G0aux(1),Gloc(1,2),Sigma(2,2)%reg,G0aux(2),Gloc(2,1)], Kfoo(3), cc_params)
        call sum_kb_contour_gf(Kfoo(1:3),[1d0,-1d0,-1d0],Ker1,cc_params)
        !
        ! K2 = g0*G*S*fbar0*Fbar*Sigma + g0*F*Sigmabar*fbar0*Fbar*Sigma +  g0*G*S*fbar0*Gbar*Sbar    +  g0*F*Sigmabar*fbar0*Gbar*Sbar
        ! K2 = Kfoo(2)*Sigma           + Kfoo(3)*Sigma                  +  Kfoo(1)*S*fbar0*Gbar*Sbar +  g0*F*Sigmabar*fbar0*Gbar*Sbar
        ! K2 = Kfoo(4)                 + Kfoo(5)                        +  Kfoo(6)                   +  Kfoo(7)        
        call convolute_kb_contour_gf(Kfoo(2),Sigma(1,1)%reg,Kfoo(4),cc_params)
        call convolute_kb_contour_gf(Kfoo(3),Sigma(1,1)%reg,Kfoo(5),cc_params)
        call convolute_kb_contour_gf([Kfoo(1),Sigma(1,2)%reg,G0aux(2),Gloc(2,2),Sigma(2,1)%reg]           ,Kfoo(6),cc_params)
        call convolute_kb_contour_gf([g0aux(1),Gloc(1,2),Sigma(2,2)%reg,g0aux(2),Gloc(2,2),Sigma(2,1)%reg],Kfoo(7),cc_params)
        call sum_kb_contour_gf(Ker2(4:7),[1d0,1d0,1d0,1d0],Ker2,cc_params)
        !
        call vie_kb_contour_gf(Ker1,Ker2,Gwf(1,1),cc_params)
        !
        ! Solve for G0(2,1) : 
        ! G0(2,1)=F0bar= f0bar*Fbar - f0bar*Fbar*Sigma*G0 - f0bar*Gbar*Sbar*G0
        ! G0(2,1)=F0bar= Ker1(1)    - Ker1(2)             - Ker1(3)
        call convolute_kb_contour_gf(G0aux(2),Gloc(2,1),Kfoo(1),cc_params)
        call convolute_kb_contour_gf( [G0aux(2),Gloc(2,1),Sigma(1,1)%reg,Gwf(1,1)] ,Kfoo(2),cc_params)
        call convolute_kb_contour_gf( [G0aux(2),Gloc(2,2),Sigma(2,1)%reg,Gwf(1,1)] ,Kfoo(3),cc_params)
        call sum_kb_contour_gf( Kfoo(1:3), [1d0,-1d0,-1d0], Gwf(2,1), cc_params)
        !
        ! get the other components using symmetries:
        call get_bar(Gwf(2,2),Gwf(1,1),cc_params)
        call get_bar(Gwf(1,2),Gwf(2,1),cc_params)
        !
        !
        !
        !CHECK CONVERGENCE
        !===================================================================================================================!
        converged = convergence_check(Gwf(1,1),cc_params)
     enddo
     !
     !
     !
     !EVALUATE AND PRINT THE RESULTS OF THE CALCULATION
     !===================================================================================================================!
     call measure_observables(Gloc,Sigma,Gk,Ham,Wt,cc_params)
     forall(ik=1:Lk)nk(itime,ik)=dimag(Gk(1,ik)%less(itime,itime))
  enddo
  !
  !
  !EVALUATE AND PRINT OTHER RESULTS OF THE CALCULATION
  call splot3d("nkVSepsikVStime.nipt",cc_params%t,epsik,nk)
  call plot_kb_contour_gf("Sigma",Sigma(1,1),cc_params)
  call plot_kb_contour_gf("Self",Sigma(1,2),cc_params)
  call plot_kb_contour_gf("Gloc",Gloc(1,1),cc_params)
  call plot_kb_contour_gf("Floc",Gloc(1,2),cc_params)
  call plot_kb_contour_gf("G0",Gwf(1,1),cc_params)
  call plot_kb_contour_gf("F0",Gwf(1,2),cc_params)
  print*,"BRAVO"


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : check convergence of the calculation:
  !+-------------------------------------------------------------------+
  function convergence_check(G,params) result(converged)
    type(kb_contour_gf)                 :: G
    type(kb_contour_params)             :: params
    logical                             :: converged
    complex(8),dimension(:),allocatable :: test_func
    integer                             :: i,N,L,Ntot
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    Ntot=2*N+L+1
    allocate(test_func(Ntot))
    test_func=zero
    do i=1,N
       test_func(i)  = G%ret(N,i)
       test_func(N+i)= G%less(N,i)
    enddo
    do i=0,L
       test_func(2*N+i+1)=G%lmix(N,i)
    enddo
    converged=check_convergence(test_func,dmft_error,Nsuccess,nloop)
    deallocate(test_func)
    !if(isnan(err))stop "Aborted convergence: error=NaN"
  end function convergence_check

end PROGRAM neqDMFT
