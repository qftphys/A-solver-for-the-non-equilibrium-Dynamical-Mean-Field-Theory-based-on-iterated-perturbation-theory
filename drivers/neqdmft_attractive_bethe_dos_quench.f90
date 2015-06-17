program neqDMFT
  USE NEQ_DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                 :: i,j,ik,itime,iloop,ix,iy,Lk
  logical                                 :: converged
  real(8)                                 :: wband
  character(len=16)                       :: finput
  type(kb_contour_gf),dimension(2,2)      :: Gloc,Gwf
  type(kb_contour_gf),dimension(2,2)      :: Sigma
  type(kb_contour_gf),dimension(2,2)      :: SigmaReg
  complex(8),dimension(:,:,:),allocatable :: SigmaHFB
  !
  type(kb_contour_gf)                     :: Kerk,Ker0(2),Ker1(3),Ker2(4),Ker(2)
  type(kb_contour_gf)                     :: G0aux(2)
  type(kb_contour_gf),allocatable         :: Gk(:,:),Gk_aux(:,:)
  type(kb_contour_dgf),allocatable        :: dGk_aux(:,:),dGk_aux_old(:,:)
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
  !=====================================================================
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
  allocate(SigmaHFB(2,2,cc_params%Ntime)) !Self-Energy function (Hartree-Fock-Bogoliubov term)
  call allocate_kb_contour_gf(SigmaReg,cc_params) !Self-Energy function (Regular or 2nd-order term)
  call allocate_kb_contour_gf(Sigma,cc_params) !Self-Energy function
  call allocate_kb_contour_gf(Gloc,cc_params)  !Local Green's function
  call allocate_kb_contour_gf(Gwf,cc_params)   !Local Weiss-Field function
  call allocate_kb_contour_gf(G0aux,cc_params)
  call allocate_kb_contour_gf(Kerk,cc_params)
  do i=1,2
     call allocate_kb_contour_gf(Ker(i),cc_params)
     call allocate_kb_contour_gf(Ker0(i),cc_params)
  enddo
  do i=1,3
     call allocate_kb_contour_gf(Ker1(i),cc_params)
  enddo
  do i=1,4
     call allocate_kb_contour_gf(Ker2(i),cc_params)
  enddo
  !
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
  allocate(ham(cc_params%Ntime,Lk))
  allocate(nk(cc_params%Ntime,Lk))
  allocate(Delta(cc_params%Ntime))
  allocate(nt(cc_params%Ntime))


  !READ OR GUESS THE INITIAL WEISS FIELD G0 (t=t'=0)
  cc_params%Nt=1
  call neq_continue_equilibirum(Gwf,Gk,Gk_aux,dGk_aux,Gloc,SigmaHFB,SigmaReg,Sigma,epsik,wt,cc_params)
  call measure_observables(Gloc,Sigma,Gk,Ham,Wt,cc_params)
  Delta(1) = measure_delta(Gloc(1,2),cc_params)
  do ik=1,Lk
     nk(1,ik)=dimag(Gk(1,ik)%less(1,1))
  enddo
  do i=1,cc_params%Ntime
     do ik=1,Lk
        ham(i,ik)=epsik(ik)
     enddo
  enddo



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
        !
        !
        !
        !IMPURITY SOLVER: IPT.
        !call neq_solve_hfb(Gloc,SigmaHFB,cc_params)
        nt(itime)   =measure_dens(Gloc(1,1),cc_params)
        delta(itime)=measure_delta(Gloc(1,2),cc_params)
        SigmaHFB(1,1,itime)=zero
        SigmaHFB(1,2,itime)=-delta(itime)
        SigmaHFB(2,1,itime)=-delta(itime)
        SigmaHFB(2,2,itime)=zero
        !get regular/2nd order term of sigma from the weiss fields
        !call neq_solve_ipt(Gwf,SigmaReg,cc_params)
        SigmaReg=zero
        call sum_kb_contour_gf(SigmaReg(1,1),1d0,SigmaHFB(1,1,:),1d0,Sigma(1,1),cc_params)
        call sum_kb_contour_gf(SigmaReg(1,2),1d0,SigmaHFB(1,2,:),1d0,Sigma(1,2),cc_params)
        call sum_kb_contour_gf(SigmaReg(2,1),1d0,SigmaHFB(2,1,:),1d0,Sigma(2,1),cc_params)
        call sum_kb_contour_gf(SigmaReg(2,2),1d0,SigmaHFB(2,2,:),1d0,Sigma(2,2),cc_params)
        !
        !
        !
        !PERFORM THE SELF_CONSISTENCY: get local GF + update Weiss-Field
        !solve KBE for the auxiliary functions:
        do ik=1,Lk
           ! solve a VIDE for g11k: id/dt gk = delta_cc + hk_11*gk + K*gk; K=Sigma(1,1)
           ! solve a VIDE for f21k=>g22k: id/dt fbark = delta_cc + hk_12*fbark + K*fbark; K=Sigma(2,2)
           call vide_kb_contour_gf( Ham(:,ik)+SigmaHFB(1,1,:),SigmaReg(1,1),Gk_aux(1,ik),dGk_aux_old(1,ik),dGk_aux(1,ik),cc_params)
           call vide_kb_contour_gf(-Ham(:,ik)+SigmaHFB(2,2,:),SigmaReg(2,2),Gk_aux(2,ik),dGk_aux_old(2,ik),dGk_aux(2,ik),cc_params)
           ! solve a VIE for G(1,1;k) and a convolution to get G(2,1;k)
           ! Ker      = gk * Sigma(1,2) * fbark * Sigma(2,1)
           ! G(1,1;k) = gk + Ker*G(1,1;k) 
           ! G(2,1;k) = fbark * Sigma(2,1) * G(1,1;k)
           call convolute_kb_contour_gf( [Gk_aux(1,ik),Sigma(1,2),Gk_aux(2,ik),Sigma(2,1)],Kerk,cc_params)
           call vie_kb_contour_gf(gk_aux(1,ik),Kerk,Gk(1,ik),cc_params)
           call convolute_kb_contour_gf( [Gk_aux(2,ik),Sigma(2,1),Gk(1,ik)] ,Gk(2,ik),cc_params)
        enddo
        call sum_kb_contour_gf(Gk(1,:),wt(:),Gloc(1,1),cc_params)
        call sum_kb_contour_gf(Gk(2,:),wt(:),Gloc(2,1),cc_params)
        call get_bar(Gloc(2,2),Gloc(1,1),cc_params)
        call get_bar(Gloc(1,2),Gloc(2,1),cc_params)
        !
        !
        !
        ! SOLVE FOR THE WEISS FIELDS: G0(1,1) & G0(2,1):
        !===================================================================================================================!
        ! solve the VIE for auxiliary weiss fields g0,fbar0
        ! g0    = 11 - K*g0    ; K=G(1,1)*SigmaReg(1,1) + G(1,2)*SigmaReg(2,1)
        ! fbar0 = 11 - K*fbar0 ; K=G(2,2)*SigmaReg(2,2) + G(2,1)*SigmaReg(1,2)
        call convolute_kb_contour_gf(Gloc(1,1),SigmaReg(1,1),Ker0(1),cc_params)
        call convolute_kb_contour_gf(Gloc(1,2),SigmaReg(2,1),Ker0(2),cc_params)
        call sum_kb_contour_gf(Ker0(1),1d0,Ker0(2),1d0,Kerk,cc_params)
        call vie_kb_contour_gf(Kerk,g0aux(1),cc_params)
        !
        call convolute_kb_contour_gf(Gloc(2,2),SigmaReg(2,2),Ker0(1),cc_params)
        call convolute_kb_contour_gf(Gloc(2,1),SigmaReg(1,2),Ker0(2),cc_params)
        call sum_kb_contour_gf(Ker0(1),1d0,Ker0(2),1d0,Kerk,cc_params)
        call vie_kb_contour_gf(Kerk,g0aux(2),cc_params)
        !
        ! Solve VIE for G0(1,1): G0(1,1) = K1 + K2*G0(1,1)
        ! K1 = g0*G - g0*G*Sreg*fbar0*Fbar - g0*F*SigmaRegbar*fbar0*Fbar
        ! K1 = Ker1(1) + Ker1(2) + Ker1(3)
        call convolute_kb_contour_gf( g0aux(1),Gloc(1,1),Ker1(1),cc_params)
        call convolute_kb_contour_gf( [Ker1(1),SigmaReg(1,2),G0aux(2),Gloc(2,1)]           , Ker1(2), cc_params)
        call convolute_kb_contour_gf( [G0aux(1),Gloc(1,2),SigmaReg(2,2),G0aux(2),Gloc(2,1)], Ker1(3), cc_params)
        call sum_kb_contour_gf(Ker1(1:3),[1d0,-1d0,-1d0],Ker(1),cc_params)
        !
        ! K2 = g0*G*S*fbar0*Fbar*Sigma + g0*F*Sigmabar*fbar0*Fbar*Sigma +  g0*G*S*fbar0*Gbar*Sbar    +  g0*F*Sigmabar*fbar0*Gbar*Sbar
        ! K2 = Ker2(1)                 + Ker2(2)                        +  Ker2(3)                   +  Ker2(4)
        ! K2 = Ker1(2)*Sigma           + Ker1(3)*Sigma                  +  Ker1(1)*S*fbar0*Gbar*Sbar +  g0*F*Sigmabar*fbar0*Gbar*Sbar
        ! K2 = K_gamma(5) + K_delta(5) + K_epsi(1)*Sbar + K_zeta(1)*Sbar
        ! K2 = K_gamma(5) + K_delta(5) + K_epsi(2) + K_zeta(2)
        call convolute_kb_contour_gf(Ker1(2),Sigma(1,1),Ker2(1),cc_params)
        call convolute_kb_contour_gf(Ker1(3),Sigma(1,1),Ker2(2),cc_params)
        call convolute_kb_contour_gf([Ker1(1),SigmaReg(1,2),G0aux(2),Gloc(2,2),Sigma(2,1)]           ,Ker2(3),cc_params)
        call convolute_kb_contour_gf([g0aux(1),Gloc(1,2),SigmaReg(2,2),g0aux(2),Gloc(2,2),Sigma(2,1)],Ker2(4),cc_params)
        call sum_kb_contour_gf(Ker2(1:4),[1d0,1d0,1d0,1d0],Ker(2),cc_params)
        !
        call vie_kb_contour_gf(Ker(1),Ker(2),Gwf(1,1),cc_params)
        !
        ! Solve for G0(2,1) : 
        ! G0(2,1)=F0bar= f0bar*Fbar - f0bar*Fbar*Sigma*G0 - f0bar*Gbar*Sbar*G0
        ! G0(2,1)=F0bar= Ker1(1)    - Ker1(2)             - Ker1(3)
        call convolute_kb_contour_gf(G0aux(2),Gloc(2,1),Ker1(1),cc_params)
        call convolute_kb_contour_gf( [G0aux(2),Gloc(2,1),Sigma(1,1),Gwf(1,1)] ,Ker1(2),cc_params)
        call convolute_kb_contour_gf( [G0aux(2),Gloc(2,2),Sigma(2,1),Gwf(1,1)] ,Ker1(3),cc_params)
        call sum_kb_contour_gf( Ker1(1:3), [1d0,-1d0,-1d0], Gwf(2,1), cc_params)
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
  call plot_kb_contour_gf("Sigma.nipt",Sigma(1,1),cc_params)
  call plot_kb_contour_gf("Gloc.nipt",Gloc(1,1),cc_params)
  call plot_kb_contour_gf("G0.nipt",Gwf(1,1),cc_params)
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
