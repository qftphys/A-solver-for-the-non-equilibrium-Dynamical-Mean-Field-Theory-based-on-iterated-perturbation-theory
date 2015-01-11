!###############################################################
!     PROGRAM  : FUNCS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!     LAST UPDATE: 07/2009
!###############################################################
module FUNX_NEQ
  !LIBRARY:  
  USE TOOLS
  USE GRIDS
  USE LATTICE
  USE PROGRESS
  USE DLPLOT
  USE FFTW
  USE MPI
  !LOCAL:
  use VARS_GLOBAL
  implicit none
  private
  public allocate_locFunx,&
       Xc,                &
       getG0guess,        &
       Build_bath_t,      &
       Build_bath_w,      &
       dysonG0,           &
       update_sigma
  save
contains
  !+----------------------------------------------------------------+
  !PROGRAM  : MASSIVE_ALLOC
  !TYPE     : subroutine
  !PURPOSE  : massive allocation of work array
  !+----------------------------------------------------------------+
  subroutine allocate_locFunx(char)
    character(len=1) :: char
    if(char=='a')then
       if(mpiID==0)call dump("Allocating local-Funcs:")
       allocate(G0gtr(0:nstep,0:nstep),G0less(0:nstep,0:nstep))!Weiss Fields      
       allocate(S0gtr(-nstep:nstep),S0less(-nstep:nstep),S0ret(-nstep:nstep))
       allocate(Sgtr(0:nstep,0:nstep),Sless(0:nstep,0:nstep),Sret(0:nstep,0:nstep)) 
       allocate(icGkless(Lk))                                  !Initial conditions:
       !allocate(Gmktau(Lk,-Ltau:Ltau))
       allocate(locGless(0:nstep,0:nstep),locGgtr(0:nstep,0:nstep))!Local GF
       if(mpiID==0)call dump("done")
       if(mpiID==0)call dump("")

    elseif(char=='d')then
       if(mpiID==0)call dump("Deallocating:")
       deallocate(G0gtr,G0less)
       deallocate(S0less,S0gtr,S0ret)
       deallocate(Sgtr,Sless)
       deallocate(icGkless)
       !deallocate(Gmktau)
       deallocate(locGless,locGgtr)
       if(mpiID==0)call dump("done")
       if(mpiID==0)call dump("")
    end if
  end subroutine allocate_locFunx
  !******************************************************************
  !******************************************************************
  !******************************************************************








  !+------------------------------------------------------------------+
  !PROGRAM  : Xc
  !TYPE     : function
  !PURPOSE  : build  the function X^c(e,t,t') equal to the term in front
  !of the exponential in the expressions for G_0(t,t') 
  !COMMENT  : The function is defined WITHOUT the imaginary unit xi
  !+------------------------------------------------------------------+
  function Xc(cc,x)
    real(8)          :: Xc,x
    character(len=1) :: cc
    Xc=0.d0    
    if(cc == '<')then
       Xc=fermi(x,beta)
    elseif(cc == '>')then
       Xc=fermi(x,beta)-1.d0
    endif
  end function Xc
  !********************************************************************
  !********************************************************************
  !********************************************************************













  !+-------------------------------------------------------------------+
  !PROGRAM  : GETG0GUESS
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial Weiss Fields G0^{<,>,R,\lceil} used to obtain
  !the self-energies functions (cf update_sigma below).
  !The functions are build starting fomr the Bath DOS obtained 
  !using the Keldysh solution of the associated equilibrium problem.
  !+-------------------------------------------------------------------+
  subroutine getG0guess()
    integer :: i,j,ik,itau,jtau,istep,it,m
    real(8) :: w,en,intE,ttau,dos
    complex(8) :: peso,A,nt,An
    real(8) :: rho,nfp,nfm
    complex(8),allocatable,dimension(:,:) :: dummy
    complex(8),dimension(-L:L):: g0tgtr,g0tless,g0tdummy
    if(mpiID==0)call dump("Get G0guess(t,t'):")    
    G0gtr= zero;G0less= zero
    do ik=1,2*L
       en=wr(ik)
       nfp=Xc('>',en)
       nfm=Xc('<',en)
       dos=-aimag(eqG0w(ik))/pi
       do i=0,nstep
          do j=0,nstep
             peso=exp(-xi*en*(t(i)-t(j)))
             G0less(i,j)=G0less(i,j)+xi*nfm*peso*dos*fmesh
             G0gtr(i,j)=G0gtr(i,j)  +xi*nfp*peso*dos*fmesh
          enddo
       enddo
    enddo
    if(mpiID==0)then
       call plot_dislin3D("G0less_t1t2","X/$\Delta t$","Y/$\Delta t$","Z", &
            t(0:nstep)/dt,t(0:nstep)/dt,G0less(0:nstep,0:nstep))
       call plot_dislin3D("G0gtr_t1t2","X/$\Delta t$","Y/$\Delta t$","Z", &
            t(0:nstep)/dt,t(0:nstep)/dt,G0gtr(0:nstep,0:nstep))
       call dump("")    
    endif
    return
  end subroutine getG0guess
  !********************************************************************
  !********************************************************************
  !********************************************************************








  !+-------------------------------------------------------------------+
  !PROGRAM  : BUILD_BATH_
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial fields G0^{<,>,R,\lceil} used to obtain
  !the self-energies functions.
  !+-------------------------------------------------------------------+
  subroutine Build_bath_t()
    integer :: i,j,ik,jk,im,itau,jtau,istep,it,m
    real(8) :: w,en,intE
    complex(8) :: peso,nt
    real(8) :: nfp,nfm,arg
    complex(8),dimension(-nstep:nstep) :: Gbless,Gbgtr

    if(mpiID==0)call dump("Get Bath:")    
    Gbless=zero
    Gbgtr=zero
    do im=1,Lmu 
       en=epsimu(im)
       nfm=Xc("<",en)
       nfp=Xc(">",en)
       do i=-nstep,nstep
          arg=t(i)!-t(j)
          arg=en*arg
          peso=exp(-xi*arg)
          Gbless(i)=Gbless(i)+ xi*nfm*peso/dble(Lmu)
          Gbgtr(i) =Gbgtr(i) + xi*nfp*peso/dble(Lmu)
       enddo
    enddo

    !Get the Self-energy contribution of the bath:
    !=================================================
    S0less=Vpd**2*Gbless
    S0gtr =Vpd**2*Gbgtr
    do i=-nstep,nstep
       S0ret(i)=heaviside(t(i))*(S0gtr(i)-S0less(i))
    enddo
    if(mpiID==0)then
       call splot("bathG0less_t.ipt",t(0:nstep),aimag(Gbless(0:nstep)),real(Gbless(0:nstep)))
       call splot("bathG0gtr_t.ipt",t(0:nstep),aimag(Gbgtr(0:nstep)),real(Gbgtr(0:nstep)))
       call splot("S0less_t.ipt",t(0:nstep),aimag(S0less(0:nstep)),real(S0less(0:nstep)))
       call splot("S0gtr_t.ipt",t(0:nstep),aimag(S0gtr(0:nstep)),real(S0gtr(0:nstep)))
       call splot("S0ret_t.ipt",t(0:nstep),aimag(S0ret(0:nstep)),real(S0ret(0:nstep)))
       call dump("done");call dump("")
    endif
    return
  end subroutine Build_bath_t
  !********************************************************************
  !********************************************************************
  !********************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  : BUILD_BATH
  !TYPE     : Subroutine
  !PURPOSE  : Build the initial fields G0^{<,>,R,\lceil} used to obtain
  !the self-energies functions.
  !COMMENT  : Construct the Bath-Sigma functions from freq. domain
  !(that is 'cause you know the analytic expression, same
  !procedure as in IPTkeldysh @equilibrium)
  !This turns out to be the same provided one takes 
  !the limit \eta\goto 0.d0
  !+-------------------------------------------------------------------+
  subroutine Build_bath_w
    integer :: i,j,im
    real(8) :: ex,A,An,w,beta0
    complex(8) :: iw,fg
    complex(8),dimension(2*L) :: gb0fless,gb0fgtr
    complex(8),dimension(-L:L):: gb0tgtr,gb0tless,exa
    if(mpiID==0)call dump("Get Bath:")    
    beta0=beta!/10
    ex=-1.d0       
    do i=-L,L
       ex=-ex
       exa(i)=ex
    enddo
    do i=1,2*L
       w = wr(i)
       iw= cmplx(w,eps)
       fg=zero
       do im=1,Lmu
          fg=fg+1.d0/(iw-epsimu(im))/dble(Lmu)
       enddo
       A = -aimag(fg)/pi
       An= A*fermi(w,beta0)
       gb0fless(i)= pi2*xi*An
       gb0fgtr(i) = pi2*xi*(An-A)
    enddo
    call cfft_rw2rt(gb0fless,gb0tless,L)  ; gb0tless=fmesh/pi2*gb0tless
    call cfft_rw2rt(gb0fgtr, gb0tgtr,L)   ; gb0tgtr =fmesh/pi2*gb0tgtr
    gb0tgtr =exa*gb0tgtr
    gb0tless=exa*gb0tless

    !Get the Self-energy contribution of the bath:
    !=================================================
    S0gtr=zero
    S0less=zero
    S0ret =zero
    S0less=Vpd**2*gb0tless(-nstep:nstep)
    S0gtr =Vpd**2*gb0tgtr(-nstep:nstep)
    do i=-nstep,nstep
       S0ret(i)=heaviside(t(i))*(S0gtr(i)-S0less(i))
    enddo
    if(mpiID==0)then
       call splot("bathG0less_t.ipt",t(0:nstep),aimag(gb0tless(0:nstep)),real(gb0tless(0:nstep)),TT)
       call splot("bathG0gtr_t.ipt",t(0:nstep),aimag(gb0tgtr(0:nstep)),real(gb0tgtr(0:nstep)),TT)
       call splot("S0less_t.ipt",t(0:nstep),aimag(S0less(0:nstep)),real(S0less(0:nstep)),TT)
       call splot("S0gtr_t.ipt",t(0:nstep),aimag(S0gtr(0:nstep)),real(S0gtr(0:nstep)),TT)
       call splot("S0ret_t.ipt",t(0:nstep),aimag(S0ret(0:nstep)),real(S0ret(0:nstep)),TT)
       call dump("done");call dump("")
    endif
    return
  end subroutine Build_bath_w
  !********************************************************************
  !********************************************************************
  !********************************************************************









  !+-------------------------------------------------------------------+
  !PROGRAM  : DYSONG0 
  !TYPE     : Subroutine
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !NOTE     :
  !+-------------------------------------------------------------------+
  subroutine dysonG0()
    integer :: i,j,l,m,itau,jtau,ik
    complex(8) :: An
    integer :: Nt
    complex(8),dimension(nstep+1,nstep+1) :: Uno,GammaRet,GammaAdv
    complex(8),dimension(nstep+1,nstep+1) :: dGless,dGgtr,dSret,dSless,dSgtr
    complex(8),dimension(nstep+1,nstep+1) :: dGret,dG0ret,dG0adv,dG0less,dG0gtr

    complex(8),dimension(-nstep:nstep)   :: g0tless,g0tgtr,g0tret
    complex(8),dimension(-nstep:nstep)   :: gtless,gtgtr,gtret


    if(mpiID==0)call dump("Update WF: Dyson")
    Nt=nstep+1

    call shiftFW(locGless,dGless)
    call shiftFW(locGgtr,dGgtr)
    call shiftFW(Sret,dSret)
    call shiftFW(Sless,dSless)
    call shiftFW(Sgtr,dSgtr)
    dGret=zero
    !Get G_loc^{A,R}, S^{A,R}
    do i=1,Nt
       do j=1,Nt
          dGret(i,j)=heaviside(t(i)-t(j))*(dGgtr(i,j) - dGless(i,j))
       enddo
    enddo
    forall(i=1:Nt,j=1:Nt)gtret(i-j)=dGret(i,j)
    if(mpiID==0)then
       call splot("dysonGret_t",t(0:nstep),aimag(gtret(0:nstep)),real(gtret(0:nstep)))
       call plot_dislin3D("dysonGret_t1t2","X/$\Delta t$","Y/$\Delta t$","Z",&
            t(0:nstep)/dt,t(0:nstep)/dt,dGret) 
    endif

    !Update G_0^R using Dyson: G_0^R= [G^R^-1 - Sigma^R]^-1 = \Gamma^R^-1\cdot G^R    
    !\Gamma_R^{-1} = [1+G_loc^R*\Sigma^R]^{-1} : \Gamma_A=[\Gamma_R]^+
    Uno=zero ; forall(i=1:Nt)Uno(i,i)=one/dt
    GammaRet = Uno+matmul(dGret,dSret)*dt  
    call InvMat(dt*GammaRet*dt,Nt) 
    dG0ret = matmul(GammaRet,dGret)*dt
    forall(i=1:Nt,j=1:Nt)g0tret(i-j)=dG0ret(i,j)
    if(mpiID==0)then
       call splot("dysonG0ret_t",t(0:nstep),aimag(g0tret(0:nstep)),real(g0tret(0:nstep)))
       call plot_dislin3D("dysonG0ret_t1t2","X","Y","Z",t(0:nstep)/dt,t(0:nstep)/dt,dG0ret)
    endif

    !Update G_0^<, G0^>
    dG0less=matmul(GammaRet,matmul(dGless,conjg(transpose(GammaRet))))*dt**2-&
         matmul(dG0ret,matmul(dSless,conjg(transpose(dG0ret))))*dt**2
    dG0gtr =matmul(GammaRet,matmul(dGgtr,conjg(transpose(GammaRet))))*dt**2 -&
         matmul(dG0ret,matmul(dSgtr,conjg(transpose(dG0ret))))*dt**2
    forall(i=1:Nt,j=1:Nt)
       g0tless(i-j)=dG0less(i,j)
       g0tgtr(i-j) =dG0gtr(i,j)
    end forall
    if(mpiID==0)then
       call splot("dysonG0less_t",t(0:nstep),aimag(g0tless(0:nstep)),real(g0tless(0:nstep)))
       call splot("dysonG0gtr_t",t(0:nstep),aimag(g0tgtr(0:nstep)),real(g0tgtr(0:nstep)))
       call plot_dislin3D("dysonG0less_t1t2","X","Y","Z",t(0:nstep)/dt,t(0:nstep)/dt,dG0less) 
       call plot_dislin3D("dysonG0gtr_t1t2","X","Y","Z",t(0:nstep)/dt,t(0:nstep)/dt,dG0gtr)
    endif
    call shiftBW(dG0less,G0less)
    call shiftBW(dG0gtr,G0gtr)
    do i=1,3
       if(mpiID==0)call dump("")
    enddo
  end subroutine dysonG0
  !********************************************************************
  !********************************************************************
  !********************************************************************








  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine update_sigma
    integer    :: i,j,itau,jtau,ndim,ik,jk
    complex(8) :: stret(-nstep:nstep),sfret(2*nstep)
    if(mpiID==0)call dump("Get Sigma(t,t')")
    !Get Sigma^>,<
    Sless =zero
    Sgtr  =zero
    Sret  =zero
    do i=0,nstep
       do j=0,nstep
          Sless(i,j)=U**2*(G0less(i,j)**2)*G0gtr(j,i)
          Sgtr (i,j)=U**2*(G0gtr(i,j)**2)*G0less(j,i)
          Sret(i,j) =heaviside(t(i)-t(j))*(Sgtr(i,j)-Sless(i,j))          
          stret(i-j)=Sret(i,j)
       enddo
    enddo
    if(mpiID==0)then
       call plot_dislin3D("Sless_t1t2","X/$\Delta t$","Y/$\Delta t$","Z",&
            t(0:nstep)/dt,t(0:nstep)/dt,Sless(0:nstep,0:nstep))
       call plot_dislin3D("Sgtr_t1t2","X/$\Delta t$","Y/$\Delta t$","Z",&
            t(0:nstep)/dt,t(0:nstep)/dt,Sgtr(0:nstep,0:nstep))
       call cfft_rt2rw(stret,sfret,nstep) ;    sfret=sfret*dt ; call swap_fftrt2rw(sfret)
       call splot("Sret_realw",wrmini,aimag(sfret),real(sfret))
       call dump("")
    endif
    return
  end subroutine update_sigma
  !********************************************************************
  !********************************************************************
  !********************************************************************



  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine shiftFW(Gin,Gout)
    complex(8),dimension(0:,0:) :: Gin
    complex(8),dimension(:,:)   :: Gout
    integer :: i,j,Ndim1,Ndim2
    Ndim1=size(Gout,1)
    Ndim2=size(Gout,2)
    do i=1,Ndim1
       do j=1,Ndim2
          Gout(i,j)=Gin(i-1,j-1)
       enddo
    enddo
  end subroutine shiftFW



  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine shiftBW(Gin,Gout)
    complex(8),dimension(0:,0:) :: Gout
    complex(8),dimension(:,:)   :: Gin
    integer :: i,j,Ndim1,Ndim2
    Ndim1=size(Gin,1)
    Ndim2=size(Gin,2)
    do i=1,Ndim1
       do j=1,Ndim2
          Gout(i-1,j-1)=Gin(i,j)
       enddo
    enddo
  end subroutine shiftBW

end module FUNX_NEQ
