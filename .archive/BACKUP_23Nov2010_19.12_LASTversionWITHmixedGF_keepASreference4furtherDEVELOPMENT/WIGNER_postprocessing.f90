program WIGNER_ROTATION
  !USE COMMON_VARS
  USE VARS_GLOBAL
  USE GRIDS
  USE DLPLOT
  USE TOOLS
  USE FFTW
  implicit none
  integer           :: i,j,ik,k,m,ia,ir,irel,iave
  complex(8)        :: delta

  !Wigner variables:
  real(8),dimension(:),allocatable       :: trel,tave
  complex(8),dimension(:,:),allocatable  :: locGret,Sret
  complex(8),dimension(:,:),allocatable  :: wgnGgtr,wgnGless,wgnGret
  complex(8),dimension(:,:),allocatable  :: wgnSgtr,wgnSless,wgnSret
  complex(8),dimension(:,:),allocatable  :: gfret_wgn,sfret_wgn,gfless_wgn
  real(8),dimension(:,:),allocatable     :: nf_wgn

  call read_input_init("inputFILE.in")
  !+ some wigner specific variables:
  allocate(trel(-nstep:nstep),tave(0:nstep))
  allocate(locGless(0:nstep,0:nstep),locGgtr(0:nstep,0:nstep),locGret(0:nstep,0:nstep))
  allocate(Sless(0:nstep,0:nstep),Sgtr(0:nstep,0:nstep),Sret(0:nstep,0:nstep))
  allocate(wgnGless(-nstep:nstep,0:nstep),&
       wgnGgtr(-nstep:nstep,0:nstep),     &
       wgnGret(-nstep:nstep,0:nstep))
  allocate(wgnSless(-nstep:nstep,0:nstep),&
       wgnSgtr(-nstep:nstep,0:nstep),     &
       wgnSret(-nstep:nstep,0:nstep))
  allocate(gfret_wgn(0:nstep,2*nstep),gfless_wgn(0:nstep,2*nstep),sfret_wgn(0:nstep,2*nstep))
  allocate(nf_wgn(0:nstep,2*nstep))


  call init_trel(t,trel,nstep)
  call init_tave(t,tave,nstep)

  !Read the functions:
  call read_data(locGless(0:nstep,0:nstep),"DATAfunx/locGless_data",nstep)
  call read_data(locGgtr(0:nstep,0:nstep),"DATAfunx/locGgtr_data",nstep)
  call read_data(Sless(0:nstep,0:nstep),"DATAfunx/Sless_data",nstep)
  call read_data(Sgtr(0:nstep,0:nstep),"DATAfunx/Sgtr_data",nstep)
  forall(i=0:nstep,j=0:nstep)locGret(i,j)=heaviside(t(i)-t(j))*(locGgtr(i,j)-locGless(i,j))
  forall(i=0:nstep,j=0:nstep)Sret(i,j)=heaviside(t(i)-t(j))*(Sgtr(i,j)-Sless(i,j))

  !Perform the Wigner Rotation:
  wgnGless= wigner_transform(locGless,nstep)
  wgnGret = wigner_transform(locGret,nstep)
  wgnSless= wigner_transform(Sless,nstep)
  wgnSret = wigner_transform(Sret,nstep)
  call system("if [ ! -d WIGNER ]; then mkdir WIGNER; fi")
  call plot_3D("wgnGless3D","X","Y","Z",trel(-nstep:nstep)/dt,tave(1:nstep)/dt,wgnGless(-nstep:nstep,1:nstep))
  call plot_3D("wgnSless3D","X","Y","Z",trel(-nstep:nstep)/dt,tave(1:nstep)/dt,wgnSless(-nstep:nstep,1:nstep))
  call system("mv wgn*3D WIGNER/")

  delta=10.0*(one+xi)/dble(nstep)
  do ia=0,nstep
     call cfft_rt2rw(wgnGret(:,ia),gfret_wgn(ia,:),nstep)  ;gfret_wgn(ia,:)=gfret_wgn(ia,:)*dt;call swap_fftrt2rw(gfret_wgn(ia,:))
     call cfft_rt2rw(wgnGless(:,ia),gfless_wgn(ia,:),nstep);gfless_wgn(ia,:)=gfless_wgn(ia,:)*dt;call swap_fftrt2rw(gfless_wgn(ia,:))
     call cfft_rt2rw(wgnSret(:,ia),sfret_wgn(ia,:),nstep)  ;sfret_wgn(ia,:)=sfret_wgn(ia,:)*dt;call swap_fftrt2rw(sfret_wgn(ia,:))
     call splot("WIGNER/wgnDOS.ipt",wr,(-aimag(gfret_wgn(ia,:) - delta*dble(ia)*pi)/pi),append=TT)
     call splot("WIGNER/wgnSigma_realw.ipt",wr,(sfret_wgn(ia,:) + delta*dble(ia)),append=TT)
     call splot("WIGNER/wgnGless_realw.ipt",wr,(gfless_wgn(ia,:) + delta*dble(ia)),append=TT)
     nf_wgn(ia,:) = -aimag(gfless_wgn(ia,:))!/aimag(gfret_wgn(ia,:))/pi2
  enddo
  call splot("wgndosVSrealwVStime.ipt",tave(0:nstep),wr(1:2*nstep),-aimag(gfret_wgn(0:nstep,1:2*nstep))/pi)
  call splot("wgnnfVSrealwVStime.ipt",tave(0:nstep),wr(1:2*nstep),nf_wgn(0:nstep,1:2*nstep))
  call system("mv wgn* plot*wgn* WIGNER/")
  !call plot_3D("wgnDOS3D","X","Y","Z",tave(0:nstep)/dt,wr(1:2*nstep),-aimag(gfret_wgn(0:nstep,1:2*nstep))/pi)



contains
  !-------------------------------------
  subroutine read_data(Gin,inFILE,M)
    character(len=*)              :: inFILE
    integer                       :: M,i,j
    complex(8),dimension(0:M,0:M) :: Gin
    open(20,file=trim(inFILE))
    do i=0,M
       do j=0,M
          read(20,*)Gin(i,j)
       enddo
    enddo
    close(20)
  end subroutine read_data
  !-------------------------------------

  !-------------------------------------
  function wigner_transform(Gin,M)
    integer                       :: M,ir,ia,i,j
    complex(8),dimension(0:M,0:M) :: Gin
    complex(8),dimension(-M:M,0:M):: Gout,wigner_transform
    do ir=-nstep,nstep
       do ia=0,nstep
          forall(j=0:nstep,i=0:nstep, i-j==ir .AND. (i+j)/2==ia)
             Gout(ir,ia)=Gin(i,j)
          end forall
       enddo
    enddo
    wigner_transform=Gout
  end function wigner_transform
  !-------------------------------------

  !-------------------------------------
  subroutine init_trel(time,tr,M)
    integer                            :: i,j,ir
    integer,intent(in)                 :: M
    real(8),dimension(-M:M),intent(in) :: time
    real(8),dimension(-M:M),intent(out):: tr
    do ir=-nstep,nstep
       forall(j=0:nstep,i=0:nstep, i-j==ir)tr(ir)=time(i)-time(j)
    enddo
    return
  end subroutine init_trel
  !-------------------------------------

  !-------------------------------------
  subroutine init_tave(time,ta,M)
    integer                            :: i,j,ia
    integer,intent(in)                 :: M
    real(8),dimension(-M:M),intent(in) :: time
    real(8),dimension(0:M),intent(out) :: ta
    do ia=0,nstep
       forall(j=0:nstep,i=0:nstep, (i+j)/2==ia)ta(ia)=(time(i)+time(j))/2.d0
    enddo
    return
  end subroutine init_tave

end program WIGNER_ROTATION
