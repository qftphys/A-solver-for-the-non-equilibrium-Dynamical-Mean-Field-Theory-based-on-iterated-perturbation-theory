PROGRAM neq_keldysh
  !########################################################################
  !     PROGRAM  : neq_keldysh
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model with applied external electric
  !                field E, on the hypercubic lattice, using neq. DMFT 
  !     AUTHORS  : Adriano Amaricci
  !     LAST UPDATE: 07/2009
  !########################################################################
  USE VARIABILI !definisce le variabili globali sharate dove serve
  USE FUNZIONI  !funzioni utili in altri parti del programma
  USE GLOC     !operazioni sulle matricione \hat G_{(0)loc}(t,t') 

  implicit none
  include "mpif.h" 
  integer :: i,j
  real(8) :: thetap,thetam

  namelist/variables/beta,U,xmu,Efield
  namelist/parameters/L,Lepsi,emax,emin,fmesh,nloop

  call MPI_INIT(mpiIERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiMYID,mpiIERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiIERR)
  call MPI_BARRIER(MPI_COMM_WORLD,mpiIERR)
  print*,'Processor ',mpiMYID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiIERR)

  open(10,file='inputFILE.in')
  read(10,nml=variables)
  read(10,nml=parameters)
  close(10)

  ts=1.d0/sqrt2 !t_hopping
  D=one
  temp=1.d0/beta

  N=2*L+1
  Nepsi=2*Lepsi+1
  allocate(t(N),e(Nepsi))

  !time array
  dt=pi2/(fmesh*dble(N))           
  tmin=-dble(L)*dt
  do i=1,N
     t(i)=tmin+dble(i-1)*dt
  enddo

  !energy-array
  de=(emax-emin)/dble(Nepsi)
  do i=1,Nepsi
     e(i)=emin+dble(i-1)*de
  enddo

  allocate(g0loc(4,N,N))
  allocate(g0hat(2*N,2*N),ghat(2*N,2*N),sighat(2*N,2*N),fgloc(2*N,2*N))
!  allocate(g0epsiloc(2*N,2*N))!,glocal(2*N,2*N),glocal2(2*N,2*N),identity(2*N,2*N))


  !Starts the DMFT iteration
  do iloop=1,nloop
     if(mpiMYID==0)then
        print*, 'iloop',iloop,'/',nloop
        if(iloop.eq.1)then
           print*,'Entro nel Guess'
           !costruisci i blocchi della matrice \hat{G0_loc}, G0loc_(c=t,<,>,tbar)
           do i=1,size(vchar)
              char=vchar(i)
              print*,''
              print*,char
              call build_g0loc(char,g0loc(i,:,:)) !ogni blocco g0loc_c
           enddo
           print*,'Gather blocks'
           call gather_blocks(g0loc,g0hat)       !mette insieme in g0hat
        else
           print*,''
           !implementa la self-consistenza \hat G0_loc = [(\hat G)^-1 + Sigma]^-1
           !inverti \hat G che deve essere definita alla fine del dmft-loop
           call InvMat(ghat,2*N) !Ottiene (\hat G)^-1
           ghat=ghat+sighat            !fa la somma
           call InvMat(ghat,2*N) !la re-inverte
           g0hat=ghat                  !Dyson equation
           !estrarre g0loc_c=>,< dall matriciona g0hat
           do i=1,size(vchar)
              char=vchar(i)
              call extract_g0loc(char,g0hat,g0loc(i,:,:))
           enddo
        endif
        print*,'Sono arrivato all Sigma(t,t`)'


        !Get Sigma(t,t')
        do i=1,N
           do j=1,N !Recall that the 2nd block is defined with minus sign
              g0loc(2,i,j)=-(U**2)*(g0loc(2,i,j)**2)*(g0loc(2,j,i))
              g0loc(3,i,j)=(U**2)*(g0loc(3,i,j)**2)*(g0loc(3,j,i))
              thetap=heaviside(t(i)-t(j))
              thetam=heaviside(t(j)-t(i))
              g0loc(1,i,j)=-thetap*g0loc(2,i,j)+thetam*g0loc(3,i,j)
              g0loc(4,i,j)=thetap*g0loc(3,i,j)-thetam*g0loc(2,i,j)
           enddo
        enddo
        call gather_blocks(g0loc,sighat)
     endif
     call MPI_BCAST(sighat,2*N*2*N,MPI_COMPLEX,0,MPI_COMM_WORLD,mpiIERR)
     print*,''
     !Get G_loc(t,t')=ghat
     !All the slaves do have a Sighat copy.
     print*,'sto per entrare in get_loc'
     call get_gloc(fgloc,N)
     call MPI_REDUCE(fgloc,ghat,2*N*2*N,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,mpiIERR)
  enddo
end PROGRAM neq_keldysh
