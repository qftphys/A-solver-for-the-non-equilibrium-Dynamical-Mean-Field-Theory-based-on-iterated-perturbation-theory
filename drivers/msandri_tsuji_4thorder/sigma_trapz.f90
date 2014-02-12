  subroutine calc_sigma_4th(G, Sigma4,params)
    type(kb_contour_gf)                   :: G
    type(kb_contour_gf)                   :: Sigma4(4)
    type(kb_contour_params)               :: params
    integer                               :: N,L
    complex(8),dimension(:,:),allocatable :: G_gtr,G_31,G_33,Sigma_21,chi_12,chi_21
    integer                               :: it,itp,it1,it2,itau,itaup,itau1,itau2
    integer                               :: Ntot
    real(8),dimension(:),allocatable      :: Ut
    complex(8)                            :: int2,int1
    real(8)                               :: dt,dtau
    complex(8),dimension(:),allocatable   :: dum1,dum2



    N = params%Nt
    L = params%Ntau
    dt= params%dt
    dtau= params%dtau
    
    Ntot = 2*N+L+1
    allocate(Ut(Ntot))
    allocate(dum1(0:max(L,N)))
    allocate(dum2(0:max(L,N)))

    Ut(1:N)           = U 
    Ut(N+1:2*N)       = U
    Ut(2*N+1:2*N+L+1) = Ui

    allocate(G_gtr(N,N),G_31(0:L,N),Sigma_21(N,N),G_33(0:L,0:L),chi_12(N,N),chi_21(N,N))

    do i=1,N
       do j=1,i
          G_gtr(i,j) = G%less(i,j)+G%ret(i,j)
       enddo
       do j=i+1,N
          G_gtr(i,j) = G%less(i,j)-conjg(G%ret(j,i))
       enddo
    enddo

    do i=0,L
       do j=1,N
          G_31(i,j) = conjg(G%lmix(j,L-i))
       enddo
    enddo

    do i=0,L
       do j=0,i
          G_33(i,j) = xi*G%mats(i-j)
       enddo
       do j=i+1,L
          G_33(i,j) = -xi*G%mats(L+i-j)
       enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SIGMA^4a!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'Get 4a'
    chi_12=zero
    chi_21=zero
    Sigma4(1) = zero
    Sigma_21 = zero

    do it=1,N
       do itp=1,N
          dum1=zero
          do it1=1,it
             dum1(it1)=Ut(it1)*(G_gtr(it,it1)**2-G%less(it,it1)**2)*G%less(it1,itp)**2
          enddo
          chi_12(it,itp) = chi_12(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          dum1=zero
          do it1=1,itp
             dum1(it1)=Ut(it1)*G%less(it,it1)**2*(G%less(it1,itp)**2-G_gtr(it1,itp)**2)
          enddo
          chi_12(it,itp) = chi_12(it,itp) + dt*kb_trapz(dum1(0:),1,itp)
       enddo
    enddo

    do it=1,N
       do itp=1,N
          do it1=1,it
             dum1(it1) = Ut(it1)*(G_gtr(it,it1)**2-G%less(it,it1)**2)*G_gtr(it1,itp)**2
          enddo
          chi_21(it,itp) = chi_21(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          do it1=1,itp
             dum1(it1) = Ut(it1)*G_gtr(it,it1)**2*(G%less(it1,itp)**2-G_gtr(it1,itp)**2)
          enddo
          chi_21(it,itp) = chi_21(it,itp) + dt*kb_trapz(dum1(0:),1,itp)
       enddo
    enddo
    
    do it=1,N
       do itp=1,N
          dum1=zero
          do it1=1,it
             dum1(it1) = Ut(it1)*(G_gtr(it,it1)**2-G%less(it,it1)**2)*chi_12(it1,itp)
          enddo
          Sigma4(1)%less(it,itp) = Sigma4(1)%less(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          dum1=zero
          do it1=1,itp
             dum1(it1) = Ut(it1)*G%less(it,it1)**2*(chi_12(it1,itp)-chi_21(it1,itp))
          enddo
          Sigma4(1)%less(it,itp) = Sigma4(1)%less(it,itp) + dt*kb_trapz(dum1(0:),1,itp)

          Sigma4(1)%less(it,itp) = Sigma4(1)%less(it,itp)*&
               Ut(it)*Ut(itp)*G%less(it,itp)
       enddo
    enddo

    do it=1,N
       do itp=1,N
          dum1=zero
          do it1=1,it
             dum1(it1) = Ut(it1)*(G_gtr(it,it1)**2-G%less(it,it1)**2)*chi_21(it1,itp)
          enddo
          Sigma_21(it,itp) = Sigma_21(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          dum1=zero
          do it1=1,itp
             dum1(it1) = Ut(it1)*G_gtr(it,it1)**2*(chi_12(it1,itp)-chi_21(it1,itp))
          enddo
          Sigma_21(it,itp) = Sigma_21(it,itp) + dt*kb_trapz(dum1(0:),1,itp)

          Sigma_21(it,itp) = Sigma_21(it,itp)&
               *Ut(it)*Ut(itp)*G_gtr(it,itp)
       enddo
    enddo


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Get retarded component
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N
       Sigma4(1)%ret(N,i) = Sigma_21(N,i) - Sigma4(1)%less(N,i)
    enddo
    



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SIGMA^4b!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'Get 4b'
    chi_12=zero
    chi_21=zero
    Sigma4(2) = zero
    Sigma_21 = zero

    do it=1,N
       do itp=1,N
          dum1=zero
          do it1=1,it
             dum1(it1)=Ut(it1)*(G_gtr(it,it1)**3-G%less(it,it1)**3)*G%less(it1,itp)
          enddo
          chi_12(it,itp) = chi_12(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          dum1=zero
          do it1=1,itp
             dum1(it1)=Ut(it1)*G%less(it,it1)**3*(G%less(it1,itp)-G_gtr(it1,itp))
          enddo
          chi_12(it,itp) = chi_12(it,itp) + dt*kb_trapz(dum1(0:),1,itp)
       enddo
    enddo

    do it=1,N
       do itp=1,N
          dum1=zero
          do it1=1,it
             dum1(it1)=Ut(it1)*(G_gtr(it,it1)**3-G%less(it,it1)**3)*G_gtr(it1,itp)
          enddo
          chi_21(it,itp) = chi_21(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          dum1=zero
          do it1=1,itp
             dum1(it1)=Ut(it1)*G_gtr(it,it1)**3*(G%less(it1,itp)-G_gtr(it1,itp))
          enddo
          chi_21(it,itp) = chi_21(it,itp) + dt*kb_trapz(dum1(0:),1,itp)
       enddo
    enddo
    
    do it=1,N
       do itp=1,N
          dum1=zero
          do it1=1,it
             dum1(it1)=Ut(it1)*(G_gtr(it,it1)-G%less(it,it1))*chi_12(it1,itp)
          enddo
          Sigma4(2)%less(it,itp) = Sigma4(2)%less(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          dum1=zero
          do it1=1,itp
             dum1(it1)=Ut(it1)*G%less(it,it1)*(chi_12(it1,itp)-chi_21(it1,itp))
          enddo
          Sigma4(2)%less(it,itp) = Sigma4(2)%less(it,itp) + dt*kb_trapz(dum1(0:),1,itp)

          Sigma4(2)%less(it,itp) = Sigma4(2)%less(it,itp)*&
               Ut(it)*Ut(itp)*G%less(it,itp)**2
       enddo
    enddo

    do it=1,N
       do itp=1,N
          dum1=zero
          do it1=1,it
             dum1(it1)=Ut(it1)*(G_gtr(it,it1)-G%less(it,it1))*chi_21(it1,itp)
          enddo
          Sigma_21(it,itp) = Sigma_21(it,itp) + dt*kb_trapz(dum1(0:),1,it)
          dum1=zero
          do it1=1,itp
             dum1(it1)=Ut(it1)*G_gtr(it,it1)*(chi_12(it1,itp)-chi_21(it1,itp))
          enddo
          Sigma_21(it,itp) = Sigma_21(it,itp) + dt*kb_trapz(dum1(0:),1,itp)

          Sigma_21(it,itp) = Sigma_21(it,itp)&
               *Ut(it)*Ut(itp)*G_gtr(it,itp)**2
       enddo
    enddo


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Get retarded component
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N
       Sigma4(2)%ret(N,i) = Sigma_21(N,i) - Sigma4(2)%less(N,i)
    enddo






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SIGMA^4c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'Get 4c'
    sigma4(3)=zero
    sigma_21=zero
    !!LESS COMPONENT
    do it=1,N
       !!$OMP PARALLEL DO PRIVATE(int1,int2,itp,it1,it2,itau1,itau2)
       do itp=1,N
          int1=zero
          
          do it1=1,it
             int2=zero
             do it2=1,it1
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it1)
             do it2=it1,it
                dum2(it2) = Ut(it2)*G%less(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1+1,it)
             do it2=it,N
                dum2(it2) = Ut(it2)*G%less(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it+1,N)
             do it2=N,itp,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,N-1)
             do it2=itp,1,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,itp-1)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = Ut(it1)*G_gtr(it,it1)*G%less(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),1,it)


          do it1=it,N
             int2=zero
             do it2=1,it
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it)
             do it2=it,it1
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it,it1)
             do it2=it1,N
                dum2(it2) = Ut(it2)*G%less(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1,N)
             do it2=N,itp,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,N)
             do it2=itp,1,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,itp)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = Ut(it1)*G%less(it,it1)*G%less(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),it,N)


          do it1=N,itp,-1
             int2=zero
             do it2=1,it
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it)
             do it2=it,N
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it,N)
             do it2=N,it1,-1
                dum2(it2) = -Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1,N)
             do it2=it1,itp,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,it1)
             do it2=itp,1,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,itp)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = -Ut(it1)*G%less(it,it1)*G%less(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),itp,N)


          do it1=itp,1,-1
             int2=zero
             do it2=1,it
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it)
             do it2=it,N
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it,N)
             do it2=N,itp,-1
                dum2(it2) = -Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,N)
             do it2=itp,it1,-1
                dum2(it2) = -Ut(it2)*G_gtr(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1,itp)
             do it2=it1,1,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it1)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = -Ut(it1)*G%less(it,it1)*G_gtr(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),1,itp)

          Sigma4(3)%less(it,itp) = Ut(it)*Ut(itp)*int1

       enddo
       !!$OMP END PARALLEL DO
    enddo
   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!GREATER COMPONENT
    do it=1,N
       !!$OMP PARALLEL DO PRIVATE(int1,int2,itp,it1,it2,itau1,itau2)
       do itp=1,N
          int1=zero
          do it1=1,itp

             int2=zero
             do it2=1,it1
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it1)
             do it2=it1,itp
                dum2(it2) = Ut(it2)*G%less(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1,itp)
             do it2=itp,N
                dum2(it2) = Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,N)
             do it2=N,it,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it,N)
             do it2=it,1,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = Ut(it1)*G_gtr(it,it1)*G%less(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),1,itp)


          do it1=itp,N
             int2=zero
             do it2=1,itp
                dum2(it2) =  Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,itp)
             do it2=itp,it1
                dum2(it2) =  Ut(it2)*G_gtr(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,it1)
             do it2=it1,N
                dum2(it2) =  Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1,N)
             do it2=N,it,-1
                dum2(it2) =  -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it,N)
             do it2=it,1,-1
                dum2(it2) =  -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it)
             do itau2=0,L
                dum2(itau2) =  -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = Ut(it1)*G_gtr(it,it1)*G_gtr(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),itp,N)


          do it1=N,it,-1
             int2=zero
             do it2=1,itp
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,itp)
             do it2=itp,N
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,N)
             do it2=N,it1,-1
                dum2(it2) = -Ut(it2)*G_gtr(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1,N)
             do it2=it1,it,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it,it1)
             do it2=it,1,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = -Ut(it1)*G_gtr(it,it1)*G_gtr(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),it,N)


          do it1=it,1,-1
             int2=zero
             do it2=1,itp
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G%less(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,itp)
             do it2=itp,N
                dum2(it2) = Ut(it2)*G_gtr(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),itp,N)
             do it2=N,it,-1
                dum2(it2) = -Ut(it2)*G_gtr(it1,it2)*G_gtr(it2,itp)*G_gtr(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it,N)
             do it2=it,it1,-1
                dum2(it2) = -Ut(it2)*G_gtr(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),it1,it)
             do it2=it1,1,-1
                dum2(it2) = -Ut(it2)*G%less(it1,it2)*G_gtr(it2,itp)*G%less(it,it2)**2
             enddo
             int2 = int2 + dt*kb_trapz(dum2(0:),1,it1)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it1,itau2)*G_31(itau2,itp)*G%lmix(it,itau2)**2
             enddo
             int2 = int2 + dtau*kb_trapz(dum2(0:),0,L)
             dum1(it1) = -Ut(it1)*G%less(it,it1)*G_gtr(it1,itp)**2*int2
          enddo
          int1 = int1 + dt*kb_trapz(dum1(0:),1,it)

          Sigma_21(it,itp) = Ut(it)*Ut(itp)*int1

       enddo
       !!$OMP END PARALLEL DO
    enddo
  
 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Get retarded component
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N
       Sigma4(3)%ret(N,i) = Sigma_21(N,i) - Sigma4(3)%less(N,i)
    enddo
    




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SIGMA^4d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'Get 4d'
    sigma4(4)=zero
    sigma_21=zero
    !!LESS COMPONENT
    do it=1,N
       !!$OMP PARALLEL DO PRIVATE(int1,int2,itp,it1,it2,itau1,itau2)
       do itp=1,N
          int1=zero
          do it1=1,it

             int2=zero
             do it2=1,it1
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it1)
             do it2=it1,it
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,it)
             do it2=it,N
                dum2(it2) = Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,N)
             do it2=N,itp,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,N)
             do it2=itp,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,itp)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,itau2)*G_31(itau2,itp)*G%lmix(it1,itau2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) = Ut(it1)*G_gtr(it,it1)*G%less(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),1,it)

          do it1=it,N
             int2=zero
             do it2=1,it
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it)
             do it2=it,it1
                dum2(it2) = Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,it1)
             do it2=it1,N
                dum2(it2) = Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,N)
             do it2=N,itp,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,N)
             do it2=itp,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,itp)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,itau2)*G_31(itau2,itp)*G%lmix(it1,itau2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) = Ut(it1)*G%less(it,it1)*G%less(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),it,N)


          do it1=N,itp,-1
             int2=zero
             do it2=1,it
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it)
             do it2=it,N
                dum2(it2) = Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,N)
             do it2=N,it1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,N)
             do it2=it1,itp,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,it1)
             do it2=itp,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,itp)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,itau2)*G_31(itau2,itp)*G%lmix(it1,itau2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) = -Ut(it1)*G%less(it,it1)*G%less(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),itp,N)


          do it1=itp,1,-1
             int2=zero
             do it2=1,it
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it)
             do it2=it,N
                dum2(it2) = Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,N)
             do it2=N,itp,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,N)
             do it2=itp,it1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,itp)
             do it2=it1,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it1)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,itau2)*G_31(itau2,itp)*G%lmix(it1,itau2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) = -Ut(it1)*G%less(it,it1)*G_gtr(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),1,itp)
  
          Sigma4(4)%less(it,itp) = Ut(it)*Ut(itp)*G%less(it,itp)*int1

       enddo
       !!$OMP END PARALLEL DO
    enddo



    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!GREATER COMPONENT
    do it=1,N
       !!$OMP PARALLEL DO PRIVATE(int1,int2,itp,it1,it2,itau1,itau2)
       do itp=1,N
          int1=zero
          do it1=1,itp

             int2=zero
             do it2=1,it1
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it1)
             do it2=it1,itp
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,itp)
             do it2=itp,N
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,N)
             do it2=N,it,-1
                dum2(it2) = -Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,N)
             do it2=it,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,it2)*G_31(it2,itp)*G%lmix(it1,it2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) = +Ut(it1)*G_gtr(it,it1)*G%less(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),1,itp)


          do it1=itp,N
             int2=zero
             do it2=1,itp
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,itp)
             do it2=itp,it1
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,it1)
             do it2=it1,N
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,N)
             do it2=N,it,-1
                dum2(it2) = -Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,N)
             do it2=it,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,it2)*G_31(it2,itp)*G%lmix(it1,it2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) = Ut(it1)*G_gtr(it,it1)*G_gtr(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),itp,N)


          do it1=N,it,-1
             int2=zero
             do it2=1,itp
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,itp)
             do it2=itp,N
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,N)
             do it2=N,it1,-1
                dum2(it2) = -Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,N)
             do it2=it1,it,-1
                dum2(it2) = -Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,it1)
             do it2=it,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,it2)*G_31(it2,itp)*G%lmix(it1,it2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) =  -Ut(it1)*G_gtr(it,it1)*G_gtr(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),it,N)


          do it1=it,1,-1
             int2=zero
             do it2=1,itp
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G%less(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,itp)
             do it2=itp,N
                dum2(it2) = Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),itp,N)
             do it2=N,it,-1
                dum2(it2) = -Ut(it2)*G_gtr(it,it2)*G_gtr(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it,N)
             do it2=it,it1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G_gtr(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),it1,it)
             do it2=it1,1,-1
                dum2(it2) = -Ut(it2)*G%less(it,it2)*G_gtr(it2,itp)*G%less(it1,it2)**2*dt
             enddo
             int2 = int2 + kb_trapz(dum2(0:),1,it1)
             do itau2=0,L
                dum2(itau2) = -xi*Ui*G%lmix(it,it2)*G_31(it2,itp)*G%lmix(it1,it2)**2*dtau
             enddo
             int2 = int2 + kb_trapz(dum2(0:),0,L)
             dum1(it1) = -Ut(it1)*G%less(it,it1)*G_gtr(it1,itp)*int2*dt
          enddo
          int1 = int1 + kb_trapz(dum1(0:),1,it)

    
          Sigma_21(it,itp) = Ut(it)*Ut(itp)*G_gtr(it,itp)*int1

       enddo
       !!$OMP END PARALLEL DO
    enddo


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Get retarded component
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N
       Sigma4(4)%ret(N,i) = Sigma_21(N,i) - Sigma4(4)%less(N,i)
    enddo
    




    deallocate(G_gtr,G_31,Sigma_21,G_33)

  end subroutine calc_sigma_4th

