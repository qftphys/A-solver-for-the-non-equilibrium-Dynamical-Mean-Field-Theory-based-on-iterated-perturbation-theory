  !=======Inversion of the Keldysh-Schwinger Matrix==========================
  if(FF)then
     !1) build time/antitime-ordered GF: G^t && G^at; S^t && S^at
     forall(i=0:nstep,j=0:nstep)
        locGtt(i,j) = heaviside(t(i-j))*locGgtr(i,j) + heaviside(t(j-i))*locGless(i,j)
        locGat(i,j) = heaviside(t(j-i))*locGgtr(i,j) + heaviside(t(i-j))*locGless(i,j)
        Stt(i,j) = heaviside(t(i-j))*Sgtr(i,j) + heaviside(t(j-i))*Sless(i,j)
        Sat(i,j) = heaviside(t(j-i))*Sgtr(i,j) + heaviside(t(i-j))*Sless(i,j)
     end forall
     ! forall(i=0:nstep)
     !    locGtt(i,i) = locGless(i,i) 
     !    locGat(i,i) = locGgtr(i,i) 
     !    Stt(i,i) = Sless(i,i) 
     !    Sat(i,i) = Sgtr(i,i) 
     ! end forall
     !call plot_3D("dGtt3D","X","Y","Z",t(0:nstep),t(0:nstep),locGtt(0:nstep,0:nstep))
     !call plot_3D("dStt3D","X","Y","Z",t(0:nstep),t(0:nstep),Stt(0:nstep,0:nstep))

     !2) Build the KS matrices: \FF = {{FF^t , FF^>}, {-FF^<, -FF^at}}
     NN=nstep+1                  !Size of the KS matrices
     allocate(locGmat(1:2*NN,1:2*NN),Smat(1:2*NN,1:2*NN))
     allocate(G0mat(1:2*NN,1:2*NN),GammaMat(1:2*NN,1:2*NN),UnoMat(1:2*NN,1:2*NN))
     forall(i=1:NN,j=1:NN)
        locGmat(i,j)       = locGtt(i-1,j-1)   !11
        locGmat(i,NN+j)    = locGgtr(i-1,j-1)  !12
        locGmat(NN+i,j)    =-locGless(i-1,j-1) !21
        locGmat(NN+i,NN+j) =-locGat(i-1,j-1)   !22
        !
        Smat(i,j)       = Stt(i-1,j-1)   !11
        Smat(i,NN+j)    = Sgtr(i-1,j-1)  !12
        Smat(NN+i,j)    =-Sless(i-1,j-1) !21
        Smat(NN+i,NN+j) =-Sat(i-1,j-1)   !22
     end forall

     !3)Begin inversion of Dyson equation: 
     !G = g + g*S*G ; G = g*[I + S*G] = g*\G ==> g = G * {\G}^1
     UnoMat=zero   ; forall(i=1:2*NN)UnoMat(i,i)=One/dt !Form the delta function
     GammaMat = UnoMat+matmul(Smat,locGmat)*dt          !Form \G operator
     GammaMat = GammaMat*dt**2                          !Prepare for the inversion     
     call mat_inversion_GJ(GammaMat)                       !Inversion
     G0mat    = matmul(locGmat,GammaMat)*dt             !Update G0 operators:

     !4) Extract the bigger&&lesser components
     forall(i=1:NN,j=1:NN)
        G0gtr(i-1,j-1)  = G0mat(i,j+NN)
        G0less(i-1,j-1) =-G0mat(NN+i,j)
     end forall

     forall(i=0:nstep,j=0:nstep)
        G0ret(i,j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
        g0tret(i-j)=G0ret(i,j)
     end forall
     call cfft_rt2rw(g0tret,g0fret,nstep) ; g0fret=g0fret*dt ; call swap_fftrt2rw(g0fret)

     deallocate(locGmat,Smat,G0mat,GammaMat,UnoMat)
  endif










  !=======Component by component inversion==========================
  if(TT)then
     forall(i=0:nstep,j=0:nstep)
        locGret(i,j)= heaviside(t(i)-t(j))*(locGgtr(i,j) - locGless(i,j))
        Sret(i,j)   = heaviside(t(i)-t(j))*(Sgtr(i,j) - Sless(i,j))
     end forall
     !forall(i=0:nstep)locGret(i,i)=-xi!locGless(i,i)
     locGadv=conjg(transpose(locGret))
     Sadv=conjg(transpose(Sret))

     !G0ret^-1 = Gret^-1 + Sret
     ! GammaRet(0:nstep,0:nstep) = locGret(0:nstep,0:nstep)*dt**2
     ! call mat_inversion(GammaRet(0:nstep,0:nstep))
     ! G0ret(0:nstep,0:nstep) = GammaRet(0:nstep,0:nstep) + Sret(0:nstep,0:nstep)
     ! G0ret(0:nstep,0:nstep)=G0mat(0:nstep,0:nstep)*dt**2
     ! call mat_inversion(G0ret(0:nstep,0:nstep))


     ! !G0ret = [\11 + Gret * Sret]^-1 * Gret
     Uno=zero   ; forall(i=0:nstep)Uno(i,i)=One/dt
     GammaRet(0:nstep,0:nstep) = Uno+matmul(locGret(0:nstep,0:nstep),Sret(0:nstep,0:nstep))*dt
     GammaRet(0:nstep,0:nstep) = GammaRet(0:nstep,0:nstep)*dt**2
     call mat_inversion(GammaRet(0:nstep,0:nstep))
     G0ret(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),locGret(0:nstep,0:nstep))*dt
     forall(i=0:nstep)G0ret(i,i)=-xi
     G0adv=conjg(transpose(G0ret))

     !G0less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sless * gA
     G0less(0:nstep,0:nstep) = matmul(GammaRet(0:nstep,0:nstep),matmul(locGless(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sless(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt

     !G0gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sgtr * gA
     G0gtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),matmul(locGgtr(0:nstep,0:nstep),&
          conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt  -&
          matmul(G0ret(0:nstep,0:nstep),matmul(Sgtr(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt

     !G0less(0:nstep,0:nstep)=0.8d0*dG0less(0:nstep,0:nstep)+(1.d0-0.8d0)*G0less(0:nstep,0:nstep)
     !G0gtr(0:nstep,0:nstep)=0.8d0*dG0gtr(0:nstep,0:nstep)+(1.d0-0.8d0)*G0gtr(0:nstep,0:nstep)
  endif





  !NOT WORKING and I finished my patience...:
  ! !=======Inversion of the Larkin-Ovchinnikov Matrix==========================
  ! if(FF)then
  !    !1) build Retarded/Advanced/Keldysh  GF: FF^R && FF^A=[FF^R]^+ && F^K=FF^> + FF^<
  !    forall(i=0:nstep,j=0:nstep)
  !       locGret(i,j)= heaviside(t(i-j))*(locGgtr(i,j) - locGless(i,j))
  !       Sret(i,j)   = heaviside(t(i-j))*(Sgtr(i,j) - Sless(i,j))
  !    end forall
  !    forall(i=0:nstep)locGret(i,i)=-xi
  !    forall(i=0:nstep)Sret(i,i)=-xi
  !    locGadv=conjg(transpose(locGret)) !;forall(i=0:nstep)locGadv(i,i)=zero
  !    Sadv=conjg(transpose(Sret))       !;forall(i=0:nstep)Sadv(i,i)=zero
  !    locGkel(0:nstep,0:nstep) = locGgtr(0:nstep,0:nstep) + locGless(0:nstep,0:nstep)
  !    Skel(0:nstep,0:nstep) = Sgtr(0:nstep,0:nstep) + Sless(0:nstep,0:nstep)
  !    NN=2*nstep+2                  !Size of the KS matrices
  !    !2) Build the KS matrices: \FF = {{FF^t , FF^>}, {-FF^<, -FF^at}}
  !    !G_loc(t,t')
  !    locGmat(0:nstep,0:nstep)       = locGret(0:nstep,0:nstep)
  !    locGmat(0:nstep,nstep+1:NN)    = locGkel(0:nstep,0:nstep)
  !    locGmat(nstep+1:NN,0:nstep)    = zero
  !    locGmat(nstep+1:NN,nstep+1:NN) = locGadv(0:nstep,0:nstep)!conjg(transpose(locGret(0:nstep,0:nstep)))
  !    !S(t,t')
  !    Smat(0:nstep,0:nstep)       = Sret(0:nstep,0:nstep)
  !    Smat(0:nstep,nstep+1:NN)    = Skel(0:nstep,0:nstep)
  !    Smat(nstep+1:NN,0:nstep)    = zero
  !    Smat(nstep+1:NN,nstep+1:NN) = Sadv(0:nstep,0:nstep)!conjg(transpose(Sret(0:nstep,0:nstep)))
  !    !3)Begin inversion of Dyson equation: 
  !    !G = g + g*S*G
  !    !G = g*[I + S*G] = g*\G
  !    !==> g = G * {\G}^1
  !    UnoMat=zero   ; forall(i=0:NN)UnoMat(i,i)=One/dt                           !Form the identity
  !    GammaMat(0:NN,0:NN) = UnoMat+matmul(Smat(0:NN,0:NN),locGmat(0:NN,0:NN))*dt !Form \G operator
  !    GammaMat(0:NN,0:NN) = GammaMat(0:NN,0:NN)*dt**2                            !Prepare for the inversion of continuous operator
  !    call mat_inversion(GammaMat(0:NN,0:NN))                                    !Inversion
  !    G0mat(0:NN,0:NN) = matmul(locGmat(0:NN,0:NN),GammaMat(0:NN,0:NN))*dt       !Update the KS non-interacting operators:
  !    !4) Extract the bigger&&lesser components
  !    !G0^< = G0^K + 1/2 * (-G^R + G^A) 
  !    !G0^> = G0^K + 1/2 * (G^A - G^R) 
  !    G0ret(0:nstep,0:nstep)  = G0mat(0:nstep,0:nstep)
  !    forall(i=0:nstep)G0ret(i,i)=-xi
  !    !G0ret(0:nstep,0:nstep)=zero + xi*aimag(G0ret(0:nstep,0:nstep))
  !    G0adv(0:nstep,0:nstep)=conjg(transpose(G0ret(0:nstep,0:nstep)))
  !    !forall(i=0:nstep)G0adv(i,i)=zero
  !    G0less(0:nstep,0:nstep) = G0mat(0:nstep,nstep+1:NN) + (-G0ret(0:nstep,0:nstep) + G0adv(0:nstep,0:nstep))/2.d0
  !    G0gtr(0:nstep,0:nstep)  = G0mat(0:nstep,nstep+1:NN) + (G0ret(0:nstep,0:nstep) - G0adv(0:nstep,0:nstep))/2.d0
  !    call plot_3D("dG0less3D","X","Y","Z",t(0:nstep),t(0:nstep),G0less(0:nstep,0:nstep))
  !    call plot_3D("dG0gtr3D","X","Y","Z",t(0:nstep),t(0:nstep),G0gtr(0:nstep,0:nstep))
  ! end if
