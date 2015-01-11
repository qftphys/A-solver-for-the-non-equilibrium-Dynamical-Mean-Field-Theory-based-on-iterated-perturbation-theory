  !Set up the delta-function-matrix
  Uno=zero  ; forall(i=0:nstep)Uno(i,i)=One/dt

  !Store the S0 TTI to  matrices X0:
  forall(i=0:nstep,j=0:nstep)
     X0less(i,j)=S0less(i-j)
     X0gtr(i,j) =S0gtr(i-j)
  end forall

  !Get the retarded/advanced G0/X component
  forall(i=0:nstep,j=0:nstep)
     X0ret(i,j)= heaviside(t(i)-t(j))*(X0gtr(i,j) - X0less(i,j))
     G0ret(i,j)= heaviside(t(i)-t(j))*(G0gtr(i,j) - G0less(i,j))
  end forall
  G0adv=conjg(transpose(G0ret))

  !Get the Gret:
  !GammaM= [\11 - G0ret * S0ret]
  !Gret  = GammaM^-1 * G0ret
  GammaM(0:nstep,0:nstep) = Uno-matmul(G0ret(0:nstep,0:nstep),X0ret(0:nstep,0:nstep))*dt
  GammaM(0:nstep,0:nstep) = GammaM(0:nstep,0:nstep)*dt**2
  call mat_inversion(GammaM(0:nstep,0:nstep))
  Gret(0:nstep,0:nstep) = matmul(GammaM(0:nstep,0:nstep),G0ret(0:nstep,0:nstep))*dt
  forall(i=0:nstep)Gret(i,i)=-xi
  Gadv=conjg(transpose(Gret))

  !Get the GammaP (used in the following):
  GammaP(0:nstep,0:nstep) = Uno+matmul(Gret(0:nstep,0:nstep),X0ret(0:nstep,0:nstep))*dt

  !Gless = GammaP * G0less * GammaP^+  +  GRet * S0less * Gret^+
  Gless(0:nstep,0:nstep) = matmul(GammaP(0:nstep,0:nstep),matmul(G0less(0:nstep,0:nstep),&
       conjg(transpose(GammaP(0:nstep,0:nstep))))*dt)*dt + &
       matmul(Gret(0:nstep,0:nstep),matmul(X0less(0:nstep,0:nstep),Gadv(0:nstep,0:nstep))*dt)*dt

  !Ggtr  = GammaP * G0gtr * GammaP^+   +  Gret * S0gtr * Gret^+
  Ggtr(0:nstep,0:nstep)  = matmul(GammaP(0:nstep,0:nstep),matmul(G0gtr(0:nstep,0:nstep),&
       conjg(transpose(GammaP(0:nstep,0:nstep))))*dt)*dt  +&
       matmul(Gret(0:nstep,0:nstep),matmul(X0gtr(0:nstep,0:nstep),Gadv(0:nstep,0:nstep))*dt)*dt

  G0less=Gless
  G0gtr =Ggtr
