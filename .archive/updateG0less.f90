  !Define auxiliary arrays
  dummy_locGret(0:nstep,0:nstep) =locGret(0:nstep,0:nstep)
  dummy_locGadv(0:nstep,0:nstep) =conjg(transpose(locGret(0:nstep,0:nstep)))
  dummy_locGless(0:nstep,0:nstep)=locGless(0:nstep,0:nstep)
  dummy_locGlceil(0:nstep,0:Ltau)=locGlceil(0:nstep,0:Ltau)
  dummy_locGrceil(0:Ltau,0:nstep)=conjg(transpose(locGlceil(0:nstep,0:Ltau)))
  G0adv=conjg(transpose(G0ret(0:nstep,0:nstep)))

  Unity=zero
  do i=0,nstep
     Unity(i,i)=one
  enddo

  !Build the inverse of the operators [1+G^r.S^r], [1+S^a.G^a]
  OBJret=Unity+matmul(dummy_locGret,Sret(0:nstep,0:nstep))!*dt
  call InvMat(OBJret(0:nstep,0:nstep),nstep+1) 

  OBJadv=Unity+matmul(Sadv(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))!*dt
  call InvMat(OBJadv(0:nstep,0:nstep),nstep+1) 

  !Op1 = Or^-1*G^<*Oa^-1
  Op1=matmul(OBJret,matmul(dummy_locGless,OBJadv)) !*dt !*dt

  !Op2 = Or^-1*[G^r.S^<.G^a]*Oa^-1
  Op2=matmul(dummy_locGret,matmul(Sless(0:nstep,0:nstep),dummy_locGadv))*dt*dt
  Op2=matmul(OBJret,matmul(Op2,OBJadv)) !*dt !*dt

  !Op3 = Or^-1*[G^lceil.S^rceil.G^a]*Oa^-1
  Op3=matmul(&
       matmul(dummy_locGlceil,Srceil(0:nstep,0:Ltau)),dummy_locGadv)*dt*dtau
  Op3=matmul(OBJret,matmul(Op3,OBJadv)) !*dt !*dt

  !Op4 = Or^-1*[G^r.S^lceil.G0^rceil]
  Op4=matmul(dummy_locGret,&
       matmul(Slceil(0:nstep,0:Ltau),G0rceil(0:Ltau,0:nstep)))*dt*dtau
  Op4=matmul(OBJret,Op4) !*dt

  !Op5 = Or^-1*[G^lceil.S^M.G0^rceil]  
  Op5=matmul(dummy_locGlceil,&
       matmul(Smatsubara(0:Ltau,0:Ltau),G0rceil(0:Ltau,0:nstep)))*dtau*dtau
  Op5=matmul(OBJret,Op5) !*dt

  G0less(0:nstep,0:nstep)=Op1-Op2-Op3-Op4!-Op5
  G0adv(0:nstep,0:nstep)=-Op5
  call plot_dislin("locG0lessOp34_t1t2","X","Y","Z",&
       t(0:nstep),t(0:nstep),G0adv(0:nstep,0:nstep))
