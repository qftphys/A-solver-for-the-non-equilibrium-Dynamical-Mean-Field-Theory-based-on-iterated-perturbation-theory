  Uno=0.d0 
  forall(i=0:nstep)Uno(i,i)=1.d0/dt

  !Get G_loc^{A,R}, S^{A,R}
  dGret=zero ; dSret=zero
  do i=0,nstep
     do j=0,nstep
        dGret(i,j)=heaviside(t(i-j))*(locGgtr(i,j) - locGless(i,j))
        dSret(i,j)=heaviside(t(i-j))*(Sgtr(i,j) - Sless(i,j))
     enddo
  enddo
  dSadv=conjg(transpose(dSret))
  dGadv=conjg(transpose(dGret))


  GammaRet = Uno+matmul(dGret(0:nstep,0:nstep),dSret(0:nstep,0:nstep))*dt
  GammaAdv = Uno+matmul(dSadv(0:nstep,0:nstep),dGadv(0:nstep,0:nstep))*dt

  forall(i=0:nstep,j=0:nstep)
     dgtret(i-j)=dGret(i,j)
     dstret(i-j)=dSret(i,j)
     dgtadv(i-j)=dGadv(i,j)
     dstadv(i-j)=dSadv(i,j)
     dggtret(i-j)=GammaRet(i,j)
     dggtret(i-j)=GammaAdv(i,j)
  end forall
  call splot("provaGret_t.ipt",t(-nstep:nstep),dgtret,append=TT)
  call splot("provaSret_t.ipt",t(-nstep:nstep),dstret,append=TT)
  call splot("provaGadv_t.ipt",t(-nstep:nstep),dgtadv,append=TT)
  call splot("provaSadv_t.ipt",t(-nstep:nstep),dstadv,append=TT)
  call splot("provaGammaRet_t.ipt",t(-nstep:nstep),dggtret,append=TT)
  call splot("provaGammaAdv_t.ipt",t(-nstep:nstep),dggtadv,append=TT)


  call mat_inversion(GammaRet(0:nstep,0:nstep))
  call mat_inversion(GammaAdv(0:nstep,0:nstep))

  dG0ret = matmul(GammaRet,dGret)*dt
  dG0adv = matmul(dGadv,GammaAdv)*dt

  !Update G_0^<, G0^>
  dG0less=matmul(GammaRet,matmul(dGless,GammaAdv))*dt**2 - matmul(dG0ret,matmul(dSless,dG0adv))*dt**2
  dG0gtr =matmul(GammaRet,matmul(dGgtr,GammaAdv))*dt**2  - matmul(dG0ret,matmul(dSgtr,dG0adv))*dt**2

  forall(i=0:nstep,j=0:nstep)
     dg0tless(i-j)=dG0less(i,j)
     dg0tgtr(i-j)=dG0gtr(i,j)
     dg0tret(i-j)=dG0ret(i,j)
  end forall
  call splot("provaG0less_t.ipt",t(-nstep:nstep),dg0tless,append=TT)
  call splot("provaG0gtr_t.ipt",t(-nstep:nstep),dg0tgtr,append=TT)  
  call splot("provaG0ret_t1.ipt",t(-nstep:nstep),dg0tret,append=TT)
  do i=-nstep,nstep
     dg0tret(i)=heaviside(t(i))*(dg0tgtr(i) - dg0tless(i))
  enddo
  call splot("provaG0ret_t2.ipt",t(-nstep:nstep),dg0tret,append=TT)
