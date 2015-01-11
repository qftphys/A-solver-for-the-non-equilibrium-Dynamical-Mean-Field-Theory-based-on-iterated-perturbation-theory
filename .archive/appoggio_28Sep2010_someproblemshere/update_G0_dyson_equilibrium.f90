  do i=0,nstep
     do j=0,nstep
        gtret(i-j) = heaviside(t(i-j))*(locGgtr(i,j)-locGless(i,j))
        stret(i-j) = heaviside(t(i-j))*(Sgtr(i,j)-Sless(i,j))
     enddo
  enddo
  if(heaviside(0.d0)==1.d0)gtret(0)=gtret(0)/2.d0
  if(heaviside(0.d0)==1.d0)stret(0)=stret(0)/2.d0
  call cfft_rt2rw(gtret,gfret,nstep) ; gfret=gfret*dt ; call swap_fftrt2rw(gfret)
  call cfft_rt2rw(stret,sfret,nstep) ; sfret=dt*sfret ; call swap_fftrt2rw(sfret)
  g0fret=one/(one/gfret + sfret)
  do i=1,2*nstep
     w = wr(i)
     A = -aimag(g0fret(i))/pi
     An= A*fermi(w,beta)
     g0fless(i)= pi2*xi*An
     g0fgtr(i) = pi2*xi*(An-A)
  enddo
  call cfft_rw2rt(g0fless,g0tless,nstep) ; g0tless=fmesh/pi2*g0tless  ; g0tless=g0tless*exa 
  call cfft_rw2rt(g0fgtr, g0tgtr,nstep)  ; g0tgtr =fmesh/pi2*g0tgtr   ; g0tgtr=g0tgtr*exa
