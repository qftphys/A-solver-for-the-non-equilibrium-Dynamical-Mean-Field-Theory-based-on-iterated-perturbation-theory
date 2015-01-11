  forall(i=0:nstep,j=0:nstep)
     g0tret(i-j)=heaviside(t(i-j))*(G0gtr(i,j) - G0less(i,j))
     stret(i-j)=heaviside(t(i-j))*(Sgtr(i,j) - Sless(i,j))
  end forall
  if(heaviside(0.d0)==1.d0)gtret(0)=gtret(0)/2.d0        ;    if(heaviside(0.d0)==1.d0)stret(0)=stret(0)/2.d0
  call cfft_rt2rw(g0tret,g0fret,nstep) ;    g0fret=g0fret*dt ; call swap_fftrt2rw(g0fret) !swap because F(t) are not oscillating in this formalism:
  call cfft_rt2rw(stret,sfret,nstep)  ;    sfret=dt*sfret   ; call swap_fftrt2rw(sfret)   !swap because F(t) are not oscillating in this formalism:
  gfret = one/(one/g0fret - sfret)
  do i=1,2*nstep
     w = wr(i)
     A=-aimag(gfret(i))/pi
     gfless(i)= pi2*xi*fermi(w,beta)*A
     gfgtr(i) = pi2*xi*(fermi(w,beta)-1.d0)*A
  enddo
  call cfft_rw2rt(gfless,gtless,nstep)  ; gtless=fmesh/pi2*gtless ;  gtless=gtless*exa 
  call cfft_rw2rt(gfgtr,gtgtr,nstep)   ; gtgtr =fmesh/pi2*gtgtr  ;  gtgtr=gtgtr*exa
