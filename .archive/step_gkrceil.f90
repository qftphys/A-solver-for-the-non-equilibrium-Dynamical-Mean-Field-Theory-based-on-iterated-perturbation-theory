!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine step_Gkrceil(ik,istep)
  integer :: ik,istep
  integer :: itau
  !Get the first collision integral I(0,0)
  !if(istep==1)call Ik_rceil(istep-1)
  do itau=0,Ltau
     Gkrceil(itau,istep)=Gkrceil(itau,istep-1)*conjg(Udelta(ik,istep-1))-&
          Ikrceil(itau,istep-1)*conjg(Vdelta(ik,istep-1))
  enddo
end subroutine step_Gkrceil

!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine restep_Gkrceil(ik,istep)
  integer :: ik,istep
  integer :: itau
  complex(8) :: avIkrceil
  !call Ik_rceil(istep) !build the new "layer" at istep+1    
  do itau=0,Ltau
     avIkrceil=Ikrceil(itau,istep)+oldIkrceil(itau,istep-1)
     avIkrceil=avIkrceil/2.d0
     Gkrceil(itau,istep)=Udelta(ik,istep-1)*Gkrceil(itau,istep-1)-&
          Vdelta(ik,istep-1)*avIkrceil
  enddo
end subroutine restep_Gkrceil

!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Function
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine Ik_rceil(istep)
  integer    :: i,itau,jtau,istep
  complex(8) :: I1,Ib
  do itau=0,Ltau
     I1=zero;Ib=zero
     do i=0,istep-1
        I1=I1+Gkrceil(itau,i)*Sadv(i,istep)*dt
     enddo
     if(istep==0)I1=zero
     do jtau=0,Ltau-1
        Ib=Ib+Gmktau(itau,jtau)*Srceil(jtau,istep)*dtau
     enddo
     Ikrceil(itau,istep)=I1+Ib
     !if(istep==0)Ikrceil(itau,istep)=Ib
  enddo
  return
end subroutine Ik_rceil
