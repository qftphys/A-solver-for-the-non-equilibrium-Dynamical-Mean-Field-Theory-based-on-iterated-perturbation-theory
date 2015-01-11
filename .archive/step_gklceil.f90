!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine step_Gklceil(ik,istep)
  integer :: ik,istep
  integer :: itau
  !Get the first collision integral I(0,0)
  !if(istep==1)call Ik_lceil(istep-1)
  do itau=0,Ltau
     Gklceil(istep,itau)=Udelta(ik,istep-1)*Gklceil(istep-1,itau)-&
          Vdelta(ik,istep-1)*Iklceil(istep-1,itau)
  enddo
end subroutine step_Gklceil

!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine restep_Gklceil(ik,istep)
  integer :: ik,istep
  integer :: itau
  complex(8) :: avIklceil
  !call Ik_lceil(istep) !build the new "layer" at istep+1    
  do itau=0,Ltau
     avIklceil=Iklceil(istep,itau)+oldIklceil(istep-1,itau)
     avIklceil=avIklceil/2.d0
     Gklceil(istep,itau)=Udelta(ik,istep-1)*Gklceil(istep-1,itau)-&
          Vdelta(ik,istep-1)*avIklceil
  enddo
end subroutine restep_Gklceil

!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Function
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine Ik_lceil(istep)
  integer    :: i,itau,jtau,istep
  complex(8) :: I1,Ib
  do itau=0,Ltau
     I1=zero;Ib=zero
     do i=0,istep-1
        I1=I1+Sret(istep,i)*Gklceil(i,itau)*dt
     enddo
     if(istep==0)I1=zero
     do jtau=0,Ltau-1
        Ib=Ib+Slceil(istep,jtau)*Gmktau(jtau,itau)*dtau
     enddo
     Iklceil(istep,itau)=I1+Ib
     !if(istep==0)Iklceil(istep,itau)=Ib
  enddo

  return
end subroutine Ik_lceil
