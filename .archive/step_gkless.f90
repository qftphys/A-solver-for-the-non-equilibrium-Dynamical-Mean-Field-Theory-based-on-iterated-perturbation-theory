!+------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+------------------------------------------------------------+
subroutine step_Gkless(ik,istep)
  integer :: ik,istep
  integer :: it
  complex(8) :: Ik12,Ik1,Ik2
  !if(istep==1)call Ik_less(istep-1)
  do it=0,istep-1
     Gkless(it,istep)=Gkless(it,istep-1)*conjg(Udelta(ik,istep-1))- &
          Ikless(it,istep-1)*conjg(Vdelta(ik,istep-1))
     !print*,"Step",Gkless(it,istep),it,istep
  enddo
  Ik1= -conjg(Ikless(istep-1,istep-1))
  Ik2=        Ikless(istep-1,istep-1)
  Ik12 = Ik1 - Ik2   
  Gkless(istep,istep)=Gkless(istep-1,istep-1) - xi*dt*dreal(Ik12)  
  !print*,"Step",Gkless(istep,istep),it,istep
  return
end subroutine step_Gkless

!+------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+------------------------------------------------------------+
subroutine restep_Gkless(ik,istep)
  integer :: ik,istep
  integer :: it
  complex(8) :: avIkless,avIk12,Ik1,Ik2
  !call Ik_less(istep) !build the new "layer" at istep (T+\Delta)  
  do it=0,istep-1
     avIkless=Ikless(it,istep)+oldIkless(it,istep-1)
     avIkless=avIkless/2.d0
     Gkless(it,istep)=Gkless(it,istep-1)*conjg(Udelta(ik,istep-1))- &
          avIkless*conjg(Vdelta(ik,istep-1))
     !print*,"ReStep",Gkless(it,istep),it,istep
  enddo

  Ik1= -conjg(Ikless(istep,istep)) -conjg(oldIkless(istep-1,istep-1))
  Ik2= Ikless(istep,istep) + oldIkless(istep-1,istep-1)
  avIk12 = Ik1 - Ik2
  avIk12=avIk12/2.d0
  Gkless(istep,istep)=Gkless(istep-1,istep-1) - xi*dt*dreal(avIk12)
  !print*,"ReStep",Gkret(istep,istep),istep,istep
  !print*,'==================================='
  return
end subroutine restep_Gkless

!+------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Function
!PURPOSE  : 
!+------------------------------------------------------------+
subroutine Ik_less(istep)
  integer :: i,it,istep
  complex(8) :: I1,I2,Ib
  do it=0,istep
     I1=zero;I2=zero;Ib=zero
     do i=0,it-1
        I1=I1+Gkret(it,i)*Sless(i,istep)*dt
     enddo
     if(it==0)I1=zero
     do i=0,istep-1
        I2=I2+Gkless(it,i)*Sadv(i,istep)*dt
     enddo
     if(istep==0)I2=zero
     do i=0,Ltau
        Ib=Ib+Gklceil(it,i)*Srceil(i,istep)*dtau
     enddo
     !print*,"Ir,t,t`",I1+I2-xi*Ib,it,istep
     Ikless(it,istep)=I1+I2-xi*Ib
  enddo
  !if(istep==0)Ikless(0,0)=-xi*Ib!zero
  return
end subroutine Ik_less

